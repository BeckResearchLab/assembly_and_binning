# Import matplotlib before seaborn
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import matplotlib

import pandas as pd
import seaborn as sns

def load_data():
    """
    load all the data, as parsed by samtools idxstats.
    """
    df_all = pd.read_csv('./results/individual_contig_summary.tsv', sep='\t')
    df_all['frac reads for contig'] = df_all['# mapped reads']/df_all['total reads (in fastq)']
    return df_all 

def load_only_binned_or_unbinned(binned=False, groupby_sample=False):
    df = load_data()
    df = df[df['binned contig'] == binned]
    if groupby_sample:
        grouped = pd.DataFrame(df.groupby('sample id')['frac reads for contig'].sum())
        grouped.reset_index(inplace=True)
        grouped.rename(columns={'frac reads for contig':'sum(frac reads for contig)'}, inplace=True)
        return grouped
    else: 
        return df

def get_sample_info():
    df_all = load_data()
    return df_all[['sample id', 'oxygen', 'replicate', 'week']].drop_duplicates()

def unbinned_contigs_appearing_with_at_least_x_percent_of_reads_in_one_sample(x):
    df = load_data()
    # remove binned contigs 
    df = df[df['binned contig'] == False]
    # contig names that have appeared with x percent in at least one sample .
    frac = x/100.
    contigs = df[df['frac reads for contig'] > frac].contig.unique()
    print('number of contigs for importance level {}%: {}'.format(x, len(contigs)))
    return df[df['contig'].isin(contigs)]

def count_contigs_by_importance(percent_list):
    for p in percent_list:
        df = unbinned_contigs_appearing_with_at_least_x_percent_of_reads_in_one_sample(p)
        contigs = len(df['contig'].unique()) 
    print(p, contigs)
    return None

def sort_contigs_by_importance(percent_cutoff, filename=None):
    df = unbinned_contigs_appearing_with_at_least_x_percent_of_reads_in_one_sample(percent_cutoff)
    max_importance = pd.DataFrame(df.groupby(['contig', 'contig length'])['frac reads for contig'].max())
    max_importance.reset_index(inplace=True)
    max_importance.rename(columns={'frac reads for contig':'max(frac reads for contig)'}, inplace=True)
    max_importance.sort_values('max(frac reads for contig)', ascending=False, inplace=True)
    if filename is not None:
        max_importance.to_csv(filename, sep='\t', index=False)
    else:
        max_importance.to_csv('./results/most_important_contigs_for_cutoff_{}_percent.tsv'.format(percent_cutoff), 
            sep='\t', index=False)
    return max_importance

def df_to_assess_impact_of_adding_contigs_to_bins(percent_cutoff):
    # get stats for binned contigs:
    bc = load_only_binned_or_unbinned(binned=True, groupby_sample=True)
    bc.rename(columns={'sum(frac reads for contig)':'frac reads accounted for by binned contigs'}, inplace=True)
    print(bc.head())
    
    # get stats for un-binned contigs:
    ubc = unbinned_contigs_appearing_with_at_least_x_percent_of_reads_in_one_sample(percent_cutoff) 
    ubc = pd.DataFrame(ubc.groupby('sample id')['frac reads for contig'].sum())
    ubc.reset_index(inplace=True)
    ubc.rename(columns={'frac reads for contig':'frac reads accounted for by top un-binned contigs'}, inplace=True)
    print(ubc.head())
    
    # merge them
    df_improvement_potential = pd.merge(bc, ubc)
    si = get_sample_info()
    df_improvement_potential = pd.merge(df_improvement_potential, si)
    df_improvement_potential.head()
    assert df_improvement_potential.shape[0] == si.shape[0], 'expected {} rows'.format(si.shape[0])

    return df_improvement_potential

def plot_assess_impact_of_adding_contigs_to_bins(percent_cutoff):
    df = df_to_assess_impact_of_adding_contigs_to_bins(percent_cutoff).copy()
    df['sum'] = df['frac reads accounted for by binned contigs'] + df['frac reads accounted for by top un-binned contigs']
    fig, axs = plt.subplots(2, 1, figsize=(7, 5))
    o2_dict = {'low': axs[0], 'high': axs[1]}
    colors = {1:'#66c2a5', 2:'#fc8d62', 3:'#8da0cb', 4:'#e78ac3'}
    for (o2, rep), plot_df in df.groupby(['oxygen', 'replicate']):
        plot_df.sort_values('week', inplace=True)
        ax = o2_dict[o2]
        color = colors[rep]
        ax.plot(plot_df['week'], plot_df['frac reads accounted for by binned contigs'], 
                linestyle='-', marker='o', color=color, alpha = 0.3)
        ax.plot(plot_df['week'], plot_df['sum'], linestyle='-', marker='o', color=color)

    labels = ['rep {}'.format(n) for n in [1, 2, 3, 4]]
    for a in axs:
        ax.set_xlabel('week')
        ax.set_ylabel('fraction of reads\naccounted for')
        ax.legend(labels, bbox_to_anchor=(1.15, 1.05))

    axs[0].set_title('low oxygen')

    axs[1].set_title('high oxygen')
    axs[1].set_ylabel('fraction of reads\naccounted for')
    plt.tight_layout()

    return fig

def determine_max_length(df_list, colname):
    # get max length for plotting purposes
    maxs = list()
    for df in df_list:
        maxs.append(df[colname].max())
    max_len = max(maxs)
    return max_len

def plot_distributions_of_contig_lengths(plot_col, logy=True):
    all_contigs = load_data()
    binned = load_only_binned_or_unbinned(binned=True)
    unbinned = load_only_binned_or_unbinned(binned=False)
    dfs = [all_contigs, binned, unbinned]
    
    ml = determine_max_length(dfs, plot_col)
    titles = ['all contigs', 'binned contigs', 'unbinned contigs']

    binwidth = 2000
    max_bin_size = max([df[plot_col].max() for df in dfs])
    bins = bins=np.arange(0, ml + binwidth, binwidth)

    fig, axs = plt.subplots(3, 1, figsize=(7, 5), sharex=True, sharey=True)    
    for (df, axnum) in zip(dfs, [0, 1, 2]):
        ax = axs[axnum]
        ax.hist(df[plot_col], bins = bins)
        ax.set_title(titles[axnum])
        if logy:
            ax.set_yscale('log')
        ax.set_xlabel('contig size (bp)')
        ax.set_ylabel(plot_col)
        ax.get_xaxis().set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
    fig.tight_layout()
    return fig

def plot_distributions_of_contig_lengths_one_plot(plot_col, logy=True):
    """
    Plot a black bar for the total number of contigs in each 2kb range, 
    and a green bar for the subset of those that are binned.
    The log scale makes it hard to see the missing bins, so see other plots.
    """
    all_contigs = load_data()
    binned = load_only_binned_or_unbinned(binned=True)
    dfs = [all_contigs, binned]
    
    ml = determine_max_length(dfs, plot_col)

    binwidth = 2000
    max_bin_size = max([df[plot_col].max() for df in dfs])
    bins = bins=np.arange(0, ml + binwidth, binwidth)

    fig, ax = plt.subplots(1, 1, figsize=(8, 3.5))    
    ax.hist(binned[plot_col], color='#000000', bins=bins)
    ax.hist(binned[plot_col], color='#31a354', bins=bins) # green
    ax.set_title('portion of contigs binned at different contig lengths')
    if logy:
        ax.set_yscale('log')
    ax.set_xlabel('contig size (bp)')
    plt.axvline(x=1000, color='k', linestyle='--')
    ax.set_ylabel('# contigs')
    ax.get_xaxis().set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
    fig.tight_layout()
    return fig


def count_contigs_by_length_range(df, bin_width):
    """
    Count number of contigs in intervals specified by bin_width.
    """
    max_length = determine_max_length([df], 'contig length')
    bin_edges = np.arange(0, max_length + bin_width, bin_width, dtype=int).tolist()
    print('bin edges (sample): {}'.format(bin_edges[0:5]))
    
    df['hist bin'] = pd.cut(df['contig length'], bins=bin_edges)
    num_contigs = pd.DataFrame(df.groupby('hist bin')['contig'].count())
    num_contigs.reset_index(inplace=True)
    num_contigs.rename(columns={'contig':'# contigs'}, inplace=True)
    num_contigs['# contigs'].fillna(0)
    return num_contigs 

def plot_frac_reads_binned_at_different_contig_lengths(bin_width, logx=True):
    """
    While plotting total numbers of contigs at different bin lengths allows
    observation of the vastly higher number of shorter contigs, the fraction
    binned at different lengths is hard to see. 
    Plotting as fractions hides the total number falling in each interval but
    reveals what length contigs failed to be binned. 
    """
    all_contigs = load_data()
    binned = load_only_binned_or_unbinned(binned=True)
    unbinned = load_only_binned_or_unbinned(binned=False)
    dfs = [all_contigs, binned, unbinned]
    acc = count_contigs_by_length_range(df=all_contigs, bin_width=bin_width)
    bc = count_contigs_by_length_range(df=binned, bin_width=bin_width)
    ubc = count_contigs_by_length_range(df=unbinned, bin_width=bin_width)
    acc.rename(columns={'# contigs':'# binned + unbinned contigs'}, inplace=True)
    bc.rename(columns={'# contigs':'# binned contigs'}, inplace=True)
    ubc.rename(columns={'# contigs':'# unbinned contigs'}, inplace=True)
    fracs = pd.merge(bc, ubc, how='outer')
    fracs = pd.merge(fracs, acc, how='outer')

    fracs['# unbinned contigs'].fillna(0, inplace=True)
    fracs['# binned contigs'].fillna(0, inplace=True)
    fracs['sum of contigs'] = fracs['# binned contigs'] + fracs['# unbinned contigs']
    fracs['sum of contigs'].fillna(0, inplace=True)
    fracs['frac binned contigs'] = fracs['# binned contigs']/fracs['# binned + unbinned contigs']
    fracs['frac unbinned contigs'] = fracs['# unbinned contigs']/fracs['# binned + unbinned contigs']
    # indepenent check:
    assert (fracs['# binned + unbinned contigs'] - fracs['sum of contigs']).fillna(0).max() < 0.000001
    
    fracs['upper bound for contig length'] = bc['hist bin'].str.extract('([0-9]+)]')
    xlabels = fracs['upper bound for contig length']
    frac_binned = fracs['frac binned contigs']
    frac_unbinned = fracs['frac unbinned contigs']

    fig, ax = plt.subplots(1, 1, figsize=(7, 3), sharex=True, sharey=True)    
    # scatter, not bar, because some of the fractions are NaN b/c of dividing by zero. 
    ax.scatter([int(i) for i in xlabels.tolist()], frac_binned.tolist(), alpha=0.5)
    plt.axvline(x=1000, color='k', linestyle='--')
    ax.get_xaxis().set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
    ax.set_xlabel('approx. contig size')
    ax.set_ylabel('fraction of contigs binned\n({} bp resolution)'.format(bin_width))
    if logx:
        ax.set_xscale('log')

    return fig
