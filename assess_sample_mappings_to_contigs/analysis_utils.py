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

def frac_reads_binned_at_different_contig_lengths(bin_width):
    """
    Prepare a dataframe with statistics about the fractions binned and unbinned
    along with the count of contigs at each contig-length interval.
    Igores # of reads mapped at different lengths.
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
    return fracs

def plot_frac_reads_binned_at_different_contig_lengths(bin_width, logx=True):
    """
    While plotting total numbers of contigs at different bin lengths allows
    observation of the vastly higher number of shorter contigs, the fraction
    binned at different lengths is hard to see. 
    Plotting as fractions hides the total number falling in each interval but
    reveals what length contigs failed to be binned. 
    """
    fracs = frac_reads_binned_at_different_contig_lengths(bin_width)
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

def plot_frac_reads_binned_at_different_contig_lengths_and_total(bin_width, logx=True):
    """
    """
    fracs = frac_reads_binned_at_different_contig_lengths(bin_width)
    xlabels = fracs['upper bound for contig length']
    frac_binned = fracs['frac binned contigs']
    frac_unbinned = fracs['frac unbinned contigs']
    x = [int(i) for i in xlabels.tolist()]

    fig, axs = plt.subplots(2, 1, figsize=(7, 6), sharex=True, sharey=False)    
    # scatter, not bar, because some of the fractions are NaN b/c of dividing by zero. 
    # Top plot is total number of contigs.  Log scale. 
    ax = axs[0]
    y = fracs['# binned + unbinned contigs'].tolist()
    y =  [np.nan if yval == 0 else yval for yval in y]
    ax.scatter(x, y, alpha=0.5, color='#636363')
    ax.set_yscale('log')
    ax.set_ylim(1, max(y)) # cropping 10^7 to 10^8 without this
    ax.set_ylabel('total number of contigs\n({} bp resolution)'.format(bin_width))

    # Bottom plot is frac mapped. 
    ax = axs[1]
    ax.scatter(x, frac_binned.tolist(), alpha=0.5, color='#3182bd')
    ax.set_xlabel('approx. contig size')
    ax.set_ylabel('fraction of contigs binned\n({} bp resolution)'.format(bin_width))

    for ax in axs:
        if logx:
            ax.set_xscale('log')
        ax.get_xaxis().set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))

    return fig

def summarise_contigs_by_contgig_size(
    size_boundaries_list=[10**2, 500, 10**3, 5*10**3, 10**4, 5*10**4, 10**5, 5*10**5, 10**6]):
    df = load_data()
    df['hist bin'] = pd.cut(df['contig length'], bins=size_boundaries_list)
    df['(# mapped reads)/(contig len)'] = df['# mapped reads']/df['contig length']
    df_orig = df.copy()
    
    # add count of reads by sample across contig length
    total_mapped_by_sample = df.groupby('sample id')['# mapped reads'].sum().to_frame().reset_index()
    total_mapped_by_sample.rename(columns={'# mapped reads':'sample sum(reads mapped)'}, 
                                  inplace=True)
    df = pd.merge(df, total_mapped_by_sample, how='outer')
    
    # add count for frac mapped reads
    binned_reads_by_sample = df[df['binned contig'] == True].groupby('sample id')['# mapped reads'].sum().to_frame().reset_index()
    binned_reads_by_sample.rename(columns={'# mapped reads':
                                           'sample sum(reads mapped) (all bins)'}, 
                                  inplace=True)
    df = pd.merge(df, binned_reads_by_sample, how='outer')
    
    num_reads = df.groupby(['sample id','hist bin'])['# mapped reads'].sum().to_frame().reset_index()
    result = df[['sample id', 'oxygen', 'replicate', 'week', 
                 'total reads (in fastq)', 'sample sum(reads mapped)',
                 'sample sum(reads mapped) (all bins)']
               ].drop_duplicates()
    result = pd.merge(result, num_reads)
    result['upper bound for contig length'] = result['hist bin'].str.extract('([0-9]+)]').astype(int)
    
    # add avg counts/bp:
    counts_per_bp = df_orig.groupby(['sample id','hist bin'])['(# mapped reads)/(contig len)'].mean().to_frame().reset_index()
    counts_per_bp.rename(columns={'(# mapped reads)/(contig len)':
                                  'mean((# mapped reads)/(contig len))'}, 
                                  inplace=True)
    result = pd.merge(result, counts_per_bp, how='outer')
    # add avg counts/bp, unbinned :
    counts_per_bp_unbinned = df_orig[df_orig['binned contig'] == False].groupby(
        ['sample id','hist bin'])['(# mapped reads)/(contig len)'
                                 ].mean().to_frame().reset_index()
    counts_per_bp_unbinned.rename(columns={'(# mapped reads)/(contig len)':
                                  'mean((# mapped reads)/(contig len)), unbinned'}, 
                                  inplace=True)
    result = pd.merge(result, counts_per_bp_unbinned, how='outer')
    
    # add avg counts/bp, binned :
    counts_per_bp_unbinned = df_orig[df_orig['binned contig'] == True].groupby(
        ['sample id','hist bin'])['(# mapped reads)/(contig len)'
                                 ].mean().to_frame().reset_index()
    counts_per_bp_unbinned.rename(columns={'(# mapped reads)/(contig len)':
                                  'mean((# mapped reads)/(contig len)), binned'}, 
                                  inplace=True)
    result = pd.merge(result, counts_per_bp_unbinned, how='outer')
    
    return result
        
def plot_num_reads_assigned_to_contigs_shorter_than_length_x(x=1500):
    df = summarise_contigs_by_contgig_size()
    print(df.head(1))
    fig, axs = plt.subplots(2, 1, figsize=(7, 5), sharex=True, sharey=True, 
                            subplot_kw = {'ylim':(0,1)})
    o2_dict = {'low': axs[0], 'high': axs[1]}
    colors = {1:'#66c2a5', 2:'#fc8d62', 3:'#8da0cb', 4:'#e78ac3'}
    
    max_contig_size = x
    dfp = df[df['upper bound for contig length'] <= max_contig_size]
    sample_sums = dfp.groupby('sample id')['# mapped reads'].sum().to_frame().reset_index()
    sample_sums.rename(columns={'# mapped reads':'sum(mapped reads) for contigs < {}bp'.format(max_contig_size)}, inplace=True)

    dfp = pd.merge(dfp, sample_sums)
    dfp['frac reads mapped for contigs < {}bp'.format(max_contig_size)] = dfp['sum(mapped reads) for contigs < {}bp'.format(max_contig_size)]/dfp['total reads (in fastq)']
    dfp['sample total reads mapped'] = dfp['sample sum(reads mapped)']/dfp['total reads (in fastq)']
    dfp['sample total reads mapped to bins'] = dfp['sample sum(reads mapped) (all bins)']/dfp['total reads (in fastq)']
    dfp['max contig size'] = max_contig_size
    
    
    fig.suptitle('number of reads assigned to contigs shorter than {} bp'.format(x))
    labels = ['rep {}'.format(n) for n in [1, 2, 3, 4]]
    
    for (o2, rep), plot_df in dfp.groupby(['oxygen', 'replicate']):
        print(o2, rep)
        plot_df.sort_values('week', inplace=True)
        ax = o2_dict[o2]
        color = colors[rep]

        # frac of reads mapped to contigs in df
        ax.plot(plot_df['week'], plot_df['frac reads mapped for contigs < {}bp'.format(max_contig_size)], 
                linestyle='-', marker="H", color=color, label="replicate {} (reads for contigs < {}bp)".format(rep, max_contig_size))
        handles, labels = ax.get_legend_handles_labels()  # "$\u2639$"
        lgd = ax.legend(handles, labels, bbox_to_anchor=(1.65, 1.05))
        
        # all the reads mapped to contigs
        ax.plot(plot_df['week'], plot_df['sample total reads mapped'],
                linestyle='--', marker="o", color=color, alpha = .3, label="")
        
        # frac of reads mapped to binned contigs
        ax.plot(plot_df['week'], plot_df['sample total reads mapped to bins'],
                linestyle='--', marker=r"$b$", color=color, alpha = .3, label="")
        
        ax.set_ylabel('fraction of reads'.format(x))       
        ax.set_xlabel('week')
        
    axs[0].set_title('low oxygen')
    axs[1].set_title('high oxygen')
    plt.subplots_adjust(hspace=0.25)
    
    return fig

def plot_log_ratio_of_coverage_like_metric():
    df = summarise_contigs_by_contgig_size()
    fig, axs = plt.subplots(2, 1, figsize=(4, 8), sharex=True, sharey=True, )
    o2_dict = {'low': axs[0], 'high': axs[1]}
    colors = {1:'#66c2a5', 2:'#fc8d62', 3:'#8da0cb', 4:'#e78ac3'}
    
    x = 'mean((# mapped reads)/(contig len)), binned'
    y = 'mean((# mapped reads)/(contig len)), unbinned'
    df['ratio'] = df[y]/df[x]
    df['log_2(ratio)'] = np.log2(df[y]/df[x])
    min_ratio = df['ratio'].min()
    max_ratio = df['ratio'].max()
    print('min: {}, max: {}'.format(min_ratio, max_ratio))
    
    fig.suptitle('log ratio of mean "coverage"')
    labels = ['rep {}'.format(n) for n in [1, 2, 3, 4]]
    
    for (o2, rep), plot_df in df.groupby(['oxygen', 'replicate']):
        #plot_df.sort_values('week', inplace=True)
        ax = o2_dict[o2]
        color = colors[rep]

        # frac of reads mapped to contigs in df

        print((plot_df[y]/plot_df[x]).head(10))
        
        ax.plot(plot_df['upper bound for contig length'], 
                plot_df['log_2(ratio)'],
                linestyle='', marker="o", color=color, label="replicate {}".format(rep), 
                alpha = 0.3)
        handles, labels = ax.get_legend_handles_labels()  
        lgd = ax.legend(handles, labels, bbox_to_anchor=(1.5, 1.05))
    
        ax.set_xlabel('contig length-ish')
        ax.set_ylabel('log2(unbinned/binned) for reads/len')       
        ax.axhline(y=0, color='gray', linestyle='-')
        ax.set_xscale('log')
        
    axs[0].set_title('low oxygen')
    axs[1].set_title('high oxygen')
    plt.subplots_adjust(hspace=0.25)
    
    return fig
    
def plot_num_reads_across_samples():
    """
    Plot number of reads in .fastq (unproccessed original data)
    """
    df = load_data()
    info = df[['sample id', 'oxygen', 'replicate', 'week', 
               'total reads (in fastq)']].drop_duplicates()
    print("info.shape: {}".format(info.shape))
    print(info.head())
    fig, axs = plt.subplots(2, 1, figsize=(4, 8), sharex=True, sharey=True, )
    
    o2_dict = {'low': axs[0], 'high': axs[1]}
    colors = {1:'#66c2a5', 2:'#fc8d62', 3:'#8da0cb', 4:'#e78ac3'}

    x = 'week'
    y = 'total reads (in fastq)'
    
    fig.suptitle("# reads in each samples' .fastq", size=16)
    
    for (o2, rep), plot_df in info.groupby(['oxygen', 'replicate']):
        plot_df.sort_values('week', inplace=True)
        ax = o2_dict[o2]
        color = colors[rep]
        
        ax.plot(plot_df[x], plot_df[y],
                linestyle='-', marker="o", color=color, label="replicate {}".format(rep),
                alpha = 1)
        handles, labels = ax.get_legend_handles_labels()
        lgd = ax.legend(handles, labels, bbox_to_anchor=(1.5, 1.05))

        ax.set_xlabel(x)
        ax.set_ylabel(y)
        ax.axhline(y=0, color='gray', linestyle='-')

    axs[0].set_title('low oxygen')
    axs[1].set_title('high oxygen')
    plt.subplots_adjust(hspace=0.25)
    
    return fig
    
def prepare_bins(binwidth):
    # got max_size from max_contig_size = load_data()['contig length'].max()
    max_contig_size=686225
    print('max contig size: {}'.format(max_contig_size))
    bins = np.arange(0, max_contig_size + binwidth, binwidth).tolist()
    print('first bins for binwidth {}: {}'.format(binwidth, bins[0:10]))
    return bins

def plot_good_vs_bad_low_o2_samples(binwidth):
    bins = prepare_bins(binwidth)
    df = summarise_contigs_by_contgig_size(bins)
    
    df['frac reads assigned to contigs this length'] = df['# mapped reads']/df['total reads (in fastq)']
    print(df.columns)
    df_extract2 = df[(df['oxygen']== 'low') & 
                     ((df['replicate']== 2) | (df['replicate']== 1)) & 
                     df['week'].isin([8, 9, 10, 11, 12, 13])]
    df_extracts = df_extract2
    
    
    fig, axs = plt.subplots(1, 1, figsize=(4, 3)) #, sharex=True, sharey=True)
    
    green = '#31a354'
    orange = '#d95f0e'
    colors = {('low', 8):orange, ('low', 9):orange, ('low', 10):orange, ('low', 11):green, ('low', 12):green, ('low', 13):green,
              ('high', 10):green, ('high', 11):orange}
    print(colors)
    
    shapes = {1:'o', 2:'^'}

    x = 'upper bound for contig length'
    y = 'frac reads assigned to contigs this length'
    
    fig.suptitle("Fraction of reads assigned to contigs of different lengths", size=14)
    
    def get_unique(df_x, col):
        uniques = df_x[col].unique()
        assert len(uniques) == 1, "too many unique values for {}: {}".format(col, uniques)
        return uniques[0]
    
    df_extracts.sort_values('week', inplace=True)
    
    for sample_id, plot_df in df_extracts.groupby(['sample id'], sort=False):
        o2 = get_unique(plot_df, 'oxygen') 
        rep = get_unique(plot_df, 'replicate') 
        week = get_unique(plot_df, 'week')
        label = "{} O2, rep {}, week {}".format(o2, rep, week)
        print(label)
        shape = shapes[rep]
        
        plot_df.sort_values('upper bound for contig length', inplace=True)
        
        ax=axs
        color = colors[(o2, week)]
        
        ax.plot(plot_df[x], plot_df[y],
                linestyle='-', marker=shape, color=color, label=label,
                alpha = 0.5)
        handles, labels = ax.get_legend_handles_labels()
        lgd = ax.legend(handles, labels, bbox_to_anchor=(1.75, 1.05))

        ax.set_xlabel(x)
        ax.set_ylabel('frac reads assigned to\ncontigs this length')
        ax.axhline(y=0, color='gray', linestyle='-')
        ax.set_xscale('log')

    plt.subplots_adjust(hspace=0.25)
    
    return fig 

def plot_good_vs_bad_low_o2_samples_all_reps(binwidth=1500):
    """
    Replaces plot_good_vs_bad_low_o2_samples().  Plot all Low2 samples for week >=8.
    """
    bins = prepare_bins(binwidth)
    df = summarise_contigs_by_contgig_size(bins)

    df['frac reads assigned to contigs this length'] = df['# mapped reads']/df['total reads (in fastq)']
    print(df.columns)
    df_extract2 = df[(df['oxygen']== 'low') &
                     #((df['replicate']== 2) | (df['replicate']== 1)) &
                     df['week'].isin([8, 9, 10, 11, 12, 13])]
    df_extracts = df_extract2


    fig, axs = plt.subplots(1, 1, figsize=(4, 3)) #, sharex=True, sharey=True)

    green = '#31a354'
    orange = '#d95f0e'
    def assign_color(rep, week):
        if rep == 3 or rep == 4:
            return orange
        # now only reps 1 and 2 remain.  These became good after week 10. 
        if week >= 11:
            return green
        else:
            return orange

    shapes = {1:'o', 2:'o', 3:'.', 4:'.'}

    x = 'upper bound for contig length'
    y = 'frac reads assigned to contigs this length'

    fig.suptitle("fraction of reads assigned to\ncontigs of different lengths", size=14)
    
    def get_unique(df_x, col):
        uniques = df_x[col].unique()
        assert len(uniques) == 1, "too many unique values for {}: {}".format(col, uniques)
        return uniques[0]

    df_extracts.sort_values(['week', 'replicate'], inplace=True)

    for sample_id, plot_df in df_extracts.groupby(['sample id'], sort=False):
        o2 = get_unique(plot_df, 'oxygen')
        rep = get_unique(plot_df, 'replicate')
        week = get_unique(plot_df, 'week')
        label = "{} O2, rep {}, week {}".format(o2, rep, week)
        print(label)
        shape = shapes[rep]

        plot_df.sort_values('upper bound for contig length', inplace=True)

        ax=axs
        color = assign_color(rep, week)

        ax.plot(plot_df[x], plot_df[y],
                linestyle='-', marker=shape, color=color, label=label,
                alpha = 0.5)
        handles, labels = ax.get_legend_handles_labels()
        lgd = ax.legend(handles, labels, bbox_to_anchor=(1.75, 1.7))

        ax.set_xlabel(x)
        ax.set_ylabel('frac reads assigned to\ncontigs this length')
        ax.axhline(y=0, color='gray', linestyle='-')
        ax.set_xscale('log')

    plt.subplots_adjust(hspace=0.25, top=0.8)

    return fig
