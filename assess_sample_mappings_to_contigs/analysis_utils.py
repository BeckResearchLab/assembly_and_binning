# Import matplotlib before seaborn
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import pandas as pd
import seaborn as sns

def load_data():
    df_all = pd.read_csv('./results/individual_contig_summary.tsv', sep='\t')
    df_all['frac reads for contig'] = df_all['# mapped reads']/df_all['total reads (in fastq)']
    return df_all 

def load_only_binned_or_unbinned(binned=False, groupby_sample=False):
    df = load_data()
    df = df[df['binned contig'] == True]
    if groupby_sample:
        grouped = pd.DataFrame(df.groupby('sample id')['frac reads for contig'].sum())
        grouped.reset_index(inplace=True)
        grouped.rename(columns={'frac reads for contig':'sum(frac reads for contig)'}, inplace=True)
        return grouped

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

    
