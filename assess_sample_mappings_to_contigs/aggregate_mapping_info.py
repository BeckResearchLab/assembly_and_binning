import os 
import re

import numpy as np

import pandas as pd

# files like 
# /work/m4b_binning/assembly/map_reads/flagstat/idxstat_results/8777.4.112218.TGCCAT.fastq.gz.mapped.sorted.bam.idxstats

def parse_idxstat(idxstat_path):
    # The output is TAB-delimited with each line consisting of 
    # reference sequence name, sequence length, # mapped reads and # unmapped reads. 
    colnames = ['contig', 'contig length', '# mapped reads', '# unmapped reads']
    df = pd.read_csv(idxstat_path, sep='\t', names=colnames)
    df['cryptic metagenome name'] = re.search('([0-9]+.[0-9]+.[0-9]+.[ACTG]+)',idxstat_path).group(1)
    return df

def coarsen_size_distribution(df, bin_edges):
    df['hist bin'] = pd.cut(df['contig length'], bins=bin_edges)
    sums = pd.DataFrame(df.groupby('hist bin')['# mapped reads'].sum())
    sums.reset_index(inplace=True)
    # clean it up
    sums.rename(columns={'# mapped reads':'sum(mapped reads)'}, inplace=True)
    sums['sum(mapped reads)'] = sums['sum(mapped reads)'].fillna(0)
    sums['sum(mapped reads)'] = sums['sum(mapped reads)'].astype(int)
    sums['cryptic metagenome name'] = df['cryptic metagenome name'].iloc[0]
    return sums

def parse_idxstats_bin_and_append_results(idxstat_dir, bins, outpath):
    idxstat_files = os.listdir(idxstat_dir)
    idxstat_paths = [os.path.join(idxstat_dir, p) for p in idxstat_files]

    sample_info = pd.read_csv('../data/sample_info/sample_info.tsv', sep='\t')
    del sample_info['project']
    sample_translation = pd.read_csv('../data/sample_info/meta4_sample_names--cryptic_to_sample_number.tsv', sep='\t')
    sample_translation = sample_translation[['sample id', 'cryptic metagenome name']]

    files_completed = 0
    for i in idxstat_paths:
        df = parse_idxstat(i)
        df_coarsened = coarsen_size_distribution(df=df, bin_edges=bins)
        df_coarsened = pd.merge(df_coarsened, sample_translation)
        df_coarsened = pd.merge(df_coarsened, sample_info)

        if files_completed == 0:
            df_coarsened.to_csv(outpath, sep='\t')
        #elif files_completed > 3:
        #    return
        else:
            # append to the already existing file
            df_coarsened.to_csv(outpath, sep='\t', mode='a', header=False)
        
        files_completed += 1
        print('done with file # {}'.format(files_completed))


if __name__ == '__main__':
    all_contig_lengths = pd.read_csv('../assembly/final.contigs.len', sep='\t', names = ['contig', 'len'])
    max_size = all_contig_lengths['len'].max()
    magnitude = 10** int(np.ceil(np.log10(max_size)))
    n_bins = 100
    # slick way to get bins for aggregation:
    # [int(n) for n in np.linspace(0, 100000, 11).tolist()]
    bin_edges = [int(n) for n in np.linspace(0, magnitude, n_bins).tolist()]

    parse_idxstats_bin_and_append_results(
        '/work/m4b_binning/assembly/map_reads/flagstat/idxstat_results/', 
        bins=bin_edges, outpath='./results--all_contigs.tsv')

