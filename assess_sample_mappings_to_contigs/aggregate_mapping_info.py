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
    sums.rename(columns={'# mapped reads':'sum(reads mapped)'}, inplace=True)
    sums['sum(reads mapped)'] = sums['sum(reads mapped)'].fillna(0)
    sums['sum(reads mapped)'] = sums['sum(reads mapped)'].astype(int)
    sums['cryptic metagenome name'] = df['cryptic metagenome name'].iloc[0]
    return sums

def merge_on_meta_info(df):
    # could be faster if I didn't read this every time!
    sample_info = pd.read_csv('../data/sample_info/sample_info.tsv', sep='\t')
    del sample_info['project']
    sample_translation = pd.read_csv('../data/sample_info/meta4_sample_names--cryptic_to_sample_number.tsv', sep='\t')
    sample_translation = sample_translation[['sample id', 'cryptic metagenome name']]
    total_reads_per_sample = pd.read_csv('../data/sample_info/sample_read_counts.tsv', 
                                         names=['cryptic metagenome name', 'total reads (in fastq)'], sep='\t')
    total_reads_per_sample['cryptic metagenome name'] = total_reads_per_sample['cryptic metagenome name'].str.strip('.fastq.gz')

    for sdf in [sample_translation, sample_info, total_reads_per_sample]:
        df = pd.merge(df, sdf)
    return df


def parse_idxstats_bin_coarsely_and_append_results(idxstat_dir, bins, outpath):
    idxstat_files = os.listdir(idxstat_dir)
    idxstat_paths = [os.path.join(idxstat_dir, p) for p in idxstat_files]

    # only has contigs that *do* map to bins
    bin_contig_mappings = pd.read_csv('../abundances/results/bin_contig_mappings.tsv', sep='\t')
    binned_contigs = bin_contig_mappings['contigName'].unique().tolist() # len = 304871
    # bin_contig_mappings.columns = ['contigName', 'file', 'contigs', 'bin', 'Bin Id', 'bin_id']
    bin_contig_mappings = bin_contig_mappings[['contigName', 'bin']]

    files_completed = 0
    for i in idxstat_paths:
        df = parse_idxstat(i)
        
        df['binned contig'] = False
        df.loc[df['contig'].isin(binned_contigs), 'binned contig'] = True

        df_coarsened = coarsen_size_distribution(df=df, bin_edges=bins)

        binned_coarsened = coarsen_size_distribution(df=df[df['binned contig'] == True].copy(), bin_edges=bins)
        binned_coarsened.rename(columns={'sum(reads mapped)':'sum(reads mapped to bins)'}, inplace=True)
        df_coarsened = pd.merge(df_coarsened, binned_coarsened[['hist bin', 'sum(reads mapped to bins)']], how='outer')

        unbinned_coarsened = coarsen_size_distribution(df=df[df['binned contig'] == False].copy(), bin_edges=bins)
        unbinned_coarsened.rename(columns={'sum(reads mapped)':'sum(reads not mapped to bins)'}, inplace=True)
        df_coarsened = pd.merge(df_coarsened, unbinned_coarsened[['hist bin', 'sum(reads not mapped to bins)']], how='outer')

        df_coarsened = merge_on_meta_info(df_coarsened)
        assert df_coarsened.shape[0] > 0, "df_coarsened has no rows left"

        if files_completed == 0:
            df_coarsened.to_csv(outpath, sep='\t', index=False)
        else:
            # append to the already existing file
            df_coarsened.to_csv(outpath, sep='\t', mode='a', header=False, index=False)
        
        files_completed += 1
        print('parse_idxstats_bin_coarsely_and_append_results(): done with file # {}'.format(files_completed))


def individual_contig_summary(idxstat_dir, outpath):
    idxstat_files = os.listdir(idxstat_dir)
    idxstat_paths = [os.path.join(idxstat_dir, p) for p in idxstat_files]

    # only has contigs that *do* map to bins
    bin_contig_mappings = pd.read_csv('../abundances/results/bin_contig_mappings.tsv', sep='\t')
    binned_contigs = bin_contig_mappings['contigName'].unique().tolist() # len = 304871
    # bin_contig_mappings.columns = ['contigName', 'file', 'contigs', 'bin', 'Bin Id', 'bin_id']
    bin_contig_mappings = bin_contig_mappings[['contigName', 'bin']]

    files_completed = 0
    for i in idxstat_paths:
        df = parse_idxstat(i)
        n_rows_before = df.shape[0]

        # drop rows without reads mapped.
        df = df[df['# mapped reads'] != 0]
        n_rows_after = df.shape[0]
        n_dropped = n_rows_before - n_rows_after
        percent_dropped = (n_dropped/n_rows_before)*100
        print('dropped {} rows ({})%'.format(n_dropped, percent_dropped))

        df = merge_on_meta_info(df)
        df['binned contig'] = False
        df.loc[df['contig'].isin(binned_contigs), 'binned contig'] = True

        if files_completed == 0:
            assert df.shape[0] > 0, "df_coarsened has no rows left"
            df.to_csv(outpath, sep='\t', index=False)
        else:
            # append to the already existing file
            df.to_csv(outpath, sep='\t', mode='a', header=False, index=False)
        

        files_completed += 1
        print('individual_contig_summary(): done with file # {}'.format(files_completed))


def mkdir(dirname):
    if not os.path.exists(dirname):
        os.mkdir(dirname)
    

if __name__ == '__main__':
    all_contig_lengths = pd.read_csv('../assembly/final.contigs.len', sep='\t', names = ['contig', 'len'])
    max_size = all_contig_lengths['len'].max()
    print('max size: {}'.format(max_size))
    # slick way to get bins for aggregation:
    # [int(n) for n in np.linspace(0, 100000, 11).tolist()]
    bin_width = 2000
    bin_edges = np.arange(0, max_size + bin_width, bin_width, dtype=int).tolist()
    print(bin_edges[0:10])

    outdir = './results'
    mkdir(outdir)

    outfile_linear_bins = os.path.join(outdir, 'results--linear_bins.tsv')
    parse_idxstats_bin_coarsely_and_append_results(
        '/work/m4b_binning/assembly/map_reads/flagstat/idxstat_results/', 
        bins=bin_edges, outpath=outfile_linear_bins)

    # Repeat with log-scale bins to zoom in on the < 10kb contigs
    outfile_log10_bins= os.path.join(outdir, 'results--log10_bins.tsv')
    # [     10,     100,    1000,   10000,  100000, 1000000]
    magnitude = 10** int(np.ceil(np.log10(max_size)))
    bin_edges = np.logspace(start=1, stop=np.log10(magnitude), num=np.log10(magnitude), base=10, dtype=int).tolist()
    print(bin_edges[0:10])
    parse_idxstats_bin_coarsely_and_append_results(
        '/work/m4b_binning/assembly/map_reads/flagstat/idxstat_results/', 
        bins=bin_edges, outpath=outfile_log10_bins)

    outfile_all_contigs = os.path.join(outdir, 'individual_contig_summary.tsv')
    individual_contig_summary('/work/m4b_binning/assembly/map_reads/flagstat/idxstat_results/', outfile_all_contigs)
    

