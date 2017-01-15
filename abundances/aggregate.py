import argparse
import os
import numpy as np
import pandas as pd

pd.set_option('display.height', 1000)
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

def check_for_needed_files(*args):
    for f in args:
        assert os.path.exists(f), "file {} doen't exist".format(f)

def aggregate(contig_path, depth_path, id_trans_path,  sample_details_path):

    contigs = pd.read_csv(contig_path, sep='\t', dtype={'bin':str})

    # get coverage
    coverage = pd.read_csv(depth_path, sep='\t') # takes a while. shape = (2,617,225, 179)
    # rename columns from 8777.4.112218.TGCCAT.fastq.gz.mapped.sorted.bam to 8777.4.112218.TGCCAT
    coverage.columns = coverage.columns.str.replace('.fastq.gz.mapped.sorted.bam', '')
    keep_cols = [c for c in coverage.columns if '-var' not in c]
    coverage = coverage[keep_cols]
    coverage = pd.melt(coverage, id_vars=['contigName', 'contigLen', 'totalAvgDepth'], var_name='cryptic metagenome name', value_name='coverage')
    coverage = pd.merge(coverage, contigs[['contigName', 'bin']], how='outer')
    # ??? DO I WANT TO DO THIS?? 
    # ?? remove contigs that were not assigned to a bin? 
    coverage['bin'].fillna('none', inplace=True)
    #coverage.loc[coverage['bin'].isnull(), 'bin'] = 'none' # reassign
    #coverage = coverage


    # Normalize & Aggregate: compute fraction of reads going to each bin, for each sample. 
    # Does not average for genome size. 
    coverage['len*coverage'] = coverage['contigLen']*coverage['coverage'] 

    # get total coverage for each sample. 
    sample_sums = pd.DataFrame(coverage.groupby('cryptic metagenome name')['len*coverage'].sum())
    sample_sums.reset_index(inplace=True)
    sample_sums.rename(columns={'len*coverage':'sum(len*coverage), sample'}, inplace=True)
    coverage = pd.merge(coverage, sample_sums, how='outer')

    # aggregate across contigs
    abundance = pd.DataFrame(coverage.groupby(['bin', 'cryptic metagenome name'])['len*coverage'].sum())
    abundance.reset_index(inplace=True)
    abundance.rename(columns={'len*coverage':'bin_sum(len*coverage)'}, inplace=True)
    # merge on the total sums:
    abundance = pd.merge(abundance, sample_sums, how='outer')
    abundance['abundance'] = abundance['bin_sum(len*coverage)']/abundance['sum(len*coverage), sample']

    # check that each sample sums to one. 
    sample_sum_errors = pd.Series(abundance.groupby(['cryptic metagenome name'])['abundance'].sum() - np.ones(88, )) 
    assert(sample_sum_errors.max() < 0.0001)  # was on order of 1E-13 when I checked by hand

    # finally, merge on the info to get out of 'cryptic metagenome name' mess
    cryptic_id_translations = pd.read_csv(id_trans_path, sep='\t')
    sample_details = pd.read_csv(sample_details_path, sep='\t')
    assert sample_details.shape[0] == 88  # there is an 83-row problem version of the file
    sample_details = pd.merge(cryptic_id_translations, sample_details, how='outer')
    assert sample_details.shape[0] == 88  # check one more time 

    # merge sample info to abundances for final results:
    results = pd.merge(abundance, sample_details, how='outer')
    assert(abundance.shape[0] == results.shape[0]), "should have only added columns when merging on sample info"
    check_sample_sums = results.groupby('sample id')['abundance'].sum()
    assert check_sample_sums.min() > 0.999, 'not all samples have sum(abundances) == 1: some are fall short'
    assert check_sample_sums.max() < 1.001, 'not all samples have sum(abundances) == 1: some have too much'
    
    return results


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                prog='summarise_bin_abundance_for_each_sample',
                description="""Script identifies which contigs belong to which bins and 
                               summarises depth.txt for each bin in each sample""")

    parser.add_argument("-contig_path", dest="contig_path", type=str, help="Specify the tsv that links contigs to bins")
    parser.add_argument("-depth_path", dest="depth_path", type=str, help="path to metabat's depth.txt file")
    parser.add_argument("-sample_info_dir", dest="sample_info_dir", type=str, help="path to sample info folder")

    args = parser.parse_args()

    if args.contig_path is None:
        print(parser.print_help())
    else:
        id_trans_path = os.path.join(args.sample_info_dir, 'meta4_sample_names--cryptic_to_sample_number.tsv')
        sample_details_path =  os.path.join(args.sample_info_dir, 'sample_info.tsv')
        # check that all the files exist:
        check_for_needed_files(args.contig_path, args.depth_path, id_trans_path,  sample_details_path)
        
        df = aggregate(contig_path = args.contig_path,
                       depth_path = args.depth_path, 
                       # need two files:   '../data/sample_info/meta4_sample_names--cryptic_to_sample_number.tsv', '../data/sample_info/sample_info.tsv',   
                       # for now, assume those paths will stay fixed and only pass the folder name as an arg. 
                       id_trans_path = id_trans_path,
                       sample_details_path = sample_details_path)
        df.to_csv('bin_abundances.tsv', sep='\t', index=False)
