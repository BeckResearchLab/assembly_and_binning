import os

import pandas as pd
#import seaborn as sns

total_reads = pd.read_csv('../data/sample_info/sample_read_counts.tsv', sep='\t', names = ['fastq filename', 'number of reads'])
total_reads['cryptic metagenome name'] = total_reads['fastq filename'].str.strip('.fastq.gz')
sample_info = pd.read_csv('../data/sample_info/sample_info.tsv', sep='\t')
sample_translation = pd.read_csv('../data/sample_info/meta4_sample_names--cryptic_to_sample_number.tsv', sep='\t')
read_mappings = pd.read_csv('./data/num_reads_mapped--can_double_count_multiple_mappings.tsv', sep='\t')

reads = pd.merge(sample_info, sample_translation) 
reads = pd.merge(reads, total_reads)
reads = pd.merge(reads, read_mappings)

out_path = 'total_num_reads_across_samples_with_sample_info.tsv'  
out_dir = './data'
reads.to_csv(os.path.join(out_dir, out_path), sep='\t', index=False)

