import pandas as pd

# sample_info.xls was copied from /dacb/meta4_iso/analysis/sample_info.xls on waffle
# contains some stuff that isn't specific to this project, and 
# one of the strings we want that maps sample names is embeded in a larger string.

df = pd.read_csv('./sample_info_orig.xls', sep='\t')
del df['mutation_locus']
del df['path_to_proteome_FASTA']
del df['is_qc']
del df['include']
del df['path_to_alignment']
del df['genes_table']
del df['path_to_GFF']
del df['path_to_genome_FASTA']
del df['path_to_workspace']

df['JGI sample id'] = df['path_to_FASTQ'].str.extract('/Raw_Data/([.0-9A-Z]+).fastq.gz')
del df['path_to_FASTQ']

df['O2'] = df['O2'].str.lower()

df.to_csv('./sample_info.tsv', sep='\t')

