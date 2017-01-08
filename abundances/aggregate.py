import pandas as pd

checkm = pd.read_csv('CheckM.tsv', header=1, sep='\t')
contigs = pd.read_csv('bin_contig_mappings.tsv', sep='\t')
info = pd.merge(checkm, contigs, how='outer').head()

# get coverage
coverage = pd.read_csv('../metabat/depth.txt', sep='\t')
import pdb; pdb.set_trace()

