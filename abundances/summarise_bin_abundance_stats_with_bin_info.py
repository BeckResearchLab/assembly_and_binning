import pandas as pd

pd.set_option('display.height', 1000)
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

abundances = pd.read_csv('./results/bin_abundances.tsv', sep='\t')
bin_details = pd.read_csv('../bin_info/metabat_bins--all_info.tsv', sep='\t')

# abundances df has bin numbers as stirngs so 'none' was allowed.  Make bin_details' match.
bin_details['bin'] = bin_details['bin number'].astype(str)

results = pd.merge(bin_details, abundances, how='outer')
assert results.shape[0] == abundances.shape[0], 'problem with merge: lost or gain some rows.'

results.to_csv('./results/abundances_with_taxonomy_and_bin_stats.tsv', sep='\t', index=False)

# Make a version of the table with one row per bin, and the maximum, minimum, and mean  abundances in the 88 samples. 
max_abundance = pd.DataFrame(results.groupby(['bin'])['abundance'].max())
min_abundance = pd.DataFrame(results.groupby(['bin'])['abundance'].min())
median_abundance = pd.DataFrame(results.groupby(['bin'])['abundance'].median())
max_abundance.reset_index(inplace=True)
min_abundance.reset_index(inplace=True)
median_abundance.reset_index(inplace=True)
max_abundance.rename(columns={'abundance':'max(abundance)'}, inplace=True)
min_abundance.rename(columns={'abundance':'min(abundance)'}, inplace=True)
median_abundance.rename(columns={'abundance':'med(abundance)'}, inplace=True)

summary_cols = ['Bin Id',  'Completeness', 'Contamination', 'Strain heterogeneity', 'Genome size (bp)',  '# contigs',  'N50 (contigs)',  
    'Mean contig length (bp)',  'Longest contig (bp)', 'GC', 'GC std (scaffolds > 1kbp)', 'Coding density',  '# predicted genes',  
    'bin number', 'bin name', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain', 'confidence', 'bin'] 
#     'max(abundance)', 'min(abundance)', 'med(abundance)']
results_summary = results[summary_cols].copy()
results_summary.drop_duplicates(inplace=True)
assert results_summary.shape[0] == bin_details.shape[0] + 1  # have bin 'none' now. 
results_summary = pd.merge(results_summary, max_abundance)
results_summary = pd.merge(results_summary, min_abundance)
results_summary = pd.merge(results_summary, median_abundance)

sort_cols = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain', 'max(abundance)', 'Completeness', 'Contamination']
sort_orders = [1]*(len(sort_cols)-3)
sort_orders.extend([0]*3) # make abundance, Completeness, and Contamination sort by descending.
results_summary.sort_values(by=sort_cols, ascending=sort_orders, inplace=True)
results_summary.to_csv('./results/abundances_stats_for_each_bin_with_taxonomy_and_bin_stats.tsv', sep='\t', index=False)

# summarise at the family level
# Need to change the nan to "none" so it isn't dropped in the groupby
family_level = results.copy()
family_level['family'].fillna('none', inplace=True)
family_sums = pd.DataFrame(family_level.groupby(['sample id', 'family'])['abundance'].sum()) 
family_sums.reset_index(inplace=True)
family_summary = pd.merge(family_level[['sample id', 'oxygen', 'replicate', 'week', 'family']].drop_duplicates(), family_sums) 
# In [49]: results[results['family'].isnull()][['bin', 'sample id', 'family']]['bin'].unique()
# Out[49]: array(['none'], dtype=object)
family_summary.to_csv('./results/summary_abundancy_by_family.tsv', sep='\t', index=False)
assert family_summary.groupby('sample id')['abundance'].sum().min() > 0.999
assert family_summary.groupby('sample id')['abundance'].sum().max() < 1.001


