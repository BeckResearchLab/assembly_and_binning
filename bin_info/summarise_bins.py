import pandas as pd

pd.set_option('display.height', 1000)
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)


phylophlan = pd.read_csv('./phylophlan_result.tsv', sep='\t')
# OLD: checkm = pd.read_csv('./CheckM_results.tsv', sep='\t')
checkm = pd.read_csv('./CheckM_qa.txt', sep='\t')
checkm['bin number'] = checkm['Bin Id'].str.extract(r'bin.([0-9]+)').astype(int)
checkm['bin number'] = checkm['Bin Id'].str.extract(r'bin.([0-9]+)').astype(int)
summary = pd.merge(checkm, phylophlan, on='bin number')
sort_cols = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain', 'Completeness', 'Contamination']
sort_orders = [1]*(len(sort_cols)-2)
sort_orders.append(0)
sort_orders.append(0)
summary.sort_values(by=sort_cols, ascending=sort_orders, inplace=True)
fav_cols = ['Bin Id', 'Marker lineage', '# genomes', '# markers', '# marker sets', 'Completeness', 'Contamination', 'Strain heterogeneity', 'bin number', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain', 'confidence']

summary.to_csv('metabat_bins--all_info.tsv', sep='\t', index=False)
summary[fav_cols].to_csv('metabat_bins.tsv', sep='\t', index=False)

    

