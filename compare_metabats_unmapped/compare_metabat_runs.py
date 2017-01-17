import pandas as pd

md = pd.read_csv('../abundances/bin_abundances.tsv', sep='\t') # default params
ms = pd.read_csv('../metabat.specific/bin_abundances.tsv', sep='\t')
ms2500 = pd.read_csv('../metabat.specific.min_contig2500/bin_abundances.tsv', sep='\t')

md['metabat run'] = 'default settings'
ms['metabat run'] = 'specific'
ms2500['metabat run'] = 'specific 2500'

results = pd.DataFrame()
for df in [md, ms, ms2500]:
    none_rows = df[df['bin'] == 'none']
    assert none_rows.shape[0] == 88, "expected 88 rows; got {}".format(df.shape[0])
    # get just the none rows.  The bins with numbers shouldn't match up anyway
    results = pd.concat([results, none_rows], axis=0)

assert results.shape[0] == 3*88, "problem with merge: expected 3*88 rows"

results.to_csv('./abundances_of_none.tsv', sep='\t', index=False)

# Filter to only the rows with bin ==' none'
