import argparse
import os
import re

import pandas as pd

def parse_phylophlan(path):
    # make flexible fore metabat or mycc
    df = pd.read_csv(path, sep='\t', names=['bin_date', 'taxonomy'])
    #df['bin name'] = df['bin_date'].str.extract('(metabat_bin_[0-9]+)_[0-9]+')
    #df['bin number'] = df['bin_date'].str.extract('metabat_bin_([0-9]+)_[0-9]+')
    #df['date'] = df['bin_date'].str.extract('metabat_bin_[0-9]+_([0-9]+)')
    df['bin name'] = df['bin_date'].str.extract('([A-z]+_bin_[0-9]+)_[0-9]+')
    df['bin number'] = df['bin_date'].str.extract('[A-z]_bin_([0-9]+)_[0-9]+')
    df['date'] = df['bin_date'].str.extract('[A-z]+_bin_[0-9]+_([0-9]+)')
    return df

def taxonomy_to_columns(df):
    taxonomy = df['taxonomy'].str.split('__').apply(pd.Series)
    taxonomy.columns = ['domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']
    # remove the '.c' type suffixes
    for c in taxonomy.columns:
        taxonomy[c] = taxonomy[c].str.replace(r'\.[a-z]', '')
    return taxonomy

def confidence_from_path(path):
    m = re.search('imputed_conf_([a-z]+)', path) # note: the low ends with '_low' and the others end with '-high', '-med'
    if len(m.groups()) == 1:
        return m.groups()[0]
    else:
        return "??"

def parse_and_tidy(path):
    df = parse_phylophlan(path)
    taxonomy = taxonomy_to_columns(df)
    df = pd.concat([df, taxonomy], axis=1)
    df['confidence'] = confidence_from_path(path) 
    return df

def parse_all_confidence_levels(path_to_files):
    dfs = list()
    taxa_files = os.listdir(path_to_files)
    taxa_files = [t for t in taxa_files if 'imputed_conf' in t]
    taxa_paths = [os.path.join(path_to_files, t) for t in taxa_files]
    print('parsing taxonomy from files {}'.format(taxa_files))
    for f in taxa_paths:
        dfs.append(parse_and_tidy(f))
    results = pd.concat(dfs, axis=0)
    #if len(results['8'].unique()) == 1:
    #    results.drop('8', 1, inplace=True)
    if len(results['domain'].unique()) == 1:
        results.drop('domain', 1, inplace=True)
    results.drop('taxonomy', 1, inplace=True)
    results.drop('bin_date', 1, inplace=True)
    # keep only the highest confidence
    results['conf score'] = 0
    results.loc[results['confidence'] == 'low', 'conf score'] = 1
    results.loc[results['confidence'] == 'medium', 'conf score'] = 2
    results.loc[results['confidence'] == 'high', 'conf score'] = 3
    return results
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(
                prog='parse_phylophlan',
                description="""Script parses phylophlan taxonomy calls to .tsv""")

    parser.add_argument("-input", dest="input", type=str, help="Specify directory with phylophlan taxonomy files")
    parser.add_argument("-output", dest="output", type=str, help="output tsv filename")

    args = parser.parse_args()

    if args.input is None:
        print(parser.print_help())
    else:
        df = parse_all_confidence_levels(args.input)
        df.to_csv(args.output, index=False, sep='\t')
    

