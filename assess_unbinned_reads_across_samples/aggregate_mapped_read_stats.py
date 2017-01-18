import os
import re

import pandas as pd

def parse_flagstat(flagstat_path):
    """
    Count the number of mapped reads.
    Double counts a read that maps two places (https://www.biostars.org/p/138116/)
    # pick out the line like `391696 + 0 mapped (86.31%:-nan%) 0 + 0 paired in sequencing`
    """
    f = open(flagstat_path, 'rt')  # *r*ead *t*ext
    for line in f:
        if 'mapped' in line:
            s = re.search('([0-9]+) \+ [0-9]+', line)
            assert len(s.groups()) == 1, 'expected one match.  Got {}'.format(s.groups())
            return int(s.group(1))

def parse_all_flagstats(flagstat_path_list):
    info = pd.DataFrame()
    for f in flagstat_path_list:
        bn = os.path.basename(f)
        row_dict = {'flagstat file': bn}
        row_dict['reads mapped (includes multiply mapped)'] = parse_flagstat(f)
        df_row = pd.DataFrame({k:[v] for k, v in row_dict.items()})
        info = pd.concat([info, df_row], axis=0)
    return info

if __name__ == '__main__':
    flagstat_dir = '../map_reads/flagstat/flagstat_results/'
    flagstat_paths = os.listdir('../map_reads/flagstat/flagstat_results/')
    flagstat_paths = [f for f in flagstat_paths if 'flagstat' in f]
    # put the dirname back on
    flagstat_paths = [os.path.join(flagstat_dir, f) for f in flagstat_paths if 'flagstat' in f]
    df = parse_all_flagstats(flagstat_paths)
    assert df.shape[0] == 88, "expected 88 rows (one per sample), but got {}".format(df.shape[0])
    # flagstat file names are like '8776.5.111705.ACTGAT.fastq.gz.mapped.sorted.bam.flagstat'
    df['cryptic metagenome name'] = df['flagstat file'].str.strip('.fastq.gz.mapped.sorted.bam.flagstat')

    outdir = './data'
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    outpath = os.path.join(outdir, 'num_reads_mapped--can_double_count_multiple_mappings.tsv')
    df.to_csv(outpath, sep='\t', index=False)

