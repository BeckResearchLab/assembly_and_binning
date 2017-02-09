import argparse
import os 
import string

from Bio import SeqIO

def new_filename(n_min=None, n_max=None):
    if n_min is None:
        return 'contigs_shorter_than_{}bp.fa'.format(n_max)
    if n_max is None:
        return 'contigs_longer_than_{}bp.fa'.format(n_min)
    if (n_min is not None) and (n_max is not None):
        return 'contigs_sized_{}_to_{}bp.fa'.format(n_min, n_max)

def omit_shorter_than_n_bp(original_path, n_min, n_max):

    assert (n_min is not None) or (n_max is not None), 'need to specify either min, max, or both min and max contig size(s)'

    outname = new_filename(n_min, n_max)
    print('save file with contigs longer than {}bp and shorter than {} as {}'.format(n_min, n_max, outname))
    outfile = open(outname, 'w')

    for record in SeqIO.parse(original_path, "fasta"):
        record_id = record.id 

        # get rid of the shorties
        if (n_min is not None) and (len(record.seq) < n_min):
            continue
        if (n_max is not None) and (len(record.seq) > n_max):
            continue

        outfile.write("> " + record.id + "\n")
        outfile.write(str(record.seq))
        outfile.write('\n')

    outfile.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Trim out contigs shorter than the specified length from a fasta file.')
    parser.add_argument('fasta', type=str, help='fasta source to pick contigs from')
    parser.add_argument('--min', dest='min', type=int, help='minimum contig size for inclusion, bp')
    parser.add_argument('--max', dest='max', type=int, help='maximum contig size for inclusion, bp')
    args = parser.parse_args()

    if args.min is None and args.max is None:
        print(parser.print_help())
    else:
        omit_shorter_than_n_bp(args.fasta, args.min, args.max)
    
