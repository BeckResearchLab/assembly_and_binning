import argparse
import os 
import string

from Bio import SeqIO

def new_filename(n_min):
    return 'contigs_longer_than_{}bp.fa'.format(n_min)

def omit_shorter_than_n_bp(original_path, n_min):
    outname = new_filename(n_min)
    print('save file with contigs longer than {}bp as {}'.format(n_min, outname))
    outfile = open(outname, 'w')

    for record in SeqIO.parse(original_path, "fasta"):
        record_id = record.id 

        # get rid of the shorties
        if len(record.seq) <= n_min:
            continue

        outfile.write("> " + record.id + "\n")
        outfile.write(str(record.seq))
        outfile.write('\n')

    outfile.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Trim out contigs shorter than the specified length from a fasta file.')
    parser.add_argument('fasta', type=str, help='fasta source to pick contigs from')
    parser.add_argument('--min', dest='min', type=int, help='minimum contig size for inclusion, bp')
    args = parser.parse_args()

    if args.min is None:
        print(parser.print_help())
    else:
        omit_shorter_than_n_bp(args.fasta, args.min)
    
