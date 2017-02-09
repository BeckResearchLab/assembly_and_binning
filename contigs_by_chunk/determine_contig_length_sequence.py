import argparse
import os

from Bio import SeqIO

def get_contig_lengths(fasta_path):
    """
    Save the sequence of contig lengths to a file.
    Allows you to check that they aren't in a particular order. 
    """
    fasta_basename = os.path.basename(fasta_path)
    fasta_basestring = os.path.splitext(fasta_basename)[0]
    outfile = fasta_basestring + '_contig_lengths_in_order.tsv'
    handle = open(outfile, "w")

    fasta = SeqIO.parse(open(fasta_path),"fasta")
    for i, contig in enumerate(fasta):
        print(i+1, len(contig.seq))
        handle.write('{}    {}\n'.format(i+1, len(contig.seq)))
    handle.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta", action='store',  help="multi-fasta to determine contig lengths for")
    args = parser.parse_args()
    print(args.fasta)

    if args.fasta is None:
        parser.print_help()
    else:
        get_contig_lengths(args.fasta)

    

