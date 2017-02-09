import argparse
import os

from Bio import SeqIO

def batch_iterator(iterator, batch_size):
    """
    From http://biopython.org/wiki/Split_large_file
    Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                #entry = iterator.next()  # python 2
                #entry = iterator.__next__()   # python 3
                entry = next(iterator)   # python 3
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch


def break_apart_fasta(fasta_path, num_contigs_per_chunk):
    """
    Split up a multi-fasta into separate files with a specified number of contigs in each.
    Based on http://biopython.org/wiki/Split_large_file
    """
    fasta_basename = os.path.splitext(os.path.basename(fasta_path))[0]
    dirname = fasta_basename + "_chunks_of_{}".format(num_contigs_per_chunk)
    if not os.path.exists(dirname):
        os.mkdir(dirname)

    record_iter = SeqIO.parse(open(fasta_path),"fasta")
    for i, batch in enumerate(batch_iterator(record_iter, num_contigs_per_chunk)):
        filename = "{}_group_{}.fa".format(fasta_basename, i + 1)
        outpath = os.path.join(dirname, filename)
        handle = open(outpath, "w")
        count = SeqIO.write(batch, handle, "fasta")
        handle.close()
        print("Wrote %i records to %s" % (count, filename))
    return dirname


def verify_sum_of_contigs_adds_up(original, folder_of_split_files):
    print('count sequences in original file: {}'.format(original))
    record_dict = SeqIO.index(original, 'fasta')
    original_num = len(record_dict)
    print("number of contigs in original file: {}".format(original_num))

    split_files = [os.path.join(folder_of_split_files, f) for f in os.listdir(folder_of_split_files)]
    print('count contigs for {}'.format(split_files))
    split_num = 0 
    for sf in split_files:
        record_dict = SeqIO.index(sf, 'fasta')
        split_num += len(record_dict)

    print("number of contigs in split up files: {}".format(split_num))
    assert original_num == split_num, \
        "the number of contigs in the split-up files does not match the number in the original!"

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("original", help="multi-fasta you wish to split up")
    parser.add_argument("num_contigs", type=int, 
        help="number of contigs to include in each reduced-size file")
    args = parser.parse_args()

    if (args.original is None) or (args.num_contigs is None):
        args.print_help()
    else:
        # Break up the files and return the dirname
        dirname = break_apart_fasta(args.original, args.num_contigs)
        # Check that the number of contigs produced is right. 
        verify_sum_of_contigs_adds_up(args.original, dirname)

    

