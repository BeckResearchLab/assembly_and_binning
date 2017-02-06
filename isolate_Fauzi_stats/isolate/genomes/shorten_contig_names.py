import os 
import random
import string

from Bio import SeqIO

def generate_n_random_chars(n):
    rand_string = ''
    for i in range(n):
        rand_string += random.choice(string.letters + string.digits)
    return rand_string

def shorten_one(original_path, dest_path, str_len):
    n_rand_chars = 3
    outfile = open(dest_path, 'w')

    for record in SeqIO.parse(original_path, "fasta"):
        record_id = record.id 
        print('shorten record.id from {} to {}:'.format(len(record_id), str_len))
        record_id_short = record_id[0:str_len - n_rand_chars]
        random_letters = generate_n_random_chars(n_rand_chars)
        print('generated random letters for uniquie contig name: {}'.format(random_letters))
        record_id_short += random_letters 
        print('{} --> {}'.format(record_id, record_id_short))

        #outfile.write(">"+record.description+"\n")
        print('record description: {}'.format(record.description))
        outfile.write("> " + record_id_short + "\n")
        outfile.write(str(record.seq))
        outfile.write('\n')

    outfile.close()
    pass

def shorten_all(fasta_dir, out_dir):
    fastas = os.listdir(fasta_dir)
    fastas = [f for f in fastas if ".fa" in f]

    for f in fastas:
        in_path = os.path.join(fasta_dir, f)
        out_path = os.path.join(out_dir, f)
        shorten_one(in_path, out_path, 20)

if __name__ == '__main__':

    originals_dir='./nucleotide'
    out_dir='/work/m4b_binning/assembly/isolate_Fauzi_stats/isolate/genomes/nucleotide--short_contig_names'

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    shorten_all(originals_dir, out_dir)
    
