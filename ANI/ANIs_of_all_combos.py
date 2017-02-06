import glob
import itertools
import os
import re
import subprocess
import pandas as pd

def basename_without_fa_type(string):
    bn = os.path.basename(string)
    bn = re.sub('\.ffn', '', bn) # * is zero or more matches
    return re.sub('\.fa*', '', bn) # * is zero or more matches

def assess_one_pair(fasta1, fasta2, results_dir):
    f1bn = basename_without_fa_type(fasta1)
    f2bn = basename_without_fa_type(fasta2)
    fname = f1bn + '--' + f2bn
    print(fname)
    #evo = '/work/software/enveomics/Scripts/ani.rb'
    jgi_ani = '/work/software/ANIcalculator_v1/ANIcalculator'
    stdout = os.path.join(results_dir, fname + '.out')
    stderr = os.path.join(results_dir, fname + '.err')
    txt = os.path.join(results_dir, fname + '.txt')
    tsv = os.path.join(results_dir, fname + '.tsv')
    out = open(stdout, 'w')
    err = open(stderr, 'w')
    out.write(fname + '\n')
    err.write(fname + '\n')
    out.flush()
    
    #subprocess.check_call([evo, '--seq1', fasta1, '--seq2', fasta2, '--tab', tsv, '--res', txt], stdout=out, stderr=err) 
    subprocess.check_call([jgi_ani, '-genome1fna', fasta1, '-genome2fna', fasta2, '-outfile', tsv], #'--res', txt], 
                            stdout=out, stderr=err) 
    out.close()
    err.close()

def genomes_from_filename(fname):
    fname = os.path.basename(fname)
    fname = re.sub('\.tsv', '', fname)
    print(fname.split('--'))
    return fname.split('--')

def parse_results(results_dir):
    results = pd.DataFrame()
    fnames = os.listdir(results_dir)
    #fnames = [f for f in fnames if '.tsv' in f]
    fnames = [os.path.join(results_dir, f) for f in fnames if '.tsv' in f]
    for fname in fnames:
        print(fname)
        f1, f2 = genomes_from_filename(fname)
        #deets = pd.read_csv(fname, sep='\t', names=['ANI', 'ANI std', 'fragments used', 'fragments in the smallest genome'])
        deets = pd.read_csv(fname, sep='\t') 
        deets['bin 1'] = f1
        deets['bin 2'] = f2
        print('Done parsing {}'.format(fname))
        results = pd.concat([results, deets], axis=0)
    return results
    

if __name__ == "__main__":
    # For some reason the tool doesn't write .tsv files if you don't delete this:
    os.rmdir('./ani.blast.dir')

    outdir = 'ANIs_by_JGI_tool'
    if os.path.exists(outdir):
        print('clear out {}'.format(outdir))
        files_to_remove = glob.glob(outdir + '/*')
        [os.remove(f) for f in files_to_remove]
    else:
        os.mkdir(outdir)

    places_to_look = ['/work/m4b_binning/assembly/isolate_Fauzi_stats/isolate/genomes/prokka/results/*/*.ffn',
                      '/work/m4b_binning/assembly_mycc_completed/assembly_mycc_2.5kb_56mer_lt_0.4_st_50/20170131_1745_56mer_0.7_cov/prokka/results/bin_*/*.ffn']
    fastas = []
    for place in places_to_look:
        fastas += glob.glob(place)

    print('gene-based fastas found by glob: {}'.format(len(fastas)))
    print('sample: {}'.format(fastas[0:3]))

    combos = list(itertools.combinations(fastas, 2))
    print('# of combinations for {} fastas: {}'.format(len(fastas), len(combos)))

    for (fasta1, fasta2) in combos:
        assess_one_pair(fasta1, fasta2, outdir)

    # Now do all the self-self to get the expected 100%
    print('---- do self:self combos ----')
    for (fasta1, fasta2) in zip(fastas, fastas):
        assess_one_pair(fasta1, fasta2, outdir)
    
    print('parse the output tsv files')
    results = parse_results(outdir)
    results.to_csv('ANI_summary.tsv', sep='\t', index=False)
    
    

