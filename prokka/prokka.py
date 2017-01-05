import os
import pandas as pd
import subprocess

files = os.listdir('../metabat')
print(len(files))
fastas = [f for f in files if (".fa" in f) and ('bin.' in f)]
print(len(fastas))

if not os.path.exists('./results'):
    os.mkdir('./results')

def bin_info(filename):
    substrings = filename.split('.')
    name = 'bin_{}'.format(substrings[1])
    info = {'filename':filename, 
            'number':substrings[1], 
            'name':name}
    info['out_dir'] = os.path.join('./results', name)
    info['path'] = os.path.join('../metabat', info['filename'])
    info['stdout'] = 'results/{}/{}.out'.format(name, name)
    info['stderr'] = 'results/{}/{}.err'.format(name, name)
    return info

summary = pd.DataFrame()

for fasta in fastas:
    bin_info_dict = bin_info(fasta)
    bin_info_df = pd.DataFrame({k:[v] for k, v in bin_info_dict.items()})
    summary = pd.concat([summary, bin_info_df], axis=0)

print(summary)
summary.to_csv('summary.csv')

for index, row in summary.iterrows():

    # what's in the destination path?
    path = row['out_dir']
    if not os.path.exists(path):
        os.mkdir(path)

    # check for previous completion:
    existing_files = os.listdir(path)
    # sample completed dir has 13 files: .gff, .out, .fna, .err, .gbk, .txt, .faa, .err, .log, .sqn, .tbl, .ffn, .fsa
    gff_exists = len([fn for fn in existing_files if '.gff' in fn]) == 1
    faa_exists = len([fn for fn in existing_files if '.faa' in fn]) == 1

    if (len(existing_files) > 10) and gff_exists and faa_exists: 
        # don't re-run that bin
        pass
    else:
        filename = row['filename']
        bin_name = row['name']
        out_dir = row['out_dir']
        path = row['path']
        bin_number = row['number']
        stdout_path = row['stdout']
        stderr_path = row['stderr']
        stdout_file = open(stdout_path, 'w+')
        stderr_file = open(stderr_path, 'w+')
    
        # replace when we get more info
        locus_tag = 'metabat_bin_{}'.format(bin_number)

        command = ['prokka', path, 
                   '--outdir', out_dir, 
                   '--force', 'ON',  # force write over previous results; done so stdout, stderr handled more easily 
                   '--locustag', locus_tag, 
                   '--genus', '"Genus TBD"',  # add more detail later
                   '--species', '"species TBD"',  # add more detail later
                   '--strain', '"strain TBD"'] # add more detail later
        print(' '.join(command))
        subprocess.check_call(command, stdout=stdout_file, stderr=stderr_file)
        stdout_file.close()
        stderr_file.close()

