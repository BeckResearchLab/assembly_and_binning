import os
import pandas as pd
import subprocess

source_dir = '../nucleotide'
files = os.listdir(source_dir)
print(len(files))
fastas = [f for f in files if ".fa" in f]
print(len(fastas))

out_path = './' 
if not os.path.exists(out_path):
    os.mkdir(out_path)

def bin_info(filename, source_dir, out_path):
    substrings = filename.split('.')
    name = '{}'.format(substrings[0])
    info = {'filename':filename, 
            'name':name}
    info['out_dir'] = os.path.join(out_path, name)
    info['path'] = os.path.join(source_dir, info['filename'])
    info['stdout'] = os.path.join(out_path, '{}/{}.out'.format(name, name))
    info['stderr'] = os.path.join(out_path, '{}/{}.err'.format(name, name))
    return info

summary = pd.DataFrame()

for fasta in fastas:
    bin_info_dict = bin_info(fasta, source_dir, out_path)
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
        isolate_name = row['name']
        isolate_name_short = row['name'][0:2]
        out_dir = row['out_dir']
        path = row['path']
        stdout_path = row['stdout']
        stderr_path = row['stderr']
        print("cwd: {}".format(os.getcwd()))
        print("stdout_file: {}".format(stdout_path))
        stdout_file = open(stdout_path, 'w+')
        stderr_file = open(stderr_path, 'w+')
    
        # replace when we get more info
        locus_tag = '{}'.format(isolate_name_short)

        command = ['prokka', path, 
                   # add arg for 
                    # [19:39:22] Loading and checking input file: ../nucleotide/Methylotenera_sp_G11.fa
                    # [19:39:22] Contig ID must <= 20 chars long: GQ51DRAFT_unitig_0_quiver_dupTrim_10755.1
                    # [19:39:22] Please rename your contigs or use --centre XXX to generate clean contig names.
                    #'--compliant', 
                    '--centre C --locustag L', 
                   #'--centre', 'X',
                   '--outdir', out_dir, 
                   '--force', 'ON',  # force write over previous results; done so stdout, stderr handled more easily 
                   #'--locustag', locus_tag # add more detail later
                    ]
        print(' '.join(command))
        subprocess.check_call(command, stdout=stdout_file, stderr=stderr_file)
        stdout_file.close()
        stderr_file.close()

