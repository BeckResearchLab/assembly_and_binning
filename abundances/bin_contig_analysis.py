import argparse 
import os 
import re

import pandas as pd

def summarise_contigs(fasta):
    ffile = open(fasta, 'rt')  # *r*ead *t*ext
    headers = list()
    # gather the headers
    for line in ffile:
        if ">" in line:
            line = line.rstrip('\n')
            line = line.strip('>')
            headers.append(line)
    #print(headers)
    # prepare data frame
    info = pd.DataFrame({'contigName':headers})
    info['file'] = os.path.basename(fasta)
    info['contigs'] = info.shape[0]  

    # get bin number
    m = re.search('bin.([0-9]+).fa', fasta)
    if m:
        bin = m.group(1)
    else:
        bin = '??'
    info['bin'] = bin

    # column to match CheckM:
    m = re.search('(bin.[0-9]+).fa', fasta)
    if m:
        bin_id = m.group(1)
    else:
        bin_id = '??'
    info['Bin Id'] = bin_id  
    info['bin_id'] = info['Bin Id'].str.replace('.', '_')
    #print(info)
    return info
    

def summarise_contigs_many_files(fasta_path_list):
    fasta_files = os.listdir(fasta_path_list)
    fasta_files = [os.path.join(fasta_path_list, f) for f in fasta_files]
    #print(fasta_files)
    fasta_files = [f for f in fasta_files if '.fa' in f]
    info = pd.DataFrame()
    for fasta_file in fasta_files:
        fasta_info = summarise_contigs(fasta_file)
        info = pd.concat([info, fasta_info], axis=0)
    return info


if __name__ == '__main__':

    parser = argparse.ArgumentParser(            
                prog='summarise_contigs',         
                description="""Script analyzes FASTA headings for many bins""")
                                                 
    parser.add_argument("-input", dest="input", help='path to fasta')
    parser.add_argument("-output", dest="output", help='output path')
                                                 
    args = parser.parse_args()                   
                                                 
    if args.input is None:                       
        print(parser.print_help())               
    else:                                        
        df = summarise_contigs_many_files(args.input)
        df.to_csv(args.output, sep='\t', index=False)

    print('args.input: {}'.format(args.input))

        
