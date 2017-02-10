import argparse
import glob
import os
import sys


def prep_filename(gff_dir):
    # strip off trailing '/'
    gff_dir = gff_dir.rstrip('/')
    return os.path.join(gff_dir, '{}_gffs_concatenated.gff'.format(gff_dir))    

def combine_gffs(gff_file_list, filename):
    # sort gff files, just in case
    gff_file_list = sorted(gff_file_list)

    with open(filename, 'w') as outfile:
        file_num = 1
        for f in gff_file_list:
            print('get the good stuff from {}'.format(f))
            with open(f) as f1:
                for line_num, line in enumerate(f1):        #keep the header from file1
                    # The first line is `##gff-version 3`.  Keep that for the first file. 
                    if (file_num == 1) and (line_num == 0):
                        outfile.write(line)
                    # The end of the file has the entire FASTA sequence glued on.  Remove it.
                    elif '##FASTA' in line:
                        break
                    # Delete subsequent lines like `##sequence-region k141_461591 1 2140`
                    elif (line_num > 0) and line.startswith('##'):
                        continue
                    # Don't include the FASTA sequences at the bottom. 
                    outfile.write(line)
                file_num += 1

if __name__ == '__main__':
    # expect python 3 to be used.  Source activate py3
    assert sys.version_info > (3, 0), 'expected python 3; got {}'.format(sys.version)    

    parser = argparse.ArgumentParser()
    parser.add_argument('gff_parent_folder', help='parent folder to get gffs from')
    args = parser.parse_args()
    
    if args.gff_parent_folder is None:
        print(args.print_help())
    else:
        print('gff parent folder: {}'.format(args.gff_parent_folder))
    
    gff_file_list = []
    search_path = args.gff_parent_folder + '/**/*.gff'
    print('look for all gff files at paths like {}'.format(search_path))
    for filename in glob.glob(search_path):
        gff_file_list.append(filename)
    print('found {} gff files: \n{}'.format(len(gff_file_list), gff_file_list))

    filename = prep_filename(args.gff_parent_folder)
    print('save concatenated gffs to {}'.format(filename))
    
    # combine the found gff files
    combine_gffs(gff_file_list, filename)

