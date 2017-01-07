# https://github.com/edgraham/BinSanity/blob/4636b40595237217011130469d3e9182124ce04e/Utils/checkm_analysis

import argparse
import pandas as pd

def checkm_to_pandas(file, output):
    file=open(file)
    checkm = list(file)
    #checkm = list(csv.reader(open(file,'rb')))
    import pdb; pdb.set_trace()
    new = []
    for list_ in checkm:
        if '-----------' in list_:
            continue
        # don't want to split "Strain heterogeneity" into multiple
        entries = list_.split('  ') # split on more than one space
        entries = [e for e in entries if len(e) > 0]
        # remove newline
        entries = [e for e in entries if e != '\n']
        new.append(entries)
    df = pd.DataFrame(new)
    #df.to_csv('out.tsv', sep='\t', index=False)
    df.to_csv(output, sep='\t', index=False)

#checkm_file = open('CheckM.txt')
#df = checkm_to_pandas(checkm_file)
#df.to_csv(df_path, sep='\t', index=False)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
                prog='checkm_to_pandas', 
                usage='%(prog)s -input Chekm.txt -output parsed.tsv',
                description="""Script parses CheckM qa results to .tsv""")

    parser.add_argument("-input", dest="input", 
        help="Specify a checkM file")
    parser.add_argument("-output", dest="output", 
        type=str, default=".fna",
        help="output tsv filename")

    args = parser.parse_args()
    print(args)
    
    if args.input is None:
        print(parser.print_help())
    else:
        checkm_to_pandas(args.input, args.output)
    

