import argparse
import pandas as pd

def checkm_to_pandas(file, output):
    file=open(file)
    checkm = list(file)
    #checkm = list(csv.reader(open(file,'rb')))
    new = []
    for list_ in checkm:
        if '-----------' in list_:
            continue
        # don't want to split "Strain heterogeneity" into multiple
        entries = list_.split('  ') # split on more than one space
        entries = [e for e in entries if len(e) > 0]
        # remove newline
        entries = [e for e in entries if e != '\n']
        entries = [e for e in entries if e != ' \n']
        new.append(entries)
    df = pd.DataFrame(new[1:], columns=new[0])
    # strip leading whitespace
    df.columns = df.columns.str.strip(r'^ ')
    df.to_csv(output, sep='\t', index=False)
    file.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
                prog='checkm_to_pandas', 
                description="""Script parses CheckM qa results to .tsv""")

    parser.add_argument("-input", dest="input", type=str, help="Specify a checkM file")
    parser.add_argument("-output", dest="output", type=str, help="output tsv filename")

    args = parser.parse_args()
    
    if args.input is None:
        print(parser.print_help())
    else:
        checkm_to_pandas(args.input, args.output)
    

