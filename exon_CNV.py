import pandas as pd
import argparse
import os

def args_parser():
    '''parser the arguments from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description=". \n\
        Usage: ")
    parser.add_argument("-exon", "--exon", help = "transcribed gene region (parsed from GFF3)")
    parser.add_argument("-CNV", "-CNV", help = "CNV for sites on gene region. ")
    parser.add_argument("-O", "--outdir", help = "output folder.")
    args = parser.parse_args()
    return args

def exon_number(exon):
    exon_df = pd.read_table(exon, sep = "\t", header = 0)
    contigs = sorted(set(exon_df["contig"]))
    exon_site = dict()
    for tig in contigs:
        sub_exon = exon_df[exon_df["contig"] == tig]
        exon_site[tig] = list()
        for row in sub_exon.itertuples():
            start = row.start
            end = row.end
            range_list = list(range(start, end))
            exon_site[tig] += range_list
    return exon_site

def CNV_exon(CNV, exon_dict):
    CNV_df = pd.read_table(CNV, sep = "\t", header = 0)
    contigs = sorted(set(CNV_df["chromosome"]))
    CNV_new = list()
    for c in contigs:
        sub_CNV = CNV_df[CNV_df["chromosome"] == c]
        sub_CNV_new = sub_CNV[sub_CNV.position.isin(exon_dict[c])]
        CNV_new.append(sub_CNV_new)
    CNV_df_new = pd.concat(CNV_new, axis = 0)
    return CNV_df_new

def main():
    args = args_parser()
    exon = args.exon
    CNV = args.CNV
    exon_dict = exon_number(exon)
    CNV_exon_df = CNV_exon(CNV, exon_dict)
    CNV_exon_df.to_csv(os.path.join(args.outdir, "CNV_exon.txt"), sep = "\t", index = False)

################
#### Run it ####
################

if __name__ == "__main__":
    main()

