import pandas as pd
import numpy as np
import argparse
import os

def args_parser():
    '''parser the arguments from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description=". \n\
        Usage: ")
    parser.add_argument("-CNV", "--CNV", help = "CNV for sites on gene region. ")
    parser.add_argument("-cutoff", "--cutoff", default = 0.05, help = "padj cutoff value. ")
    parser.add_argument("-O", "--output", help = "output prefix. ")
    args = parser.parse_args()
    return args

def cnv_parse(args):
    cnv = args.CNV
    cutoff = args.cutoff
    cnv_df = pd.read_table(cnv, sep = "\t", header = 0)
    cnv_df["CNV"] = np.where(cnv_df.padj > cutoff, 1, np.where((cnv_df.padj < cutoff) & (cnv_df.norm_depth > 0.1) & (cnv_df.norm_depth < 1), 1, np.where((cnv_df.padj < cutoff) & (cnv_df.norm_depth < 0.1), 0, cnv_df.norm_depth.round()))) # ***
    cnv_gene = cnv_df[["gene", "CNV"]]
    cnv_sum = cnv_gene.groupby(by = "gene").sum().round(2)
    cnv_sum = cnv_sum.rename(columns = {"CNV":"CNV_sum"})
    cnv_size = cnv_gene.groupby(by = "gene").count()
    cnv_size = cnv_size.rename(columns = {"CNV":"site_sum"})
    ###
    single_site = {i[0]:len(i[1][i[1].CNV == 1]) for i in cnv_gene.groupby(by = "gene")}
    ###
    gene_info = pd.concat([cnv_sum, cnv_size], axis = 1)
    gene_info["single_site"] = gene_info.index.map(single_site)
    gene_info["gene_CNV"] = np.where(gene_info.single_site/gene_info.site_sum > 0.6, 1, gene_info.CNV_sum/gene_info.site_sum).round(2)
    ###
    return gene_info

def main():
    args = args_parser()
    output = args.output
    gene_cnv = cnv_parse(args)
    gene_cnv.to_csv(output+".txt", sep = "\t", index = True)

################
#### Run it ####
################

if __name__ == "__main__":
    main()


