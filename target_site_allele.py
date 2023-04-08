import os
import pysam
import argparse
from Bio import SeqIO
import pandas as pd
import numpy as np
from collections import Counter

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="Usage: ")
    parser.add_argument("bam", help = "bam files. ")
    parser.add_argument("fof", help = "file of filename with target site for checking (always in triplet unit, 1-based). ")
    parser.add_argument("-btype", "--btype", help = "assign as DNA-seq or RNA-seq BAM alignment type")
    parser.add_argument("-out", "--output", help = "output prefix")
    args = parser.parse_args()
    return args

def bam_allele_parse_RNA(args):
    fof = args.fof
    bam = args.bam
    out = args.output
    df = pd.read_table(fof, sep = "\t", header = 0, comment = "#")
    bam_handle = pysam.AlignmentFile(bam, "rb")
    fw = open(out + "_allele.txt", "w")
    fw.write("read_name\ttig\tposition\tallele\n")
    ###
    for row in df.itertuples():
        rname = row.name
        rtig = row.tig
        rstart = row.genome_start - 1 # make it 0-based for bam retrieving
        rend = row.genome_end
        rstrand = row.strand
        range_list = list(range(rstart, rend))
        for read_column in bam_handle.pileup(rtig, rstart, rend):
            if read_column.reference_pos in range_list:
                for bam_read in read_column.pileups:
                    # uniquely mapped
                    # no deletion
                    if bam_read.query_position != None:
                        read_name = bam_read.alignment.query_name
                        read_allele = bam_read.alignment.query_sequence[bam_read.query_position]
                        ref_pos = read_column.reference_pos
                        fw.write(f"{read_name}\t{rtig}\t{ref_pos+1}\t{read_allele}\n") # 1-based position
    fw.close()

def bam_allele_parse_DNA(args):
    fof = args.fof
    bam = args.bam
    out = args.output
    df = pd.read_table(fof, sep = "\t", header = 0, comment = "#")
    bam_handle = pysam.AlignmentFile(bam, "rb")
    fw = open(out + "_allele.txt", "w")
    fw.write("read_name\ttig\tposition\tallele\n")
    ###
    for row in df.itertuples():
        rname = row.name
        rtig = row.tig
        rstart = row.genome_start - 1 # make it 0-based for bam retrieving
        rend = row.genome_end
        range_list = list(range(rstart, rend))
        for read_column in bam_handle.pileup(rtig, rstart, rend, min_base_quality = 20, min_mapping_quality = 10):
            if read_column.reference_pos in range_list:
                for bam_read in read_column.pileups:
                    if bam_read.query_position != None:
                        read_name = bam_read.alignment.query_name
                        read_allele = bam_read.alignment.query_sequence[bam_read.query_position]
                        ref_pos = read_column.reference_pos
                        fw.write(f"{read_name}\t{rtig}\t{ref_pos+1}\t{read_allele}\n") # 1-based position
    fw.close()

def read_call(args, out):
    fof = args.fof
    df = pd.read_table(fof, sep = "\t", header = 0, comment = "#")
    read_df = pd.read_table(out + "_allele.txt", sep = "\t", header = 0)
    # remove any duplicate rows (overlapped paired reads)
    read_df.drop_duplicates(keep = "first", inplace = True)
    fw = open(out + "_code.txt", "w")
    fw.write("name\tread_name\tcode\taa_code\n")
    # still, translate on individual row
    for row in df.itertuples():
        rname = row.name
        rtig = row.tig
        rstart = row.genome_start # still 1-based to get the range of 1-based position
        rend = row.genome_end
        range_list = list(range(rstart, rend+1))
        rstrand = row.strand
        sub_read = read_df[(read_df["tig"] == rtig) & (read_df["position"] >= rstart) & (read_df["position"] <= rend)]
        sub_read_name = sorted(set(sub_read["read_name"]))
        for rn in sub_read_name:
            rn_sub = sub_read[sub_read["read_name"] == rn]
            rn_sub = rn_sub[["position", "allele"]].sort_values(by = "position")
            if len(rn_sub) == len(range_list):
                rn_code = rn_sub.allele.str.cat()
                if rstrand == "-":
                    rn_code = reverse(complement(rn_code))
                aa_code = table_na[rn_code]
                fw.write(f"{rname}\t{rn}\t{rn_code}\t{aa_code}\n")
    fw.close()

def reverse(seq):
    """Returns a reversed string"""
    return seq[::-1]

def complement(seq):
    """Returns a complement DNA sequence"""
    complement_dict = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'N': 'N'}
    seq_list = list(seq)
    seq_list = [complement_dict[base] for base in seq_list]
    return ''.join(seq_list)

def code_prop(out):
    code_df = pd.read_table(out + "_code.txt", sep = "\t", header = 0)
    cname = sorted(set(code_df["name"]))
    fw = open(out + "_ratio.txt", "w")
    fw.write("name\taa_code\tfreq\n")
    for c in cname:
        sub_code = code_df[code_df["name"] == c]
        sub_code_dict = Counter(list(sub_code["aa_code"]))
        for ss in sub_code_dict:
            ss_freq = sub_code_dict[ss]
            fw.write(f"{c}\t{ss}\t{ss_freq}\n")
    fw.close()

def main():
    args = args_parser()
    if args.btype == "DNA":
        bam_allele_parse_DNA(args)
    if args.btype == "RNA":
        bam_allele_parse_RNA(args)
    if not args.btype in ["DNA", "RNA"]:
        print("BAM type is not correct! ")
        exit()
    out = args.output
    read_call(args, out)
    code_prop(out)
    
table_na = { 
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R', 
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W', 
}

##############
### Run it ###
##############

if __name__ == "__main__":
    main()

