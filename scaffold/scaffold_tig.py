
from Bio import SeqIO
import itertools
import pandas as pd
import argparse
import statistics
import collections

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="Takes the assigned contigs file, mis-assembled contigs file, and file to break the mis-assembled contigs, return ordered and orientated scaffold contigs (*.fa) given the source contig file. contig order and orientation were recorded in *.log file. \n\
        Usage: python tig_order.py -asn <asn> -contig_fa <fa>  -O <prefix> ")
    parser.add_argument("-asn", "--asn", required = True, help = "contig assignment on chromosome with strand. ")
    parser.add_argument("-contig_fa", "--contig_fa", help = "fasta file for the contigs. label contain sequence length")
    parser.add_argument("-O", "--output", help="output prefix (generate one scaffold fasta file). ")
    args = parser.parse_args()
    return args

def build_seqdict(args):
    output = args.output
    ## contig information acquire
    contig = args.contig_fa
    seq = SeqIO.parse(contig, "fasta")
    seq_dict = SeqIO.to_dict(seq)
    #
    asn = args.asn
    asn = pd.read_table(asn, sep = "\t", header = 0)
    asn = asn.sort_values(by = ["chromosome", "tig_median"], ascending = True) # make sure chromosome in order
    # 
    tig_seq_dict = dict()
    # collect contigs sequences for scaffolding (reverse complementing contigs if negative strand)
    for row in asn.itertuples():
        t = row.contig
        std = row.strand
        tseq = str(seq_dict[t].seq)
        if std == "+":
            tig_seq_dict[t] = tseq
        else:
            tig_seq_dict[t] = complement(reverse(tseq))
    return seq_dict,tig_seq_dict

def nogene_tigs(seq_dict, tig_seq_dict):
    tig_no = [seq_dict[entry].id for entry in seq_dict if seq_dict[entry].id not in tig_seq_dict]
    return tig_no

def order_contigs(args, seq_dict, tig_seq_dict, tig_no):
    #
    output = args.output
    fw = open(output+".fa", "w")
    #
    asn = args.asn
    asn = pd.read_table(asn, sep = "\t", header = 0)
    chrs = sorted(set(asn["chromosome"]))
    #
    N = itertools.repeat("N", 100)
    N = "".join(list(N)) ### N gap
    for c in chrs:
        subdf = asn[asn["chromosome"] == c]
        psum = sum(subdf["pgene"])
        chr_tigs = [row.contig for row in subdf.itertuples()]
        tigg_n = len(chr_tigs)
        chr_scaffold = f"{N}".join([tig_seq_dict[t] for t in chr_tigs])
        chr_len = len(chr_scaffold)
        fw.write(f">{c} len={chr_len} pgene={psum} tigs={tigg_n}\n{chr_scaffold}\n")
    for t in tig_no:
        seqt = str(seq_dict[t].seq)
        seql = len(seqt)
        fw.write(f">{t} len={seql}\n{seqt}\n")
    fw.close()

def no_tig_write(args, seq_dict, tig_no):
    fw = open(args.output + "_un-asn.fa", "w")
    for t in tig_no:
        tid = seq_dict[t].id
        tseq = seq_dict[t].seq
        fw.write(f">{tid}\n{tseq}\n")
    fw.close()

def reverse(seq):
    """Returns a reversed string"""
    return seq[::-1]

def complement(seq):
    """Returns a complement DNA sequence"""
    complement_dict = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'N':'N'}
    seq_list = list(seq)
    seq_list = [complement_dict[base] for base in seq_list]
    return ''.join(seq_list)

def main():
    args = args_parser()
    # from contig assignment to contig dictionary
    seq_dict,tig_dict = build_seqdict(args)
    # get contigs with no genes on it
    tig_no = nogene_tigs(seq_dict, tig_dict)
    no_tig_write(args, seq_dict, tig_no)
    # concatenate contigs based on its order
    order_contigs(args, seq_dict, tig_dict, tig_no)

##############
### Run it ###
##############

if __name__ == "__main__":
    main()
