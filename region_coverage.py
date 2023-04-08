""" Updated on Dec-29-2021. 
Developed by Meiyuan Ji (for questions please contact: meiyuan.ji@utah.edu).
# modified on Jan-20-2022. 
# random pick position within some region range and each picking window size of step size (default: 100 bp). 
"""

# v1: add truncate = True
# v2: add flag_filter = 2308
# v3: change flag_filter = 256

import pandas as pd
import argparse
import os
import pysam
from Bio import SeqIO
import random 

def args_parser():
    '''parser the arguments from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description=". \n\
        Usage: python gene_coverage.py -ref <ref> -gtf <gtf> -bam <bam> -O <outputdir>")
    parser.add_argument("ref", help = "reference fasta file. ")
    parser.add_argument("bam", help = "the bam file for read counts (sorted). ")
    parser.add_argument("-chr", "--chromosome", help = "assign chromosome of interest. ", nargs = "?", required = False)
    parser.add_argument("-range", "--range", nargs = 2, help = "contig range for read-coverage screening. ")
    parser.add_argument("-step", "--stepsize", type = int, default = 100, help = "step size (default 100 bp)")
    parser.add_argument("-quality", "--quality", type = int, required = False, default = 0, help = "mapping quality (default: 0)")
    parser.add_argument("-O", "--output", help = "output prefix. ")
    args = parser.parse_args()
    return args

def bam_coverage(args, seq_dict):
    '''takes gene position table and bam file, measure the read-coverage on given position with known step-size. '''
    # out is rank named
    bam = args.bam
    tig = args.contig
    range1 = int(args.range[0])
    range2 = int(args.range[1])
    step = args.stepsize
    mq = args.quality
    output = args.output
    ###
    bam_handle = pysam.AlignmentFile(bam, "rb")
    fw = open(output+".txt", "w")
    fw.write("chromosome\tposition\tread_depth\n")
    ###
    for p in range(range1, range2+step, step):
        pos = random.randint(p, p+step-1) 
        if (seq_dict[tig][pos] != "N"):
            total_read = set()
            for read_column in bam_handle.pileup(tig, pos, pos+1, min_mapping_quality = mq, ignore_orphans = False, truncate = True, flag_filter = 256): # parse the position
                # truncate (bool) set True, only columns in the exact region specified are returned
                # stepper (string) "all" to skip reads in which any of the following are set (unmapped, secondary, qc failure, PCR duplication)
                # ignore_orphans (bool) set False, reads not properly paired are counted
                # min_base_quality (int) default 13 for base-quality
                # min_mapping_quality (int) default 0 for the minimum mapping quality (quality of 10, equals to expected error of 1 in 10)
                # flag_filter (int) ignore reads where any of the bits in the flag are set (2308 for unmapped, not primary alignment, supplementary alignment)
                # flag_filter (int) ignore reads where any of the bits in the flag are set (2308 for unmapped, not primary alignment) 
                # supplementary alignment should be included as chimeric reads often observed at the junction region
                if read_column.reference_pos == pos:
                    for bam_read in read_column.pileups:
                        read_name = bam_read.alignment.query_name
                        total_read.add(read_name)
            total_read_count = len(total_read)
            fw.write(f"{tig}\t{pos}\t{total_read_count}\n")
    fw.close()

def main():
    args = args_parser()
    ref = args.reference
    seq = SeqIO.parse(ref, "fasta")
    seq_dict = SeqIO.to_dict(seq)
    ###
    bam_coverage(args, seq_dict)

################
#### Run it ####
################

if __name__ == "__main__":
    main()
