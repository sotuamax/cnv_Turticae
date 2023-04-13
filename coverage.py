""" Updated on Oct-04-2021. 
Program developed by Meiyuan Ji (for questions please contact: meiyuan.ji@utah.edu).
This program takes genome fasta file, gene annotation (GTF) and DNA-read alignment (BAM) files, and measuring the read-coverage on gene region with assigned step size. 

The purpose of this program is for providing a big-picture of read-coverage on gene region, which provide cues for estimating single-copy coverage.

Update on Mar-28-2023.
Mapping quality control with parameters setting (MQ)
Update on Jan-30-2023.
Combine bam_parse and gene_stat running in order;
go through gene exon region for read-coverage calculation (reduce running time);
Update on Jan-31-2023.
Not write temporary file in the folder. Instead, store data in a dataframe and collect all from all running cores. 
update on Feb-01-2023.
Do not collect the median and mean coverage for each gene, use the raw site coverage on gene CDS region. 
"""

import pandas as pd
import argparse
import os
import pysam
from statistics import mean
from statistics import median
from Bio import SeqIO
from mpi4py import MPI
import numpy as np 
import statistics
import random 
import time

def endtime(start):
    end = time.time()
    t = end-start
    if t < 60:
        print('{:.2f} seconds elapsed'.format(t))
    elif t < 3600:
        print('{:.2f} minutes elapsed'.format(t/60))
    else:
        print('{:.2f} hours elapsed'.format(t/3600))
    return end

def args_parser():
    '''parser the arguments from terminal command'''
    parser=argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter, description="Count read-coverage across gene coding region.")
    parser.add_argument("ref", help = "reference FASTA file. ")
    parser.add_argument("bam", help = "BAM alignment file of the reference genome (indexed and sorted). ")
    parser.add_argument("-chr", "--chromosome", help = "assign chromosome of interest. ", nargs = 1, required = False)
    parser.add_argument("-range", "--range", nargs = 2, help = "contig range for read-coverage screening. ", required = False)
    parser.add_argument("-step", "--stepsize", default = 100, type = int, help = "step size for coverage in BAM (default: 100 bp)")
    parser.add_argument("-MQ", "--mapping_quality", default = 0, type = int, help = "minimum mapping quality required (default: 0). ")
    parser.add_argument("-O", "--out", help = "output name. ", required = True)
    args = parser.parse_args()
    return args

def bam_parse(args, pos_df, seq_dict, bam_handle):
    '''Count read-coverage of given position. '''
    # pos_df is gene CDS position
    cov_list = list()
    for row in pos_df.itertuples():
        rchr = row.chromosome
        rpos = row.position
        rpos_cov = len({read.query_name for read in bam_handle.fetch(rchr, rpos, rpos+1) if seq_dict[rchr][rpos] != "N" and read.is_secondary == False and read.mapping_quality >= args.mapping_quality})
        cov_list.append(rpos_cov)
    pos_df["read_depth"] = cov_list
    return pos_df

def main():
    '''Split work based on gene number and thread number, run for read-coverage measuring via reading BAM file.'''
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    if rank == 0:
        start = time.time()
    # all threads parsing bam file
    args = args_parser()
    bam = args.bam
    out = args.out
    step = args.stepsize
    # all threads parsing reference fasta 
    ref = args.ref
    seq = SeqIO.parse(ref, "fasta")
    seq_dict = SeqIO.to_dict(seq)
    if args.chromosome != None:
        tig = args.chromosome
        range1 = int(args.range[0])
        range2 = int(args.range[1])
        seq_dict = {s:seq_dict[s] for s in seq_dict if s in tig}
    bam_handle = pysam.AlignmentFile(bam, "rb")
    # all threads parsing position dataframe
    df_list = list()
    for s in seq_dict:
        seq_len = len(seq_dict[s])
        if args.chromosome == None:
            pos_list = list(range(0, seq_len, step))
            df_s = pd.DataFrame({"chromosome":s, "position":pos_list})
        else:
            df_s = pd.DataFrame({"chromosome":s, "position":range(range1, range2)})
        print(df_s)
        df_list.append(df_s)
    df_pos = pd.concat(df_list, axis = 0)
    df_pos.index = list(range(0, len(df_pos)))
    # split work based on gene number
    if rank == 0:
        worker_tasks = {w:[] for w in range(size)} # a dictionary
        w_idx = 0
        for i in df_pos.index.tolist():
            worker_tasks[w_idx].append(i)
            w_idx = (w_idx + 1) % size
    else:
        worker_tasks = None
    worker_tasks = comm.bcast(worker_tasks, root = 0)
    # running coverage screen on BAM file
    i = worker_tasks[rank] # i stores gene names
    #
    print("Process BAM on core ", rank)
    subpos = df_pos[df_pos.index.isin(worker_tasks[rank])]
    sub_bam_count = bam_parse(args, subpos, seq_dict, bam_handle) # read-coverage of exon positions by parsing BAM file
    bam_count_gather = comm.gather(sub_bam_count, root = 0)
    if rank == 0:
        print("Gather bam count data from all processors ......")
        bam_count = pd.concat(bam_count_gather, axis = 0)
        bam_count.to_csv(out + ".txt", sep = "\t", index = False)
        endtime(start)


################
#### Run it ####
################

if __name__ == "__main__":
    main()
