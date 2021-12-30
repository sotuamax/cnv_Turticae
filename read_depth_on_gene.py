""" Updated on Oct-04-2021. 
Developed by Meiyuan Ji (for questions please contact: meiyuan.ji@utah.edu).
This program takes genome fasta file, gene annotation (GTF) and DNA-read alignment (BAM) files, and measuring the read-coverage on gene region with assigned step size. 

The purpose of this program is for providing a big-picture of read-coverage on gene region, which provide hints for estimating single-copy number of read-coverage genome wide.
"""


import pandas as pd
import argparse
import os
import pysam
from Bio import SeqIO
from mpi4py import MPI

def args_parser():
    '''parser the arguments from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description=". \n\
        Usage: python gene_coverage.py -ref <ref> -gtf <gtf> -bam <bam> -O <outputdir>")
    parser.add_argument("-ref", "--reference", help = "reference fasta file. ")
    parser.add_argument("-gtf", "--gtf", help = "gtf file matched to the reference genome (sorted and compressed). ")
    parser.add_argument("-pos", "--position", help = "gene position table. ", required = False)
    parser.add_argument("-bam", "--bam", help = "the bam file for read counts (sorted). ")
    parser.add_argument("-contig", "--contig", help = "assign contig of interest. ", nargs = "?", required = False)
    parser.add_argument("-step", "--stepsize", default = 100, help = "step size (default 100 bp)")
    parser.add_argument("-O", "--outdir", help = "output folder (create a new one when assigned not exists). ")
    args = parser.parse_args()
    return args

def GTF_pos(gtf):
    '''Takes gtf file and write gene position in a tab-seperated table (contig, start, end). '''
    gtf_handle = pysam.TabixFile(gtf, parser = pysam.asGTF())
    gene_pos = open("gene_pos.txt", "w")
    gene_pos.write(f"gene\tcontig\tstart\tend\n")
    for gene in gtf_handle.fetch():
        if gene.feature == "gene":
            attr = gene.attributes
            gid = attr.split(";")[0].split(":")[-1]
            tig = gene.contig
            start = gene.start
            end = gene.end
            gene_pos.write(f"{gid}\t{tig}\t{start}\t{end}\n")
    gene_pos.close()
    gene_pos = pd.read_table("gene_pos.txt", sep = "\t", header = 0)
    return gene_pos

def gene_pos_fun(args):
    '''Takes contigs of interest, and return genes on interested contigs as a position table'''
    contig = args.contig # if provided, it should be 1 or more contigs
    pos = args.position
    gtf = args.gtf
    if pos == None:
        gene_pos = GTF_pos(gtf)
    else:
        gene_pos = pd.read_table(pos, sep = "\t", header = 0, index_col = None)
    if contig != None:
        selected_contig = list(contig)
        gene_pos = gene_pos[gene_pos["contig"].isin(selected_contig)]
    return gene_pos

def bam_parse(args, pos_df, out, seq_dict, bam_handle): # ****
    '''takes gene position table and bam file, measure the read-coverage on given position with known step-size. '''
    # out is rank named
    step = args.stepsize
    outdir = args.outdir
    fw = open(os.path.join(outdir, out+".txt"), "w")
    fw.write("gene\tchromosome\tposition\tread_depth\n")
    for row in pos_df.itertuples():
        gid = row.gene # gene label
        chr = row.contig # contig
        s = row.start # gene start 0-based
        e = row.end # gene end
        # when reference position is "N", read-count will be skipped; 
        # when read mapping quality < 20, read will not be counted.
        for p in range(s, e+1, step):
            if (seq_dict[chr][p] != "N"):
                total_read = set()
                for read_column in bam_handle.pileup(chr, p, p+1, min_mapping_quality = 20): # parse the position
                    if read_column.reference_pos == p:
                        for bam_read in read_column.pileups:
                            read_name = bam_read.alignment.query_name
                            total_read.add(read_name)
                total_read_count = len(total_read)
                fw.write(f"{gid}\t{chr}\t{p}\t{total_read_count}\n")
    fw.close()
    #
    df = pd.read_table(os.path.join(outdir, out+".txt"), sep = "\t", header = 0)
    return df

def mean_median_gene(args, df):
    '''Caculate mean and median of coverage of each gene based on the collected sites' coverage. '''
    outdir = args.outdir
    # fw = open(os.path.join(outdir, out+".txt"), "w")
    df_gene = df[["gene", "read_depth"]]
    gene_mean = df_gene.groupby(by = "gene").mean().round(2)
    gene_median = df_gene.groupby(by = "gene").median().round(1)
    gene_size = df_gene.groupby(by = "gene").size()
    ###
    gene_mean = gene_mean.rename(columns = {"read_depth":"depth_mean"})
    gene_median = gene_median.rename(columns = {"read_depth":"depth_median"})
    gene_stat = pd.concat([gene_mean, gene_median], axis = 1)
    gene_stat["size"] = gene_size
    gene_stat.to_csv(os.path.join(outdir, "depth_per_gene.txt"), sep = "\t", index = True)

def main():
    '''Split work based on gene number and thread number, run for read-coverage measuring via reading BAM file.'''
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    # all threads parsing bam file
    args = args_parser()
    bam = args.bam
    outdir = args.outdir
    if outdir not in os.listdir():
        os.mkdir(outdir)
    bam_handle = pysam.AlignmentFile(bam, "rb")
    # all threads parsing reference fasta 
    ref = args.reference
    seq = SeqIO.parse(ref, "fasta")
    seq_dict = SeqIO.to_dict(seq)
    # all threads parsing position dataframe
    pos_df = gene_pos_fun(args)
    if rank == 0:
        print("Gene position file ready ......")
    pos_df.index = range(len(pos_df))
    # split work based on gene number (not gene size)
    if rank == 0:
        gid_num = len(pos_df)
        genes = sorted(list(pos_df["gene"]))
        worker_tasks = {w:[] for w in range(size)} # a dictionary
        w_idx = 0
        for i in range(gid_num):
            worker_tasks[w_idx].append(genes[i])
            w_idx = (w_idx + 1) % size
    else:
        worker_tasks = None
    worker_tasks = comm.bcast(worker_tasks, root = 0)
    # running job
    i = worker_tasks[rank] # i stores gene names
    subpos = pos_df[pos_df.gene.isin(i)]
    print("Process BAM on ", rank)
    subdf = bam_parse(args, subpos, "rank"+str(rank), seq_dict, bam_handle)
    # 
    df = comm.gather(subdf, root = 0) # in a list
    if rank == 0:
        df = pd.concat(df, axis = 0)
        df = df.sort_values(by = ["chromosome", "position"])
        df.to_csv(os.path.join(outdir, "pos_depth.txt"), sep = "\t", index = False)
    else:
        df = None
    df = comm.bcast(df, root = 0)
    if rank == 0:
        mean_median_gene(args, df)
    # 
    if rank == 0:
        os.remove(os.path.join(outdir, "rank" + str(rank) + ".txt"))

################
#### Run it ####
################

if __name__ == "__main__":
    main()
