# Gene copy number variation (CNV) in generalist herbivore of <i>Tetranychus urticae</i> (the two-spotted spider mite) 

## Table of Contents
- [Gene CNV estimation](#Gene-CNV-estimation)
- [Sequencing depth for genomic regions of interest](Sequencing-depth-for-genomic-regions-of-interest)
- [Amino acid at target-site](Amino-acid-at-target-site)
- [Genome scaffolding](Genome-scaffolding)

## Programs used / Dependencies 
- python3+ (modules: pysam, pandas, numpy, Biopython, mpi4py)
- R v4.1
- 

[NOTE] To enable parallel processing, python model mpi4py need to be installed. 

## Gene CNV estimation 
This pipeline was developed for estimating CNV that focus on gene coding regions. 
<br>
Input data: <br>
1. Reference genome in fasta format <br>
2. GTF annotation file for the reference genome (sorted and indexed) <br>
3. Illumina DNA-seq alignment in BAM format (sorted and indexed) <br>

All the scripts under "gCNV" folder: 
<br>
Step 1: count coverage on gene coding region (default stepsize 1 bp)
```bash
mpiexec -n 10 python gene_coverage.py [ref] [gtf] [bam] -step 5 -O [out]
# -n: core number 
# ref: reference genome in fasta file <br>
# gtf: gtf file for the reference genome <br>
# bam: bam of DNA-read aligned to reference genome <br>
# -step: stepsize assigned as 5 bp
# out: output folder
```
Step 2: estimate single-copy coverage depth for the BAM file
```bash
Rscript single_depth.R [out]/pos_depth.txt [cov_est] -O [out]
# cov_est: the start estimation value for BAM file coverage  <br>
```

Step 3: estimate gene CNV based on gene coding region coverage depth
```bash
mpiexec -n 10 python gene_CNV.py [out]/pos_depth.txt [out]/single_cov.txt -O [out]
# -n: core number 
```
Outputs: <br>
pos_depth.txt: raw data for coverage at gene coding positions <br>
single_cov.txt: single copy coverage value <br>
histogram.pdf: coverage distribution using a histogram <br>
<b>gene_cnv.txt: gene CNV report file</b> 
<br>
<br>
To run all three steps together: 
```bash
mpiexec -n 10 python gene_coverage.py [ref] [gtf] [bam] -step 5 -O [out] && Rscript single_depth.R [out]/pos_depth.txt [cov_est] -O [out] && mpiexec -n 10 python gene_CNV.py [out]/pos_depth.txt [out]/single_cov.txt -O [out]
```

## Sequencing coverage for genomic regions of interest 
To count coverage at a specific regions of interest: <br>
Input data: <br>
1. Reference genome in fasta format <br>
2. Illumina DNA-seq alignment in BAM format (sorted and indexed) <br>
<br>
Run the script "region_coverage.py": <br>
```bash
python region_coverage.py [ref] [bam] -chr [chr] -range [start] [end] -step [size] -O [out]
# ref: reference genome in fasta file <br>
# bam: bam of DNA-read aligned to reference genome <br>
# chr: chromosome of interest <br>
# start: start position for coverage counting <br>
# end: end position for coverage counting <br>
# size: step size (default 100 bp) <br>
# out: output file (prefix) <br>
```

## Genome scaffolding
Scaffold contigs to pseudo-chromosome level assembly using reference based approach 
```bash

```

#### Amino acid at target-site 

```bash
python target_site_allele.py [bam] [fof] -btype <DNA/RNA> -out [out]
# bam: bam file used for recovering <br>
# fof: formated target site position <br>
# btype: either be "DNA" or "RNA" <br>
# out: output file name <br>

```
For example about how to format "fof" the target-site position file, see folder "data"





