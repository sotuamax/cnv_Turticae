# Gene copy number variation (CNV) in generalist herbivore of <i>Tetranychus urticae</i> (the two-spotted spider mite) 

## Table of Contents
- [Gene CNV estimation](#Gene-CNV-estimation)
- []

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
<br>
All the scripts under "gCNV" folder

```bash
# count coverage depth on gene coding region (default stepsize 1 bp)
python gene_coverage.py [ref] [gtf] [bam] -O [out]
# ref: reference genome in fasta file <br>
# gtf: gtf file for the reference genome <br>
# bam: bam of DNA-read aligned to reference genome <br>
# out: output folder
<br>
# report the single-copy coverage depth for the BAM file
```
- File named "pos_depth.txt" (coverage at gene coding positions) will be generated under the \[out\] folder. <br>

.  <br>
<code> Rscript single_depth.R [out]/pos_depth.txt [cov_est] -O [out] </code> <br>
cov_est: the estimated coverage for the BAM file <br>
out: output folder (it should be the same as in "Step 1") <br>
- File named "single_cov.txt" (single copy coverage) will be generated under the \[out\] folder <br>
- File named "histogram.pdf" (coverage distribution) will be generated under the \[out\] folder <br>

<b>Step 3</b>: estimate gene CNV based on gene coding region coverage depth <br>
<code> python gene_CNV.py [out]/pos_depth.txt [out]/single_cov.txt -O [out] </code> <br>
out: output folder (it should be the same as in "Step 1" and "Step 2") <br>
- File named "gene_cnv.txt" will be generated under the \[out\] folder. <br>

To count coverage at a specific region of interest: <br>
Run the script "region_coverage.py" <br>
<code> python region_coverage.py [ref] [bam] -chr [chr] -range [start] [end] -step [size] -O [out] </code> <br>
chr: chromosome of interest <br>
start: start position for coverage counting <br>
end: end position for coverage counting <br>
size: step size <br>
out: output file (prefix) <br>
NOTE: this script works only for BAM file with <b>illumina-read</b> aligned against to the reference genome (not works on long-read sequencing data). 

#### To count coverage at a specific region of interest: <br>
Run the script "coverage.py" <br>
<code> mpiexec -n 10 python coverage.py [ref] [bam] -chr [chr] -range [start] [end] -step [stepsize] -O [out] </code> <br>
chr: chromosome of interest <br>
start: start position for coverage counting <br>
end: end position for coverage counting <br>
size: step size <br>
out: output file (prefix) <br>
NOTE: this script works only for BAM file with <b>illumina-read</b> aligned against to the reference genome (not works on long-read sequencing data). 
NOTE: this script works both for <b>illumina-read</b> and <b>pacbio-read</b> (long-read sequencing). You can assign multiple cores to support parallel running (e.g., use 10 cores, "mpiexec -n 10"). 

#### To recover the amino acid allele at given target-site: <br>
Run the script "target_site_allele.py" <br>
<code> python target_site_allele.py [bam] [fof] -btype <DNA/RNA> -out [out] </code> <br>
bam: bam file used for recovering <br>
fof: formated target site position <br>
btype: either be "DNA" or "RNA" <br>
out: output file name <br>
For example about how to format "fof" the target-site position file, see folder "data". <br>


