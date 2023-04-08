# Copy number variation in <i>Tetranychus urticae</i> (the two-spotted spider mite)
## This is a repo for the study of gene copy number variation (CNV) of the generalist spider mite herbivore <i>Tetranychus urticae</i>. 

#### The pipeline was developed for estimating gene CNV that specificially focus on gene coding region.
To start, you need to prepare the following files: <br>
1. Reference genome in fasta format <br>
2. GFF annotation file for the reference genome (sorted and indexed) <br>
3. illumina DNA-seq alignment in BAM format (sorted and indexed, you should have an estimation regarding to the coverage of the BAM file) <br>

#### The pipeline was written in Python and R, the dependencies include: <br>
Python v3.7: pandas, numpy, mpi4py <br>
R v4.1 <br>

#### Run the following commands for gene CNV estimation: <br>
All the scripts under "gCNV" folder. <br>
<b>Step 1</b>: count coverage depth on gene coding region (default stepsize 1 bp). <br>
<code> python gene_coverage.py [ref] [gtf] [bam] -O [out] </code> <br>
ref: reference genome in fasta file <br>
gtf: gtf file for the reference genome <br>
bam: bam of reads aligned to reference genome <br>
out: output folder
- A new "out" folder will be created (if not exist) and all output files will be written under the folder. <br>
- File named "pos_depth.txt" (coverage at gene coding positions) will be generated under the \[out\] folder. <br>

<b>Step 2</b>: report the single-copy coverage depth for the BAM file.  <br>
<code> Rscript single_depth.R [out]/pos_depth.txt [cov_est] -O [out] </code> <br>
cov_est: the estimated coverage for the BAM file <br>
out: output folder (it should be the same as in "Step 1") <br>
- File named "single_cov.txt" (single copy coverage) will be generated under the \[out\] folder <br>
- File named "histogram.pdf" (coverage distribution) will be generated under the \[out\] folder <br>

<b>Step 3</b>: estimate gene CNV based on gene coding region coverage depth <br>
<code> python gene_CNV.py [out]/pos_depth.txt [out]/single_cov.txt -O [out] </code> <br>
out: output folder (it should be the same as in "Step 1" and "Step 2") <br>
- File named "gene_cnv.txt" will be generated under the \[out\] folder. <br>

#### To count coverage at a specific region of interest: <br>
Run the script "region_coverage.py" <br>
<code> python region_coverage.py [ref] [bam] -chr [chr] -range [start] [end] -step [size] -O [out] </code> <br>
chr: chromosome of interest <br>
start: start position for coverage counting <br>
end: end position for coverage counting <br>
size: step size <br>
out: output file (prefix) <br>

#### 

