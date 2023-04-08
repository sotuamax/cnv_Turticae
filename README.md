# Copy number variation in <i>Tetranychus urticae</i> (the two-spotted spider mite)
## This is a repo for the study of gene copy number variation (CNV) of the generalist spider mite herbivore <i>Tetranychus urticae</i>. 

The pipeline was developed for estimating gene CNV that specificially focus on gene coding region.
To start, you need to prepare the following files:
1. Reference genome in fasta format
2. GFF annotation file for the reference genome (sorted and indexed)
3. illumina DNA-seq alignment in BAM format (sorted and indexed, you should have an estimation regarding to the coverage of the BAM file)

The pipeline was written in Python and R, the dependencies include:
Python v3.7: pandas, numpy, mpi4py
R v4.1

Run the following commands for gene CNV estimation:
Step 1: count coverage depth on gene coding region (default stepsize 1 bp).
<code> python gene_coverage.py [ref] [gtf] [bam] -O [out] </code> <br>
- A new "out" folder will be created (if not exist) and all output files will be written under the folder.
- File named "pos_depth.txt" will be generated in the folder.
<b>Step 2</b>: report the single-copy coverage depth for the BAM file.  
<code> Rscript single_depth.R [out]/pos_depth.txt [cov_est] -O [out] </code> <br>
cov_est: the estimated coverage for the BAM file
out: output folder (it should be the same as in "Step 1")
Step 3: generate gene CNV 
<code> python gene_CNV.py [out]/pos_depth.txt [out]/single_cov.txt -O [out] </code> <br>
- The pos_depth.txt file is read-depth for sites on gene region with default step size of 100 bp. <br>
- The depth_per_gene.txt file is calculated read-depth for each gene based on the mean and median value of depth on collected gene sites. The number of sites per gene was indicated as "size" in the output file. <br> 

### Estimate read-depth for single-copy region
Using the median value of read-depth per gene as observations to fit a gaussian mixture model, and model parameters (mean and standard deviation) are estimated from it. <br>
<code> Rscript single_copy_depth.R -RD depth_per_gene.txt -O [dir] </code> <br>
For inputs: <br>
The "depth_per_gene.txt" file is the output from the last step; <br>
The [dir] folder should keep consistent as it is in the last step. <br>
Output files: <br>
- fitting parameters are in the file "depth_fit.txt" (mean value used as the single-copy read-depth); <br>
- PDF file of observations used for modeling and where the mean-value sit in the observations. <br>

### Estimate copy number for each site
The read-depth on each site was devided by the estimated single-copy read-depth to get estimated copy number (column: norm_depth). <br>
<code> Rscript pos_CNV.R -RD [read_depth.txt] -fit [depth_fit.txt] -O [output] </code> <br>
Most of sites with estimated norm_depth of ~1 <br>


### Filter sites on transcribed gene region
For gene copy number variation, we are interested in the variation happened on exon region. To achieve it, sites outside of gene exon are filtered out. <br>
To prepare exon_table for this step, see my repo (gff_script) for detail. 
<code> Rscript gene_CNV.R -RD pos_depth.txt -fit depth_fit.txt -O dir </code> <br>
NOTE: pos_depth.txt and depth_fit.txt are output from the first and second step, respectively. 

<code> python exon_CNV.py -exon [exon_tab] -CNV [CNV.txt] -O dir </code> <br>
NOTE: CNV.txt is output from last step. You need to make sure dir is the same for separate steps. 

### Estimate copy number variation on gene-basis

<code> python CNV_estimate.py -CNV <exon.txt> -O <output> </code> <br>
The default cutoff of padj is 0.05. When more than 65% sites for each gene show single-copy, gene was assigned as single-copy number. <br>  
