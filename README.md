# Copy number variation for the two-spotted spider mites (T. urticae)

This pipeline was developed for estimating the copy nubmer variation of genes (only transcribed region of interest)
To start, you need to prepare the following files:
1. Reference genome in fasta format
2. GFF annotation for the reference (sorted and indexed)
3. DNA-seq alignment in BAM format (sorted and indexed)

### We started from reading BAM files for read-coverage information, and sites on gene region with window size of 100 bp are used for collecting. 
<code> python gene_coverage.py -ref <ref> -gtf <gtf> -bam <bam> -O dir </code>
### The median value of read-coverage per gene was used to fit a gaussian mixture model, and parameters (mean and standard deviation) are estimated based on the observations. 
<code> Rscript single_copy_coverage.R -RD mean_depth.txt -O dir </code> <br>
NOTE: the "mean_depth.txt" file is the output from the first step.
### Then, parameters was utilized for the estimation of copy-number variation of given site. 
<code> Rscript gene_CNV.R -RD pos_depth.txt -fit depth_fit.txt -O dir </code> <br>
NOTE: pos_depth.txt and depth_fit.txt are output from the first and second step, respectively. 
### To really focus on the copy variation on transcribed gene region, we filter sites based on exon-table (see my repo of gff_script for detail)
<code> python exon_CNV.py -exon <exon_tab> -CNV CNV.txt -O dir </code> <br>
NOTE: CNV.txt is output from last step. You need to make sure dir is the same for separate steps. 
  
