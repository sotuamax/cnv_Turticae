# Copy number variation for the two-spotted spider mites (T. urticae)

This pipeline was developed for estimating the copy nubmer variation of genes (only transcribed region of interest)
To start, you need to prepare the following files:
1. Reference genome in fasta format
2. GFF annotation for the reference (sorted and indexed)
3. DNA-seq alignment in BAM format (sorted and indexed)

### We started from reading BAM files for read-coverage information, and sites on gene region with window size of 100 bp are used for collecting. 
<code> python gene_coverage.py -ref <ref> -gtf <gtf> -bam <bam> -O <dir> </code>
- If the <dir> not exists, it will create a new folder. Under the folder, output files will be written into, including pos_depth.txt and mean_depth.txt files. The pos_depth.txt file is read-coverage for sites on gene region with step size of 100 bp. [Example]![Screen Shot 2021-12-30 at 12 13 00 PM](https://user-images.githubusercontent.com/63678158/147781629-e1a2d72b-6672-4304-9a55-a79423ea243c.png)

### The median value of read-coverage per gene was used to fit a gaussian mixture model, and parameters (mean and standard deviation) are estimated based on the observations. 
<code> Rscript single_copy_coverage.R -RD mean_depth.txt -O dir </code> <br>
NOTE: the "mean_depth.txt" file is the output from the first step.
### Then, parameters was utilized for the estimation of copy-number variation of given site. 
<code> Rscript gene_CNV.R -RD pos_depth.txt -fit depth_fit.txt -O dir </code> <br>
NOTE: pos_depth.txt and depth_fit.txt are output from the first and second step, respectively. 
### To really focus on the copy variation on transcribed gene region, we filter sites based on exon-table (see my repo of gff_script for detail)
<code> python exon_CNV.py -exon <exon_tab> -CNV CNV.txt -O dir </code> <br>
NOTE: CNV.txt is output from last step. You need to make sure dir is the same for separate steps. 
### Based on the normalized read-coverage on gene exon region, copy number variation per gene was estimated. 
<code> python CNV_estimate.py -CNV <exon.txt> -O <output> </code> <br>
The default cutoff of padj is 0.05. When more than 65% sites for each gene show single-copy, gene was assigned as single-copy number. <br>  
