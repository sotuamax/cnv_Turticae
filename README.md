# Copy number variation for the two-spotted spider mites (T. urticae)

This pipeline was developed for estimating the copy nubmer variation of genes (only transcribed region of interest)
To start, you need to prepare the following files:
1. Reference genome in fasta format
2. GFF annotation for the reference (sorted and indexed)
3. DNA-seq alignment in BAM format (sorted and indexed)

### Collect read-depth information from BAM file 
Based on the reference fasta and GFF annotation file, read-depth of sites on gene region is collected. <br>
<code> python read_depth_ongene.py -ref [ref] -gtf [gtf] -bam [bam] -O [dir] </code> <br>
- If the [dir] not exists, it will create a new folder. <br>
- Under the folder, output files will be written into, including pos_depth.txt and mean_depth.txt files. <br>
- The pos_depth.txt file is read-depth for sites on gene region with default step size of 100 bp. <br>
  ![Screen Shot 2021-12-30 at 12 13 00 PM](https://user-images.githubusercontent.com/63678158/147781629-e1a2d72b-6672-4304-9a55-a79423ea243c.png)
- The depth_per_gene.txt file is calculated read-depth for each gene based on the mean and median value of depth on collected gene sites. The number of sites per gene was indicated as "size" in the output file. <br> 
  ![Screen Shot 2021-12-30 at 12 17 55 PM](https://user-images.githubusercontent.com/63678158/147781946-80e304e0-18df-4d83-945c-684ae22d5115.png)

### Estimate read-depth for single-copy region
The median value of read-coverage per gene was used to fit a gaussian mixture model, and parameters (mean and standard deviation) are estimated based on the observations. 
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
