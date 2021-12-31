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
Using the median value of read-depth per gene as observations to fit a gaussian mixture model, and model parameters (mean and standard deviation) are estimated from it. <br>
<code> Rscript single_copy_depth.R -RD depth_per_gene.txt -O [dir] </code> <br>
For inputs: <br>
The "depth_per_gene.txt" file is the output from the last step; <br>
The [dir] folder should keep consistent as it is in the last step. <br>
Output files: <br>
- fitting parameters are in the file "depth_fit.txt" (mean value used as the single-copy read-depth); <br>
- PDF file of observations used for modeling and where the mean-value sit in the observations. <br>
![Screen Shot 2021-12-30 at 1 17 19 PM](https://user-images.githubusercontent.com/63678158/147785488-522c60a3-8e8d-4350-b2f6-38e32fdc4a08.png)

### Estimate copy number for each site
The read-depth on each site was devided by the estimated single-copy read-depth to get estimated copy number (column: norm_depth). <br>
<code> Rscript pos_CNV.R -RD [read_depth.txt] -fit [depth_fit.txt] -O [output] </code> <br>
Most of sites with estimated norm_depth of ~1 <br>
![Screen Shot 2021-12-30 at 1 36 12 PM](https://user-images.githubusercontent.com/63678158/147786442-a0e4d48c-ad56-4e4f-8128-bf3b6f363b2f.png)

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
