# PHTScfDNA
[![DOI](https://zenodo.org/badge/642143482.svg)](https://zenodo.org/doi/10.5281/zenodo.10372574)

The repository contains all the bash, R, and Rmd scripts necessary for data processing, analysis, and visualization for each figure of the paper 

The code should be ran as follows: 

## cfDNA Fragment Processing
	(1) 01-filter_bam.sh & 01-filter_bam.R:
		- This could have been done all in bash, but had a previous R script for filtering bam files using Rsamtools
	(2) 01.5-filter_bam.sh
        	- Filter out fragments <= 1000 bp
   	(3) 02-bamtobed.sh: 
        	- Use bedtools to convert filtered bam files to bedpe files to obtain fragment start, end, size
    	(4) 03-frags_gc.sh: 
        	- Use bedtools to filter out blacklisted regions of genome, and calculate gc content of fragments
	(5) 03.5-frags_5mb.sh: Intersect with 5 Mb bins to get both fragment interval as well as bin intervals 
        	- Note: When intersecting with 5Mb, fragments overlapping with two bins would be duplicated.
        	- Resolved duplicates by comparing base pair overlap and filtering out duplicate with the smallest overlap
 ## cfDNA Fragmentation Profile
 
 	(1) 10-fragsGC.sh & 10-fragsGC.R
		- GC correction of fragments including GC bias plots
	(2) 11-combine_bin.sh & 11_combine_bin.R
		- Combined tibble of fragments for each 5Mb genomic bin 
	(3) 12-size_distr.Rmd 
		- Perform cfDNA size distribution analysis and generate figure plots 
	(4) 13-fragratio.Rmd
		- Visualize fragment ratio, compare fragment length variability
  	(5) 14-mlg_loocv.RmD
   		- Perform multivariable logistic regression and LOOCV
     		- Note: CSV file with meta data and median fragment ratios in Zenodo repository called "14-mlg_loocv.csv"
 
 ## cfDNA End Motif Profiling 
	(1) 04-filter_bedpe.sh: 
		- Using awk, intersected/filtered fragments 
	(2) 05-end_motif_bed.sh:
		- For getfasta work, I needed to convert the bedpe file to create two separate bed files (read1, read2) with gc content of full fragment (needed for GC correction)
	(3) 06-motif_fasta.sh
        	- Using bedtools getfasta, get the nucleoptide sequence of the 5' fragment end 
	(4) 07-motif_merge.sh
        	- Merge the separate bed files for read1/read2 into a single bed file
	(5) 20-motif_gc_correct.sh & 20-motif_gc_correct.R
        	- GC correction of fragments 
        	- Note: this is identical to the GC correction from 10-fragsGC.R
	(6) 21-combine_motif.sh & 21-combine_motif.R
		- Create combined tibble of end motifs 
	(7) 22-end_motifs.Rmd
		- Summarize and visualize 4-mer, 2-mer, and 1-mer end motifs
  		- Note: 4-mer and 2-mer end motifs were calculated separately by editing the value within the two awk statements in 05-end_motif_bed.sh. 1-mer end motifs were then calculated from the 2-mer end motif frequencies in this script. 


        
        
      
