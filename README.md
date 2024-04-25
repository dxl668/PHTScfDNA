# PHTScfDNA
[![DOI](https://zenodo.org/badge/642143482.svg)](https://zenodo.org/doi/10.5281/zenodo.10372574)

The repository contains all the bash and R scripts neecssary for data processing, analysis, and visualization associated with the paper: "Cell-free DNA Fragmentomics and Second Malignant Neoplasm Risk in Patients with PTEN Hamartomat Tumor Syndrome". 

The scripts should be ran sequentially as follows: 

## Post-alignment process of BAM files 
	(1) 01-filter_bam.sh 
 		- Apply filtering criteria (i.e., properly paired, no PCR duplicates, no secondary alignments, MAPQ >= 30) using SAMtools 'samtools view' command. 
   	(2) 02-bamtobed.sh
        	- Convert filtered BAM files to BEDPE format using BEDtools 'bamtobed' command.
	 	- Note: Modified script (02-bamtobed_par.sh) using GNU parallel to parallelize tasks. 
	(3) 03-filter_frags_gc.sh
	 	- Obtain fragment GC content using 'bedtools nuc' command.
  		- Note: The ‘bedtools nuc’ does not support BEDPE file format. Thus, to calculate the GC content, we convert the BEDPE file to a BED file, where the entire 			cfDNA fragment is represented using the start position of the first read and the end position of the second read. 
        	- Filter out cfDNA fragments from problematic regions of the genome (i.e., ENCODE blacklisted region and genomic gaps) using 'bedtools substract' command.
	(4) 04-bins5mb.sh
 		- Intersect fragments into non-overlapping 5-Mb genomic bins using 'bedtools intersect' command. 
     		- Note: The ‘bedtools intersect’ command inadvertently duplicates fragments at the boundaries of adjacent bins. To address this, the script generates two 			temporary files: one for unique fragments (i.e., those not found at bin boundaries) and de-duplicated fragments. During de-duplication, the 				duplicated fragment with the largest base pair overlap within its respective bin is retained. These two files are then concatenated.
   	(5) 05-filter_bedpe.sh
    		- Obtain paired-end alignments associated with filtered/binned fragments
      		- Note: This step is necessary as the resultant BED files no longer contain paired-end alignments that are necessary for downstream end motif analysis after 			the filtering and binning processing steps. Paired-end alignments from the original BEDPE files (i.e., output of 02-bamtobed.sh) corresponding to 			filtered/binned fragments are retrieved using the unique identifier (QNAME) from paired-end reads. 
      
 ## cfDNA Size Distribution Analysis and Genome-wide Fragmentation Profile 
 
 	(1) 10-frags_gc.sh & 10-frags_gc.R
		- Filter cfDNA fragments by size (i.e., <100 bp and >650 bp) and perform GC-bias correction. 
	(2) 11-combine_bin.sh & 11_combine_bin.R
		- Consolidate data from individual samples into single tibble for downstream analysis. 
	(3) 12-size_distr.Rmd 
		- cfDNA size distribution analysis - Figure 1A-I and Figure S1-S3. 
  	(4) 12.5-size_distr.Rmd
   		- Quantification and statistical comparison of major cfDNA peaks. 
	(4) 13-fragratio.Rmd
		- Genome-wide fragmentation profile and correlation analysis - Figure 2-3. 
  	(5) 14-mlg_loocv.RmD
   		- Multivariable logistic regression and leave-one-out cross-validation - Table2 and Figure S4. 
 
 ## cfDNA Fragment End Motif Analysis  
	(1) 06-endmotif_bed.sh
		- Generate separate BED files for each 4-mer end motif from paired-end read. 
	(2) 07-motif_fasta.sh
        	- Extract 4-mer end motif nucleotide sequencecs using 'bedtools getfasta' command. 
	(3) 20-motif_gc.sh & 20-motif_gc.R
		- Filter cfDNA fragments by size (i.e., <100 bp and >650 bp) and perform GC-bias correction. 
	(4) 21-combine_motif.sh & 21-combine_motif.R
		- Consolidate data from individual samples into single tibble for downstream analysis. 
	(5) 22-motif_analysis.Rmd
 		- Quantification, visualization, and statistical testing of 4-mer, 2-mer, and 1-mer end motifs (Figure 4A-E)
	(6) 22.5-motif_analysis.Rmd
 		- Additional script to perform overall and parwise statistical testing for 2-mer and 1-mer end motifs. 
   
   
 	
	

        
        
      
