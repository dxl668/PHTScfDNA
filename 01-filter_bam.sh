#!/bin/bash

#load R
module load R

#Define working directoy
CWD=/home/liud3/beegfs/cfDNA/data/cfdna_pipeline
cd $CWD

#Define directory for bam/bai files
bamdir=/home/liud3/beegfs/cfDNA/data/files

#Make directory for output of 01-filter_bam.R script 
fbam_dir=/home/liud3/beegfs/cfDNA/data/cfdna_pipeline/fbam_dir

#Get file path for all bam and bai files
bam_file=$(find $bamdir -maxdepth 1 -name '*.bam')
bai_file=$(find $bamdir -maxdepth 1 -name '*.bai')
echo $bam_file 

#Define id variable as basename .bam for input for Rscript 
id=$(basename -a $bam_file) #add -a for multiple inputs  

#Run R script 
Rscript 01-filter_bam.R --id $id --bamdir $bamdir --fbam_dir $fbam_dir
echo Finished filtering original bam files!
