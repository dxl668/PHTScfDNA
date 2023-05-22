#!/bin/bash

#load R
module load R

#Define my working directory 
CWD=/home/liud3/beegfs/cfDNA/data/cfdna_pipeline
cd $CWD


#Define directory for 
motifdir=/home/liud3/beegfs/cfDNA/data/cfdna_pipeline/10-fragsGC

#Define 
outdir=/home/liud3/beegfs/cfDNA/data/cfdna_pipeline/11-combine_bin

#Run R script 
Rscript 11-combine_bin.R --motifdir $fragdir --outdir $outdir

echo Done combine dataframes with gc-corrected framgnets into single dataframe!
