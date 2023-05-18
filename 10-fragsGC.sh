#!/bin/bash

#load R
module load R

#Define my working directory 
CWD=/home/liud3/beegfs/cfDNA/data/cfdna_pipeline
cd $CWD

#Define directory for frags
motifdir=/home/liud3/beegfs/cfDNA/data/cfdna_pipeline/03.5-frags_5mb

#Directory for RDS files 
outdir=/home/liud3/beegfs/cfDNA/data/cfdna_pipeline/10-fragsGC		

#Directory for GC bias plots 
imagedir=/home/liud3/beegfs/cfDNA/data/cfdna_pipeline/10-fragsGC_plots

#Run R script 
Rscript 10-fragsGC.R --motifdir $motifdir --outdir $outdir --imagedir $imagedir
echo Done reading frags, performing GC correction, and calculating ratios 
