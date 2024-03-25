#!/bin/bash

# Define project directory
PROJECT_DIR=/home/liud3/beegfs/cfDNA/protocol
cd "$PROJECT_DIR"

# Define input/output directories 
fragdir="${PROJECT_DIR}/04-bins5mb"
outdir="${PROJECT_DIR}/10-frag_gc"
plotdir="${PROJECT_DIR}/10-gcbias_plots" #GC bias plots
statdir="${PROJECT_DIR}/10-frag_gc_stats" #Filtering statistics 

mkdir -p "$outdir" "$plotdir" "$statdir" #Create directories if needed

module load R

# Run R script 
Rscript 10-frags_gc.R \
--motifdir $fragdir \
--outdir $outdir \
--plotdir $plotdir \
--statdir $statdir
echo Done performing GC correction and summrazing fragments 
