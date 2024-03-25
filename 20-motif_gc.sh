#!/bin/bash

# Define project directory
PROJECT_DIR=/home/liud3/beegfs/cfDNA/protocol
cd "$PROJECT_DIR"

# Define input/output directories 
motifdir="${PROJECT_DIR}/08-motif_merge"
outdir="${PROJECT_DIR}/20-motif_gc"
plotdir="${PROJECT_DIR}/20-gcbias_plots" #GC bias plots
statdir="${PROJECT_DIR}/20-motif_gc_stats" #Filtering statistics 

mkdir -p "$outdir" "$plotdir" "$statdir" #Create directories if needed

module load R

#Run R script 
Rscript 20-motif_gc.R \
--motifdir $motifdir \
--outdir $outdir \
--plotdir $plotdir \
--statdir $statdir
echo Done performing GC correction and summarizing end motifs
