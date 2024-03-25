#!/bin/bash

# Define project directory
PROJECT_DIR=/home/liud3/beegfs/cfDNA/protocol
cd "$PROJECT_DIR"

# Define input/output directories 
motifdir="${PROJECT_DIR}/20-motif_gc"
outdir="${PROJECT_DIR}/21-combine_motif"
mkdir -p "$outdir"  #Create directories if needed

module load R

# Initialize R script
Rscript 21-combine_motif.R \
--motifdir $motifdir \
--outdir $outdir
echo Done combining end motif data!
