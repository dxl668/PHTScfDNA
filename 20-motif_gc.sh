#!/bin/bash
#SBATCH --job-name=LiuD_cfDNA_20-motif_gc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liud3@ccf.org
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH -p defq
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem 128000 # Memory request (128 GB)
#SBATCH -t 2-2:00 # Maximum execution time (D-HH:MM)
#SBATCH -o 20-motif_gc.out
#SBATCH -e 20-motif_gc.err

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