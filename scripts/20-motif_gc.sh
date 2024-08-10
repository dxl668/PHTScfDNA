#!/bin/bash
#SBATCH --job-name=20-motif_gc
#SBATCH --mail-type=ALL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH -p defq
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem 4000 # Memory request (4 GB)
#SBATCH -t 0-08:00 # Maximum execution time (D-HH:MM) 
#SBATCH -o 20-motif_gc.out
#SBATCH -e 20-motif_gc.err

# Uncomment line below and add email to receive notifications
# #SBATCH --mail-user=your.email@example.com

# Define project directory
PROJECT_DIR=/user/project
cd "$PROJECT_DIR"

# Define input/output directories 
motifdir="${PROJECT_DIR}/08-motif_merge"
outdir="${PROJECT_DIR}/20-motif_gc"
plotdir="${PROJECT_DIR}/20-gcbias_plots" #GC bias plots
statdir="${PROJECT_DIR}/20-motif_gc_stats" #Filtering statistics 

mkdir -p "$outdir" "$plotdir" "$statdir" #Create directories if needed

module load R/4.2.3

# Define the path to the R script
R_SCRIPT="${PROJECT_DIR}/scripts/20-motif_gc.R"

#Run R script 
Rscript "$R_SCRIPT" \
--motifdir $motifdir \
--outdir $outdir \
--plotdir $plotdir \
--statdir $statdir
echo Done performing GC correction and summarizing end motifs
