#!/bin/bash
#SBATCH --job-name=10-frag_gc
#SBATCH --mail-type=ALL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH -p bigmem
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem 6000 # Memory request (6 GB)
#SBATCH -t 0-12:00 # Maximum execution time (D-HH:MM) 
#SBATCH -o 10-frags_gc.out
#SBATCH -e 10-frags_gc.err

# Uncomment line below and add email to receive notifications
# #SBATCH --mail-user=your.email@example.com

# Define project directory
PROJECT_DIR=/user/project
cd "$PROJECT_DIR"

# Define input/output directories 
fragdir="${PROJECT_DIR}/04-bins5mb"
outdir="${PROJECT_DIR}/10-frags_gc"
plotdir="${PROJECT_DIR}/10-gcbias_plots" #GC bias plots
statdir="${PROJECT_DIR}/10-frags_gc_stats" #Filtering statistics 

mkdir -p "$outdir" "$plotdir" "$statdir" #Create directories if needed

module load R/4.2.3

# Define the path to the R script
R_SCRIPT="${PROJECT_DIR}/scripts/10-frags_gc.R"

# Run R script 
Rscript "$R_SCRIPT" \
--motifdir $fragdir \
--outdir $outdir \
--plotdir $plotdir \
--statdir $statdir
echo Done performing GC correction and summrazing fragments 
