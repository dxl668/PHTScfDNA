#!/bin/bash
#SBATCH --job-name=LiuD_cfDNA_10-frag_gc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liud3@ccf.org
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH -p bigmem
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem 128000 # Memory request (128 GB)
#SBATCH -t 2-2:00 #Maximum execution time (D-HH:MM) - 2 days-2hours
#SBATCH -o 10-frags_gc.out
#SBATCH -e 10-frags_gc.err

# Define project directory
PROJECT_DIR=/home/liud3/beegfs/cfDNA/protocol
cd "$PROJECT_DIR"

# Define input/output directories 
fragdir="${PROJECT_DIR}/04-bins5mb"
outdir="${PROJECT_DIR}/10-frags_gc"
plotdir="${PROJECT_DIR}/10-gcbias_plots" #GC bias plots
statdir="${PROJECT_DIR}/10-frags_gc_stats" #Filtering statistics 

mkdir -p "$outdir" "$plotdir" "$statdir" #Create directories if needed

module load R

# Run R script 
Rscript 10-frags_gc.R \
--motifdir $fragdir \
--outdir $outdir \
--plotdir $plotdir \
--statdir $statdir
echo Done performing GC correction and summrazing fragments 