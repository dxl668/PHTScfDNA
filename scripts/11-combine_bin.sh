#!/bin/bash
#SBATCH --job-name=11-combine_bin
#SBATCH --mail-type=ALL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH -p bigmem
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem 128000 # Memory request (128 GB)
#SBATCH -t 2-2:00 #Maximum execution time (D-HH:MM)
#SBATCH -o 11-combine_bin.out
#SBATCH -e 11-combine_bin.err

# Uncomment line below and add email to receive notifications
# #SBATCH --mail-user=your.email@example.com

# Define project directory
PROJECT_DIR=/home/liud3/beegfs/cfDNA/protocol
cd "$PROJECT_DIR"

# Define input/output directories 
fragdir="${PROJECT_DIR}/10-frags_gc"
outdir="${PROJECT_DIR}/11-combine_bin"
mkdir -p "$outdir"  #Create directories if needed

module load R

# Initialize R script
Rscript 11-combine_bin.R \
--fragdir $fragdir \
--outdir $outdir
echo Done combining end motif data!
