#!/bin/bash
#SBATCH --job-name=LiuD_cfDNA_11-combine_bin
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
#SBATCH -o 11-combine_bin.out
#SBATCH -e 11-combine_bin.err

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