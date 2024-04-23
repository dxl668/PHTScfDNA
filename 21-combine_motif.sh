#!/bin/bash
#SBATCH --job-name=LiuD_cfDNA_21-combine_motif
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
#SBATCH -o 21-combine_motif.out
#SBATCH -e 21-combine_motif.err

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