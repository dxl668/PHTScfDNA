#!/bin/bash
#SBATCH --job-name=21-combine_motif
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
#SBATCH -o 21-combine_motif.out
#SBATCH -e 21-combine_motif.err

# Uncomment line below and add email to receive notifications
# #SBATCH --mail-user=your.email@example.com

# Define project directory
PROJECT_DIR=/user/project
cd "$PROJECT_DIR"

# Define input/output directories 
motifdir="${PROJECT_DIR}/20-motif_gc"
outdir="${PROJECT_DIR}/21-combine_motif"
mkdir -p "$outdir"  #Create directories if needed

module load R/4.2.3

# Define the path to the R script
R_SCRIPT="${PROJECT_DIR}/scripts/21-combine_motif.R"

# Initialize R script
Rscript "$R_SCRIPT" \
--motifdir $motifdir \
--outdir $outdir
echo Done combining end motif data!
