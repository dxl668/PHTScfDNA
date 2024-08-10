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
#SBATCH --mem 4000 # Memory request (4 GB)
#SBATCH -t 0-08:00 # Maximum execution time (D-HH:MM)
#SBATCH -o 11-combine_bin.out
#SBATCH -e 11-combine_bin.err

# Uncomment line below and add email to receive notifications
# #SBATCH --mail-user=your.email@example.com

# Define project directory
PROJECT_DIR=/user/project
cd "$PROJECT_DIR"

# Define input/output directories 
fragdir="${PROJECT_DIR}/10-frags_gc"
outdir="${PROJECT_DIR}/11-combine_bin"
mkdir -p "$outdir"  #Create directories if needed

module load R/4.2.3

# Define the path to the R script
R_SCRIPT="${PROJECT_DIR}/scripts/11-combine_bin.R"

# Initialize R script
Rscript "$R_SCRIPT"  \
--fragdir $fragdir \
--outdir $outdir
echo Done combining end motif data!
