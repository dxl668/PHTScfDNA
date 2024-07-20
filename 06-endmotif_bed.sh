#!/bin/bash
#SBATCH --job-name=06-endmotif_bed
#SBATCH --mail-type=ALL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH -p defq
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem 128000 # Memory request (128 GB)
#SBATCH -t 2-2:00 # Maximum execution time (D-HH:MM) 
#SBATCH -o 06-endmotif_bed.out
#SBATCH -e 06-endmotif_bed.err

# Uncomment line below and add email to receive notifications
# #SBATCH --mail-user=your.email@example.com

# Define project directory
PROJECT_DIR=/home/liud3/beegfs/cfDNA/protocol
cd "$PROJECT_DIR"

# Define input/output directories 
bedpe_dir="${PROJECT_DIR}/05-filter_bedpe"
out_dir="${PROJECT_DIR}/06-endmotif_bed"
mkdir -p "$out_dir" #Create directory does not exist

# Define n-mer length
nmer=4

# Iterate through each filtered BEDPE file 
for f in $(find $bedpe_dir -maxdepth 1 -iname "*.bedpe" -type f)
do	
  id=$(basename -a -s _filtered.bedpe $f)
  echo $id
  # Create two separate bed files for each 5' n-mer end motif 
  awk -v nmer="$nmer" 'OFS="\t" {print $1, $2, $2+nmer, $7, $8, $9, $11, $12}' $f > $out_dir/${id}_${nmer}bp_r1.bed
  awk -v nmer="$nmer" 'OFS="\t" {print $4, $6-nmer, $6, $7, $8, $10, $11, $12}' $f > $out_dir/${id}_${nmer}bp_r2.bed
done 
echo Done creating fragment end motif BED files

