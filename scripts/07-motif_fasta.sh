#!/bin/bash
#SBATCH --job-name=07-motif_fasta
#SBATCH --mail-type=ALL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH -p defq
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem 6000 # Memory request (6 GB)
#SBATCH -t 0-12:00 # Maximum execution time (D-HH:MM) 
#SBATCH -o 07-motif_fasta.out
#SBATCH -e 07-motif_fasta.err

# Uncomment line below and add email to receive notifications
# #SBATCH --mail-user=your.email@example.com

# Define project directory
PROJECT_DIR=/user/project
cd "$PROJECT_DIR"

# Define input/output directories 
bed_dir="${PROJECT_DIR}/06-endmotif_bed"
out_dir="${PROJECT_DIR}/07-motif_fasta" #Intermediate file for troubleshooting 
out2_dir="${PROJECT_DIR}/08-motif_merge"
hg19="${PROJECT_DIR}/files/hg19.fa"
mkdir -p "$out_dir" "$out2_dir" #Create directories if needed

module load bedtools/2.29.0

# Iterate through each end motif BED file 
for f in $(find $bed_dir -maxdepth 1 -iname "*.bed" -type f)
do	
  id=$(basename -a -s .bed $f)
  echo $id
  # Extract nucleotide sequence for end motifs 
  bedtools getfasta -fi $hg19 -bed $f -s -bedOut -fo |\
  # Filter bed file to only have chr, start, end, strand, gc, motif sequence
  awk 'OFS="\t" {print $1, $2, $3, $6, $7, $8, toupper($9)}' - > $out_dir/${id}_fa.bed
done 
echo Done extracting end motif

# Define n-mer length (same as 06-endmotif_bed.sh)
nmer=4

# Iterate through 07-motif_fasta and merge 
for f in $(find $out_dir -maxdepth 1 -iname "*${nmer}bp_r1_fa.bed" -type f)
do  
  id=$(basename "$f" "_${nmer}bp_r1_fa.bed")
  echo $id
  # Merge read pair files together 
  cat $out_dir/${id}_${nmer}bp_r1_fa.bed $out_dir/${id}_${nmer}bp_r2_fa.bed > $out2_dir/${id}_${nmer}bp_motif.bed
done 
echo Done merging fa files 
