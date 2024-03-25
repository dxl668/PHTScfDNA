#!/bin/bash

# Define project directory
PROJECT_DIR=/home/liud3/beegfs/cfDNA/protocol
cd "$PROJECT_DIR"

# Define input/output directories 
bed_dir="${PROJECT_DIR}/06-endmotif_bed"
out_dir="${PROJECT_DIR}/07-motif_fasta" #Intermediate file for troubleshooting 
out2_dir="${PROJECT_DIR}/08-motif_merge"
hg19="${PROJECT_DIR}/files/hg19.fa"
mkdir -p "$out_dir" "$out2_dir" #Create directories if needed

module load bedtools

# Iterate through each end motif BED file 
for f in $(find $bed_dir -maxdepth 1 -iname "*.bed" -type f)
do	
  id=$(basename -a -s .bed $f)
  echo "id:" $id
  # Extract nucleotide sequence for end motifs 
  bedtools getfasta -fi $hg19 -bed $f -s -bedOut -fo |\
  # Filter bed file to only have chr, start, end, strand, gc, motif sequence
  awk 'OFS="\t" {print $1, $2, $3, $6, $7, $8, toupper($9)}' - > $out_dir/${id}_fa.bed
  echo Processing $f complete 
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
