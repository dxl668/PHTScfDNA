#!/bin/bash

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

