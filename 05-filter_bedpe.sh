#!/bin/bash

# Define project directory
PROJECT_DIR=/home/liud3/beegfs/cfDNA/protocol
cd "$PROJECT_DIR"

# Define input/output directories 
bedpe_dir="${PROJECT_DIR}/02-bamtobed"
frags_dir="${PROJECT_DIR}/04-bins5mb"
out_dir="${PROJECT_DIR}/05-filter_bedpe"
mkdir -p "$out_dir" #Create directory if they do not exist

# Load bedtools 
module load bedtools

# Iterate through BEDPE files 
for f in $(find $bedpe_dir -maxdepth 1 -iname "*.bedpe" -type f)
do	
  # Obtain sample name
  id=$(basename -a -s .bedpe $f)
  echo "id:" $id
  echo "f:" $f 
  # Extract fragments from BEDPE file from filtered fragments that were intersected with the 5Mb bins 
  # Appends GC content to 10th column  
  awk '
    BEGIN {
      # Set the input and output field separators to tab
      FS=OFS="\t"
    }
    FNR==NR {
      # For the first, create an array with QNAME as the key and GC content as the value
      arr[$8]=$10
      # Skip to the next record without executing the rest of the code
      next
    }
    ($7 in arr) {
      # For the second file (BEDPE), if the QNAME exists in the array,
      # print the current line and append the GC content from the array
      print $0, arr[$7]
    }
  ' $frags_dir/${id}_frags_5mb.bed $f > $out_dir/${id}_filtered.bedpe
done 
echo Extracting filtered/intersected fragments complete 

