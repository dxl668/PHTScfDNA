#!/bin/bash

# Define project directory
PROJECT_DIR=/home/liud3/beegfs/cfDNA
cd "$PROJECT_DIR"

# Define/make input/output directories 
bam_dir="${PROJECT_DIR}/data/files" #Aligned BAM files
fbam_dir="${PROJECT_DIR}/protocol/01-filter_bam" #Filtered BAM files
stat_dir="${PROJECT_DIR}/protocol/stat_reports" #Alignment statistics
mkdir -p "$bam_dir" "$fbam_dir" "$stat_dir" #Create directories if they do not exist

# Obtain file path for all BAM files 
bam_files=$(find $bam_dir -maxdepth 1 -name '*.bam')

# Load samtools module
module load samtools

# Iterate through each BAM file in directory
for bam_file in $bam_files; do
  fname=$(basename $bam_file)
  out_file="${fbam_dir}/${fname%.bam}_f1.bam"
  echo "Output file will be written to: $out_file"

  # Use samtools view to apply filtering criteria 
  samtools view -bh -f 2 -F 3844 -q 30 $bam_file $(seq 1 22)> $out_file
  
  # Summary statistics for the filtered BAM file
  stat_file="${stat_dir}/${fname%.bam}_flagstat.txt"
  echo "Generating stat report: $stat_file"
  samtools flagstat $out_file > $stat_file

done
echo "Finished filtering bam files"
