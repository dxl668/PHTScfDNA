#!/bin/bash
#SBATCH --job-name=LiuD_cfDNA_04-bins5mb
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
#SBATCH -t 2-2:00 #Maximum execution time (D-HH:MM) - 2 days-2hours
#SBATCH -o 04-bins5mb.out
#SBATCH -e 04-bins5mb.err

# Define project directory 
PROJECT_DIR=/home/liud3/beegfs/cfDNA/protocol
cd "$PROJECT_DIR"

# Define directories for input and output files
in_dir="${PROJECT_DIR}/03-filter_frags_gc"
out_dir="${PROJECT_DIR}/04-bins5mb"
bins5mb="${PROJECT_DIR}/files/bins5mb_filtered.bed" 
mkdir -p "$out_dir" #Create directories if they do not exist

module load bedtools 

# Iterate through all bed files 
for f in $(find $in_dir -maxdepth 1 -iname "*.bed" -type f)
do	
  # Obtain ID from filename 
  id=$(basename -a -s _frags_gc.bed $f)
  echo $id
  echo $f

  # Intersect fragments with 5Mb bins and unique intersections
  bedtools intersect -a $bins5mb -b $f -wo |\
  # Moved base overlap ($11) to $5 to avoid issues with uniq syntax 
  awk 'OFS="\t" {ov=$11; print $1, $2, $3, $4, ov, $5, $6, $7, $8, $9, $10}' - |\
  # Sort by QNAME ($9) and obtain unique fragments 
  sort -k9,9 |\
  uniq -u -f8 |\
  # Remove $5 (base pair overlap)
  awk 'OFS="\t" {print $1, $2, $3, $4, $6, $7, $8, $9, $10, $11}' - > tmp_${id}_uniq.bed 
  echo "temp_uniq:" tmp_${id}_uniq.bed 

  # Intersect fragments with 5Mb bins and obtain duplicated fragments at bin boundaries
  bedtools intersect -a $bins5mb -b $f -wo |\
  # Moved base overlap ($11) to $5 to avoid issues with uniq syntax 
  awk 'OFS="\t" {ov=$11; print $1, $2, $3, $4, ov, $5, $6, $7, $8, $9, $10}' - |\
  sort -k9,9 |\
  # Obtain duplicate fragment based QNAME ($9)
  uniq -D -f8 |\
  # Sort by QNAME and base pair overlap ($5)
  sort -k9,9 -k5,5nr |\
  # Compare base pair overlap of duplicate fragment & filter out fragment with smallest overlap  
  # If tie, will use the first bin
  awk 'BEGIN {
    FS=OFS="\t"; prev=""; maxval=0} {
  	if ($9==prev) {
  		if ($5>maxval) {maxval=$5; line=$0}
  	} else {
  		if (prev!="") print line; 
  			maxval=$5; line=$0; prev=$9
  	}
  }
  END {if (prev!="") print line}' - |\
  # Remove $5 (base pair overlap)
  awk 'OFS="\t" {print $1, $2, $3, $4, $6, $7, $8, $9, $10, $11}' - > tmp_${id}_dupl.bed 
  echo "temp_dupl:" tmp_${id}_dupl.bed 

  # Combine unique and de-duplicated fragments
  cat tmp_${id}_uniq.bed tmp_${id}_dupl.bed > $out_dir/${id}_frags_5mb.bed

  # Delete temporary files 
  rm tmp_${id}_uniq.bed tmp_${id}_dupl.bed
  echo done with $f 

done 
echo Done intersecting fragments with 5Mb bins and resolving duplicates
