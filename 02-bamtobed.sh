#!/bin/bash

#load sam/bedtools 
module load bedtools
module load samtools 

#Define working directory
CWD=/home/liud3/beegfs/cfDNA/data/cfdna_pipeline
cd $CWD

#Directory of filtered bam files 
fbam_dir=/home/liud3/beegfs/cfDNA/data/cfdna_pipeline/fbam_dir

#Directory for output bedpe files for cfDNA fragments
beddir=/home/liud3/beegfs/cfDNA/data/cfdna_pipeline/02-bamtobed

#Loop to convert bam to bedpe
for f in $(find $fbam_dir -maxdepth 1 -iname "*.bam" -type f)
do	
  #get basename without bam suffix 
  id=$(basename -a -s f1.bam $f | cut -d "_" -f 1)
  echo $id
  #convert bam to bedpe in beddir directory 
  #bedpe files require to sort bam files 
  samtools sort -n $f | \
  #convert to bedpe
  bedtools bamtobed -i - -bedpe |\
  #calculate fragment size 
  awk 'OFS = "\t" {print $1="chr"$1, $2, $3, $4="chr"$4, $5, $6, $7, $8, $9, $10, $11=$6-$2}' - > $beddir/${id}.bedpe 
done 

echo Done creating bedpe files from filtered bam files  
