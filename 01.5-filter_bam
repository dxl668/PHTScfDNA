#!/bin/bash

##Use samtools to filter bam files for fragments <= 1000 bp

#load samtools 
module load samtools

#Directory for input bam files
bam_dir=/home/liud3/beegfs/cfDNA/data/cfdna_pipeline/fbam_dir

#Directory for output bam files (filterd)
outdir=/home/liud3/beegfs/cfDNA/data/cfdna_pipeline/01.5-filter_bam

for f in $(find $bam_dir -maxdepth 1 -iname "*.bam" -type f)
do	
  #get basename without bam suffix 
  id=$(basename -a -s f1.bam $f | cut -d "_" -f 1)
  echo $id
  samtools view -h $f | \
  awk 'substr($0,1,1)=="@" || ($9<=1000) || ($9>=-1000)' | \
  samtools view -b - > $outdir/${id}_1000bp.bam
  echo $f completed
done
echo Done filtering bam files for fragments less than 1000 bp 
