#!/bin/bash

#load sam/bedtools 
module load bedtools

#Define my working directory
CWD=/home/liud3/beegfs/cfDNA/data/cfdna_pipeline
cd $CWD

#Directory for input bed files of full fragment
bedpedir=/home/liud3/beegfs/cfDNA/data/cfdna_pipeline/02-bamtobed

#DIrectory for frag files 
beddir=/home/liud3/beegfs/cfDNA/data/cfdna_pipeline/03.5-frags_5mb

#Direcotry for output files 
outdir=/home/liud3/beegfs/cfDNA/data/cfdna_pipeline/04-filter_bedpe

#test=/home/liud3/beegfs/cfDNA/data/cfdna_pipeline/test

#for loop for all bed files
for f in $(find $bedpedir -maxdepth 1 -iname "*.bedpe" -type f)
do	
  #get basename without bam suffix 
  id=$(basename -a -s .bedpe $f)
  echo "id:" $id
  echo "f:" $f 
  #extract fragments from bedpe file that match with filtered frag file from 03.5 by matching RNAME field in both files
  #appends 10th column containing gc content 
  awk 'BEGIN{FS=OFS="\t"} FNR==NR{arr[$8]=$10;next} ($7 in arr){print $0, arr[$7]}' $beddir/${id}_frags_5mb.bed $f > $outdir/${id}_filtered.bedpe
done 
echo Done getting filtered bedpe files 
