#!/bin/bash

#load sam/bedtools 
module load bedtools

#Define my working directory 
CWD=/home/liud3/beegfs/cfDNA/data/cfdna_pipeline
cd $CWD

#Directory for input bed files of full fragment
beddir=/home/liud3/beegfs/cfDNA/data/cfdna_pipeline/06-motif_fasta

#Direcotry for output files 
outdir=/home/liud3/beegfs/cfDNA/data/cfdna_pipeline/07-motif_merge

#for loop for merging 2bp files 
#changed this to .bed rather than .txt bc i named files wrong by accident 
for f in $(find $beddir -maxdepth 1 -iname "*2bp_r1_fa.bed" -type f)
do  
  #echo $f
  id=$(basename -a -s .bed $f | cut -d "_" -f 1)
  echo $id
  #merge pos/neg strand into one bed file
  cat $beddir/${id}_2bp_r1_fa.bed $beddir/${id}_2bp_r2_fa.bed > $outdir/${id}_2bp_motif.bed
done 

echo Done merging fa files 
