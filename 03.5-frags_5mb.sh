#!/bin/bash

#load sam/bedtools 
module load bedtools

#Define my working directory
CWD=/home/liud3/beegfs/cfDNA/data/cfdna_pipeline
cd $CWD

#Direcotry for fragments (stard, end, gc content information)
beddir=/home/liud3/beegfs/cfDNA/data/cfdna_pipeline/03-frags_gc

#Direcotry for frag files 
outdir=/home/liud3/beegfs/cfDNA/data/cfdna_pipeline/03.5-frags_5mb

#5mb bins 
bins5mb=/home/liud3/beegfs/cfDNA/data/cfdna_pipeline/files/bins5mb_filtered.bed

for f in $(find $beddir -maxdepth 1 -iname "*.bed" -type f)
do	
  #get basename without bam suffix 
  id=$(basename -a -s _frags_gc.bed $f)
  echo "id:" $id
  echo "f:" $f
  
  #intersect fragments with 5mb bins and obtain unique intersections
  bedtools intersect -a $bins5mb -b $f -wo |\
  #due to uniq syntax, base pair overlaps $11 was causing issues0, moved to column $5
  awk 'OFS="\t" {ov=$11; print $1, $2, $3, $4, ov, $5, $6, $7, $8, $9, $10}' - |\
  sort -k9,9 | uniq -u -f8 |\
  awk 'OFS="\t" {print $1, $2, $3, $4, $6, $7, $8, $9, $10, $11}' - > tmp_${id}_unique.bed 
  echo "temp_uniq:" tmp_${id}_unique.bed 

  #intersect fragments with 5mb and obtain duplicated fragments due to those lying on boudnary of two bins
  bedtools intersect -a $bins5mb -b $f -wo |\
  #due to uniq syntax, base pair overlaps $11 was causing issues0, moved to column $5
  awk 'OFS="\t" {ov=$11; print $1, $2, $3, $4, ov, $5, $6, $7, $8, $9, $10}' - |\
  sort -k9,9 | uniq -D -f8 |\
  sort -k9,9 -k4,4nr |\
  
  #compare base pair overlap of duplicate read; filter out fragment with smallest overlap  
  awk 'BEGIN{FS=OFS="\t"; prev=""; maxval=0} {if ($9==prev) {if ($4>maxval) {maxval=$4; line=$0}} else {if (prev!="") print line; maxval=$4; line=$0; prev=$9}}
     END {if (prev!="") print line}' - |\
  #remove $5 (base overlap)
  awk 'OFS="\t" {print $1, $2, $3, $4, $6, $7, $8, $9, $10, $11}' - > tmp_${id}_dupl.bed 
  echo "temp_uniq:" tmp_${id}_dupl.bed 

  #combine tmp files
  cat tmp_${id}_unique.bed tmp_${id}_dupl.bed > $outdir/${id}_frags_5mb.bed

  #delete tmp files 
  rm tmp_${id}_unique.bed tmp_${id}_dupl.bed
  echo done with $f 

done 
echo Done intersecting fragments with 5mb and resolving duplicates
