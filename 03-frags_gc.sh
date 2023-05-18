#!/bin/bash

#load sam/bedtools 
module load bedtools

#Define my working directory 
CWD=/home/liud3/beegfs/cfDNA/data/cfdna_pipeline
cd $CWD

#Directory for input bedpe files 
beddir=/home/liud3/beegfs/cfDNA/data/cfdna_pipeline/02-bamtobed

#Direcotry for output files 
outdir=/home/liud3/beegfs/cfDNA/data/cfdna_pipeline/03-frags_gc

#Fasta file for hg19 genome 
hg19=/home/liud3/beegfs/cfDNA/data/end_motif/files/hg19.fa

#File for filteredd blacklistedd region
filter=/home/liud3/beegfs/cfDNA/data/cfdna_pipeline/files/filter_combined.bed

#File for 5mb bins with gc/mappability 
#bins5mb=/home/liud3/beegfs/cfDNA/data/cfdna_pipeline/files/bins5mb_gcmap.bed

for f in $(find $beddir -maxdepth 1 -iname "*.bedpe" -type f)
do	
  #get basename without bam suffix 
  id=$(basename -a -s .bedpe $f)
  echo "id:" $id
  echo "f:" $f
  #create bed4 file format using fragment start/end interval/name/frag size 
  awk 'OFS = "\t" {print $1, $2, $3=$6, $7, $11}' $f |\
  #filter out black listed regions (-A, removes fragments with ANY overlap)
  bedtools subtract -a - -b $filter -A |\
  #get gc content of fragments 
  bedtools nuc -fi $hg19 -bed - | \
  #format output get gc content ($7)
  awk 'OFS = "\t" {print $1, $2, $3, $4, $5, $7}' - > $outdir/${id}_frags_gc.bed
  #intersect with 5mb bins (already filtered by gc/map); -wa and -wb option to get both 
  echo done with $f 
done 
echo Done getting fragment gc content 
