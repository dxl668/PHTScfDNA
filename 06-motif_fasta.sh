#!/bin/bash

#load sam/bedtools 
module load bedtools

#Define my working directory 
CWD=/home/liud3/beegfs/cfDNA/data/cfdna_pipeline
cd $CWD

#Directory for input bed files of full fragment
beddir=/home/liud3/beegfs/cfDNA/data/cfdna_pipeline/05-end_motif_bed

#Direcotry for output files 
outdir=/home/liud3/beegfs/cfDNA/data/cfdna_pipeline/06-motif_fasta

#Fasta file for hg19 genome 
hg19=/home/liud3/beegfs/cfDNA/data/end_motif/files/hg19.fa

#for loop for all bed files
for f in $(find $beddir -maxdepth 1 -iname "*.bed" -type f)
do	
  #get basename without bed suffix 
  id=$(basename -a -s .bed $f)
  echo "id:" $id
  echo "f:" $f
  #get nucleotide sequence for 2bp end motifs 
  bedtools getfasta -fi $hg19 -bed $f -s -bedOut -fo |\
  #filter bed file to only have chr, start, end, strand, gc, end mofit
  awk 'OFS="\t" {nuc=$9; print $1, $2, $3, $4=$6, $5=$7, $6=$8, $7=toupper(nuc)}' - > $outdir/${id}_fa.bed
  echo done with $f 
done 

echo Done with getting gc content and motif sequences
