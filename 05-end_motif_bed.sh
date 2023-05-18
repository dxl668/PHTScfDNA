#!/bin/bash

#load sam/bedtools 
module load bedtools

#Define my working directory 
CWD=/home/liud3/beegfs/cfDNA/data/cfdna_pipeline
cd $CWD

#Directory for input bedpe files with fragment info
beddir=/home/liud3/beegfs/cfDNA/data/cfdna_pipeline/04-filter_bedpe

#Direcotry for output files 
outdir=/home/liud3/beegfs/cfDNA/data/cfdna_pipeline/05-end_motif_bed

#For loop separating out bedpe file into separate bed files with read1/read2 with gc content of full fragment
for f in $(find $beddir -maxdepth 1 -iname "*.bedpe" -type f)
do	
  id=$(basename -a -s _filtered.bedpe $f)
  echo $id
  #bedfile for read 1/read2
  awk 'OFS="\t" {print $1, $2, $3=$2+2, $7, $8, $9, $11, $12, $13}' $f > $outdir/${id}_2bp_r1.bed 
  awk 'OFS="\t" {print $4, $5=$6-2, $6, $7, $8, $10, $11, $12, $13}' $f > $outdir/${id}_2bp_r2.bed
done 
echo Done creating bed files for 2bp end motifs for read1 and read2
