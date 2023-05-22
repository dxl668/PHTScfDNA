#load R
module load R

#Define my working directory 
CWD=/home/liud3/beegfs/cfDNA/data/cfdna_pipeline
cd $CWD


#Define directory for sample end motifs 
motifdir=/home/liud3/beegfs/cfDNA/data/cfdna_pipeline/20-motif_gc_correct

#Define out directory
outdir=/home/liud3/beegfs/cfDNA/data/cfdna_pipeline/21-combine_motif

#Run R script 
Rscript 21-combine_motif.R --motifdir $motifdir --outdir $outdir
echo Done combine dataframes for motif summaries!
