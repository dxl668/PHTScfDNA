#!/bin/bash
#SBATCH --job-name=03-filter_frags_gc
#SBATCH --mail-type=ALL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH -p defq
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem 6000 # Memory request (6 GB)
#SBATCH -t 0-12:00 # Maximum execution time (D-HH:MM) 
#SBATCH -o 03-filter_frags_gc.out
#SBATCH -e 03-filter_frags_gc.err

# Uncomment line below and add email to receive notifications
# #SBATCH --mail-user=your.email@example.com

# Define project directory
PROJECT_DIR=/user/project
cd "$PROJECT_DIR"

# Define directories for input and output files
bed_dir="${PROJECT_DIR}/02-bamtobed" #BEDPE files
out_dir="${PROJECT_DIR}/03-filter_frags_gc" #filtered fragments
hg19="${PROJECT_DIR}/files/hg19.fa" #h19 reference genome
filter="${PROJECT_DIR}/files/hg19-blacklist.v2.bed" #blacklisted region
gap="${PROJECT_DIR}/files/gaps.hg19.bed" #genomic gaps
mkdir -p "$out_dir" #Create directories if needed

module load bedtools

# Iterate through all BED files
for f in $(find $bed_dir -maxdepth 1 -iname "*.bedpe" -type f)
do	
  # Get sample ID
  id=$(basename -a -s .bedpe $f)
  echo "id:" $id
  # BED file format: fragment start/end interval/QNAME/frag size 
  awk 'OFS = "\t" {print $1, $2, $6, $7, $11}' $f |\
  # Sort by chromosome and start position
  sort -k1,1 -k2,2n |\
  # Filter blacklisted/gap regions
  bedtools subtract -a - -b $filter -A |\
  bedtools subtract -a - -b $gap -A |\
  # Obtain fragment GC content 
  bedtools nuc -fi $hg19 -bed - | \
  # Add GC content($7)
  awk 'OFS = "\t" {print $1, $2, $3, $4, $5, $7}' - > $out_dir/${id}_frags_gc.bed
  echo $f Complete
done 
echo Done filtering and calculating fragment GC content


