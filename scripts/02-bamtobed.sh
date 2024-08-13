#!/bin/bash
#SBATCH --job-name=02-bamtobed
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
#SBATCH -o 02-bamtobed.out
#SBATCH -e 02-bamtobed.err

# Uncomment line below and add email to receive notifications
# #SBATCH --mail-user=your.email@example.com

# Define project directory
PROJECT_DIR=/user/project
cd "$PROJECT_DIR"

# Define directories for input and output files
fbam_dir="${PROJECT_DIR}/01-filter_bam" #Filtered BAM files
bed_dir="${PROJECT_DIR}/02-bamtobed" #BEDPE files
mkdir -p "$bed_dir" #Create directories if they do not exist

# Load sam/bedtools 
module load bedtools/2.29.0
module load samtools/1.16.1 

# Convert BAM files BEDPE files using bedtools
for f in $(find $fbam_dir -maxdepth 1 -iname "*.bam" -type f)
do	
  id=$(basename -a -s f1.bam $f | cut -d "_" -f 1)
  echo "id:" $id 
  # Files must be sorted bam before converting to BEDPE 
  samtools sort -n $f | \
  # Convert to BEDPE file
  bedtools bamtobed -i - -bedpe |\
  # Calculate insert size and add 'chr' to chromosome field
  awk 'OFS = "\t" {print $1="chr"$1, $2, $3, "chr"$4, $5, $6, $7, $8, $9, $10, $6-$2}' - > $bed_dir/${id}.bedpe 
done 
echo Done converting BAM files to BEDPE files
