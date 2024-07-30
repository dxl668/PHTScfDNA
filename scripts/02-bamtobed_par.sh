#!/bin/bash
#SBATCH --job-name=02-bamtobed_par
#SBATCH --mail-type=ALL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH -p defq
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem 128000 # Memory request (128 GB)
#SBATCH -t 2-2:00 # Maximum execution time (D-HH:MM)
#SBATCH -o 02-bamtobed_par.out
#SBATCH -e 02-bamtobed_par.err

# Uncomment line below and add email to receive notifications
# #SBATCH --mail-user=your.email@example.com

# Define project directory
PROJECT_DIR="/home/liud3/beegfs/cfDNA/protocol"
cd "$PROJECT_DIR"

# Define directories for input and output files
fbam_dir="${PROJECT_DIR}/01-filter_bam" # Filtered BAM files
bed_dir="${PROJECT_DIR}/02-bamtobed_par" # BEDPE files
mkdir -p "$bed_dir" # Create directories if they do not exist
 
module load bedtools
module load samtools 
module load parallel

bam_to_bedpe() {
    f="$1"
    id=$(basename "$f" .bam | cut -d "_" -f 1)
    echo "id: $id"
    # Files must be sorted bam before converting to BEDPE 
    samtools sort -n "$f" | \
    bedtools bamtobed -i - -bedpe | \
    # Calculate insert size and add 'chr' to chromosome fields
    awk 'OFS="\t" {print "chr"$1, $2, $3, "chr"$4, $5, $6, $7, $8, $9, $10, $6-$2}' > "$bed_dir/${id}.bedpe"
}

export -f bam_to_bedpe
export bed_dir

# Using GNU Parallel to process the files in parallel
find "$fbam_dir" -maxdepth 1 -iname "*.bam" -type f | parallel bam_to_bedpe {}

echo "Done converting BAM files to BEDPE files"
