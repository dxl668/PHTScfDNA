#!/bin/bash
#SBATCH --job-name=05-filter_bedpe
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
#SBATCH -o 05-filter_bedpe.out
#SBATCH -e 05-filter_bedpe.err

# Uncomment line below and add email to receive notifications
# #SBATCH --mail-user=your.email@example.com

# Define project directory
PROJECT_DIR=/user/project
cd "$PROJECT_DIR"

# Define input/output directories 
bedpe_dir="${PROJECT_DIR}/02-bamtobed"
frags_dir="${PROJECT_DIR}/04-bins5mb"
out_dir="${PROJECT_DIR}/05-filter_bedpe"
mkdir -p "$out_dir" #Create directory if they do not exist

# Load bedtools 
module load bedtools/2.29.0

# Iterate through BEDPE files 
for f in $(find $bedpe_dir -maxdepth 1 -iname "*.bedpe" -type f)
do	
  # Obtain sample name
  id=$(basename -a -s .bedpe $f)
  # Extract fragments from BEDPE file from filtered/binned fragments
  # Appends GC content to 10th column  
  awk '
    BEGIN {
      # Set the input and output field separators to tab
      FS=OFS="\t"
    }
    FNR==NR {
      # For the first, create an array with QNAME as the key and GC content as the value
      arr[$8]=$10
      # Skip to the next record without executing the rest of the code
      next
    }
    ($7 in arr) {
      # For the second file (BEDPE), if the QNAME exists in the array,
      # print the current line and append the GC content from the array
      print $0, arr[$7]
    }
  ' $frags_dir/${id}_frags_5mb.bed $f > $out_dir/${id}_filtered.bedpe
done 
echo Extracting filtered/intersected fragments complete 

