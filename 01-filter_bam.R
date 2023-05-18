library(Rsamtools)
library(getopt)
library(GenomicAlignments)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(RCurl)

### Used for getting information from shell script
args <- commandArgs(trailingOnly = TRUE)
hh <- paste(unlist(args), collapse = " ")
listoptions <- unlist(strsplit(hh, "--"))[-1]
options.args <- sapply(listoptions, function(x) {
  unlist(strsplit(x, " "))[-1]
})
options.names <- sapply(listoptions, function(x) {
  option <- unlist(strsplit(x, " "))[1]
})
names(options.args) <- unlist(options.names)
id <- options.args[1]
bamdir <- options.args[2]
fbam_dir <- options.args[3] 

#Get sample ID
id <- unlist(id) #Need to unlist for paths to be correct
print(id)

### Create paths to bam/bai files
bamfile <- file.path(bamdir, id) 
indexed.bam <- gsub("bam", "bai", bamfile)

# Create bam file for output file destination
sample <- gsub(".bam", "", id)

out.file <- file.path(fbam_dir, paste0(sample, "f1.bam"))
print(out.file)

#creates GRanges object of chr1-22 with their seqlengh 
chromosomes <- GRanges(paste0(1:22),
                       IRanges(1, seqlengths(Hsapiens)[1:22]))

#Define what to be filtered by ScanBamParam
which <- chromosomes

# Get filtering parameters
param <- ScanBamParam(which = which,
                      flag = scanBamFlag(isPaired = TRUE, 
                                         isProperPair = TRUE, 
                                         isDuplicate = FALSE,
                                         isSecondaryAlignment = FALSE,
                                         isUnmappedQuery = FALSE),
                      mapqFilter = 30)

#For loop to filter all bam files 
for (x in 1:length(bamfile)) {
  filtered_bam <- filterBam(bamfile[x], out.file[x], indexed.bam[x],param = param)

}
