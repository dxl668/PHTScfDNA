library(tidyverse)

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
motifdir <- options.args[1]
outdir <- options.args[2]

#get motif data 
files <- list.files(motifdir, full.names=TRUE)

#read in files as list
files.list <- lapply(files, readRDS)

#get basename of files for file ids
id <- basename(files)
id <- gsub("_bin_5mb.rds", "", id)

#convert list into tibles and names with sample ID
tib.list <- lapply(files.list, as_tibble)
names(tib.list) <- id

#combined tibble
tib.list <- map2(tib.list, .y = names(tib.list), ~ mutate(.x, id = .y)) %>%
  bind_rows() %>% select(id, everything())

#create rds file for output 
out.file <- file.path(outdir, paste0("binfrag_summary.rds"))
  
#save as RDS
saveRDS(tib.list, out.file)
