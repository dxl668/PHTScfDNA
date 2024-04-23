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

# Get file names for end motif data
files <- list.files(motifdir, full.names=TRUE)
files.list <- lapply(files, readRDS) # Read RDS

# Get sample ID
id <- basename(files)
id <- gsub("_4bpmotif_gc.rds", "", id)

# Combined tibble
tib.list <- lapply(files.list, as_tibble)
names(tib.list) <- id

tib.list <- map2(tib.list, .y = names(tib.list), ~ mutate(.x, id = .y)) %>%
  bind_rows() %>% select(id, everything())

# Output file 
out.file <- file.path(outdir, paste0("endmotif_4bp_gc_summary.rds"))
  
# Save as RDS
saveRDS(tib.list, out.file)
