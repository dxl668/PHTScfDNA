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
fragdir <- options.args[1]
outdir <- options.args[2]

## For testing locally
#setwd("W:/Shared/LRI/Labs/engclab/LAB_MEMBER_CENTRAL/Darren_Liu/3_cf-DNA/protocol_ms")
#fragdir <- "data/10-frags_gc"

# Get file names for fragment/bin data
files <- list.files(fragdir, full.names=TRUE)
files.list <- lapply(files, readRDS)

# Get sample ID 
id <- basename(files)
id <- gsub("_bins5mb.rds", "", id)

# Convert list into tibbles 
tib.list <- lapply(files.list, as_tibble)
names(tib.list) <- id

# Combine tibble
tib.list <- map2(tib.list, .y = names(tib.list), ~ mutate(.x, id = .y)) %>%
  bind_rows() %>% select(id, everything())

# Output file path 
out.file <- file.path(outdir, paste0("binfrag_summary.rds"))

#save as RDS
saveRDS(tib.list, out.file)