library(tidyverse)
library(ggplot2)


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
plotdir <- options.args[3] 
statdir <- options.args[4]

### For testing on local drive
setwd("W:/Shared/LRI/Labs/engclab/LAB_MEMBER_CENTRAL/Darren_Liu/3_cf-DNA/protocol_ms")
motifdir <- "data/08-motif_merge"

### Functions
# Add last element of vector in i for GC-correction model 
seqlast <- function (from, to, by) {
  vec <- do.call(what = seq, args = list(from, to, by))
  if ( tail(vec, 1) != to ) {
    return(c(vec, to))
  } else {
    return(vec)
  }
}

# Function for GC bias correction 
gc.correct <- function(coverage, bias) {
  # Get min/max of GC - vector in increment of 0.01 
  i <- seqlast(min(bias, na.rm=TRUE), max(bias, na.rm=TRUE), by = 0.01)
  # Obtain counts correlated to GC content
  coverage.trend <- loess(coverage ~ bias)
  # Create a model of counts correlated with GC content
  coverage.model <- loess(predict(coverage.trend, i) ~ i)
  # Get predicted counts correlated with GC 
  coverage.pred <- predict(coverage.model, bias)
  # Subtract predicted coverage explained by GC from model
  coverage.corrected <- coverage - coverage.pred + median(coverage)
}

# Get end motif files
files <- list.files(motifdir, full.names=TRUE)

# Get sample ID
id <- basename(files)
id <- gsub("_4bp_motif.bed", "", id)

# File path for output 
out.file <- file.path(outdir, paste0(id, "_4bpmotif_gc.rds"))
stat.file <- file.path(statdir, paste0(id, "_4bpmotif_stat.rds"))

# Iterate through each BED file
for (x in 1:length(files)){
  df_motif <- read.table(files[1],header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
  df_motif <- df_motif %>% 
    dplyr::rename(
      "chr" = "V1",
      "start" = "V2",
      "end" = "V3",
      "strand" = "V4", 
      "frag_size" = "V5", 
      "gc" = "V6",
      "end_motif" = "V7"
    )
  
  # Calculate total # end motifs prior to filtering
  n_motifs <- nrow(df_motif)
  df_stats <- tibble(n_motifs) # Create df_stats to store statistics
  
  # Calculate # and % of sequences with "N" and "NA"
  n_N <- sum(grepl("N", df_motif$end_motif, fixed = TRUE))
  n_NA <- sum(is.na(df_motif$end_motif))
  N_percent <- (n_N / n_motifs) * 100
  NA_percent <- (n_NA / n_motifs) * 100
  N_NA_percent_filter <- N_percent + NA_percent
  df_stats <- df_stats %>% add_column(n_N, n_NA, N_NA_percent_filter)
  
  # Print counts and percentages for sequences with "N" and "NA"
  cat("Number of sequences with 'N':", n_N,"(",N_percent,"%)\n")
  cat("Number of sequences with NA:", n_NA, "(",NA_percent, "%)\n")
  
  # Filter out sequences with "N" and "NA"
  df_motif <- df_motif %>% ungroup() %>% 
    filter(!grepl("N", end_motif), !is.na(end_motif)) 
  
  # Calculate total # end motifs prior to filtering by frag size
  n_motifs <- nrow(df_motif) 
  df_stats <- df_stats %>% add_column(n_motif_filter = n_motifs)
  
  # Filter motifs by frag size 
  df_motif <- df_motif %>% ungroup() %>% 
    filter(frag_size >=100 & frag_size <= 650) %>%
    # Define nucleosome fractions 
    mutate(
      frac = case_when(
        frag_size >= 100 & frag_size <=250 ~ "mono", 
        frag_size >= 251 & frag_size <=450 ~ "di", 
        frag_size >= 451 & frag_size <=650 ~ "tri")
    ) 
  
  # Statistics for fragment size filtering
  n_filter_size <- n_motifs - nrow(df_motif)
  n_filter_size_percent <- (n_filter_size / n_motifs) * 100
  df_stats <- df_stats %>% add_column(n_filter_size, n_filter_size_percent, n_motif_size = nrow(df_motif))
  
  cat("# of end motifs before N/NA filtering:", n_motifs, "\n")
  cat("# of end motifs filtered out by frag size:", n_filter_size, "\n")
  cat("% end motifs filtered out:", n_filter_size_percent, "%\n")
  
  # Obtain end motif counts grouped by cucleosome Fraction and GC content 
  df_motif <- df_motif %>% 
    select(end_motif, frag_size, gc, frac) %>%
    mutate(gc = round(gc, 2)) %>% 
    group_by(end_motif, frac, gc) %>%
    summarise(
      count = n()
    ) %>%
    ungroup() 
  
  # Obtain GC-correction scalar calc. for each GC strata per each nuc. fraction  
  df_gc <- df_motif %>%
    group_by(frac, gc) %>% 
    summarise(
      frac.count = sum(count) 
    ) %>%
    ungroup() %>%
    group_by(frac) %>%
    mutate(
      frac.total = sum(frac.count), # Total count per nuc. fraction
      cum_count = cumsum(frac.count), # Cumulative count
      percentile = cum_count/frac.total*100 # Cumulative %
    ) %>%
    filter(percentile >= 5.00 & percentile <= 95.00) %>% 
    mutate(
      frac.count.corrected = gc.correct(frac.count, gc),
      total.corrected = sum(frac.count.corrected),
      gc.scale = frac.count.corrected / frac.count # Calculate weight for each GC strata
    ) %>%
    select(frac, gc, frac.count, frac.count.corrected, gc.scale)
  
  
  # Calculate corrected end motif counts 
  df_motif <- inner_join(df_motif, df_gc, by=c("frac", "gc")) %>% 
    mutate(
      count.corrected = gc.scale*count #GC correction 
    ) %>% ungroup() %>%
    group_by(frac) %>%
    mutate( #calculate total counts per nuc fraction
      frac.total = sum(count), 
      frac.total.corrected = sum(count.corrected)
    ) %>% ungroup() %>%
    mutate( #calculate total couts overall 
      total.count = sum(count), 
      total.corrected = sum(count.corrected)
    ) %>% ungroup() 
  
  
  print("GC-bias plot for each nucleosome fraction")  
  
  # GC-bias plot 
  gc.bias.overall <- df_gc %>%
    group_by(frac) %>% 
    group_map(~ggplot(.) + 
                aes(x = gc, y = frac.count) +
                geom_line() + 
                ggtitle(.y[[1]])+ 
                geom_line(aes(x = gc, y = frac.count.corrected, color = "red"))+
                labs(color = "corrected.counts"))
  
  frac <- unique(df_motif$frac) # For labeling 
  
  # Save GC-bias plot for each nucleosome fraction 
  for(i in 1:length(gc.bias.overall)){
    print(frac[i])
    out.gc_plot <- file.path(plotdir, paste0(id[x], "/", id[[x]], frac[i], "_4bpmotif_gc.png"))
    print(out.gc_plot)
    ggsave(filename = out.gc_plot, plot = gc.bias.overall[[i]], width = 7, height = 5)
  }
  
  # save files
  saveRDS(df_motif, out.file[x])
  saveRDS(df_stats, stat.file[x])
}


