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
fragdir <- options.args[1]
outdir <- options.args[2]
plotdir <- options.args[3] 
statdir <- options.args[4]

# Turn off scientific notation
options(scipen = 999)

# Turn on sci notiation
# Options(scipen = 0)

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

# Get filtered fragment files 
files <- list.files(fragdir, full.names=TRUE)

# Get sample ID
id <- basename(files)
id <- gsub("_frags_5mb.bed", "", id)

# File path for output 
out.file <- file.path(outdir, paste0(id, "_bins5mb.rds"))
stat.file <- file.path(statdir, paste0(id, "_bins5mb_stat.rds"))

for(x in 1:length(files)) {
  df.frags <- read.table(files[x],header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
  
  df.frags <- df.frags %>% 
    select(-V5,-V6, -V7, -V8) %>%
    rename( 
      "chr" = "V1",
      "start" = "V2",
      "end" = "V3",
      "arm" = "V4", 
      "frag_size" = "V9", 
      "gc" = "V10" 
    ) %>% 
    mutate(
      gc = round(gc, 2)
    ) %>%
    group_by(chr, start) %>%
    arrange(chr, start) %>% 
    mutate(
      bin = cur_group_id()
    ) %>% 
    relocate(bin, .before = "chr")
  
  # Calculate total # frags prior to filtering
  n_frags <- nrow(df.frags)
  df_stats <- tibble(n_frags) # Create df_stats to store statistics
  cat("Number of fragments:", n_frags)
  
  # Filter motifs by frag size 
  df.frags <- df.frags %>% ungroup() %>% 
    filter(frag_size >=100 & frag_size <= 650) 
  
  # Statistics for fragment size filtering
  n_filter_size <- n_frags - nrow(df.frags)
  n_filter_size_percent <- (n_filter_size / n_frags) * 100
  df_stats <- df_stats %>% add_column(n_filter_size, n_filter_size_percent, n_frags2 = nrow(df.frags))
  cat("# of fragments filtered out:", n_filter_size, "(",n_filter_size_percent,"%) \n")
  
  # Obtain fragment counts per bin 
  df.frags <- df.frags %>% 
    group_by(bin, frag_size, gc, chr, start, end, arm) %>%
    summarise(
      counts = n()
    ) %>%
    mutate(
      frac = case_when( # Define nucleosome fraction
        frag_size >= 100 & frag_size <=250 ~ "mono", 
        frag_size >= 251 & frag_size <=450 ~ "di", 
        frag_size >= 451 & frag_size <=650 ~ "tri")) %>%
    ungroup()
  
  # Obtain GC-correction scalar calc. for each GC strata per each nuc. fraction 
  df_gc <- df.frags %>% 
    group_by(frac, gc) %>%
    summarise(
      count = n()
    ) %>% ungroup()%>%
    group_by(frac) %>%
    mutate(
      frac.total = sum(count), 
      cum_count = cumsum(count),
      percentile = cum_count/frac.total*100
    ) %>%
    filter(percentile >= 5.00 & percentile <= 95.00)%>%
    mutate(
      count.corrected = gc.correct(count, gc),
      total.corrected = sum(count.corrected),
      gc.scale = count.corrected / count
    ) %>%
    select(frac, gc, count, count.corrected, gc.scale) %>% 
    ungroup() 
  
  # GC-bias plot 
  gc.bias.overall<- df_gc %>%
    group_by(frac) %>% 
    group_map(~ggplot(.) + 
                aes(x = gc, y = count) +
                geom_line() + 
                ggtitle(.y[[1]])+ 
                geom_line(aes(x = gc, y = count.corrected, color = "red"))+
                labs(color = "corrected.counts"))
  
  frac <- unique(df_gc$frac) # For labeling
  
  # Save GC-bias plot for each nucleosome fraction
  for(i in 1:length(gc.bias.overall)){
    #create png file for output images 
    print(frac[i]) 
    out.gc_plot <- file.path(plotdir, paste0(id[x], "/", id[[x]], frac[i], "_5mb_frags_gc.png"))
    print(out.gc_plot)
    ggsave(filename = out.gc_plot, plot = gc.bias.overall[[i]], width = 7, height = 5, create.dir = TRUE)
  }
  
  
  # Calculate corrected fragment counts 
  df.frags <- df.frags %>% 
    inner_join(df_gc, by=c("frac", "gc")) %>%
    mutate(
      counts.corrected = counts*gc.scale
    ) %>% ungroup() %>%
    group_by(bin, frag_size, chr, start, end, arm) %>% 
    summarise(
      counts = sum(counts), 
      counts.corrected = sum(counts.corrected), 
    ) %>% ungroup()
  
  # Reshape to wide format 
  df.frags1 <- df.frags %>% 
    pivot_wider(names_from = frag_size, values_from = c(counts, counts.corrected),values_fill = 0)
  
  df.frags2 <- df.frags %>% 
    mutate(
      size = case_when(
        frag_size >= 100 & frag_size <=150 ~ "short1", 
        frag_size >= 151 & frag_size <=250 ~ "long1",
        frag_size >= 251 & frag_size <=300 ~ "short2",
        frag_size >= 301 & frag_size <=450 ~ "long2",
        frag_size >= 451  & frag_size <=500 ~ "short3",
        frag_size >= 501 & frag_size <=650 ~ "long3")
    ) %>%
    group_by(bin, size, chr, start, end, arm) %>%
    summarise(
      counts = sum(counts),
      counts.corrected = sum(counts.corrected),
    ) %>% 
    ungroup() %>% 
    pivot_wider(names_from = size, 
                values_from = c(counts, counts.corrected),
                names_sep = "_",
                values_fill = 0)
  
  # Calculate ratio of short/long 
  df.frags2 <- df.frags2 %>%
    mutate(
      ratio1 = counts_short1/counts_long1,
      ratio1.corrected = counts.corrected_short1/counts.corrected_long1,
      ratio2 = counts_short2/counts_long2, 
      ratio2.corrected = counts.corrected_short2/counts.corrected_long2,
      ratio3 = counts_short3/counts_long3,
      ratio3.corrected = counts.corrected_short3/counts.corrected_long3
    )
  
  # Merge 
  df.frags2 <- df.frags2 %>% 
    inner_join(df.frags1, by=c("bin", "chr", "start", "end", "arm")) 
  
  # Save files
  saveRDS(df.frags2, out.file[x])
  saveRDS(df_stats, stat.file[x])
}




