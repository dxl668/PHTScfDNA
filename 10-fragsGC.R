library(getopt)
library(GenomicRanges)
library(rtracklayer)
library(Homo.sapiens)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Rsamtools)
class(Homo.sapiens)
library(devtools)
library(biovizBase)
library(tidyverse)
library(ggplot2)


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
imagedir <- options.args[3] 


#turn off scientific notation
options(scipen = 999)

#turn on sci notiation
#options(scipen = 0)


#Function for GC correction bias
#Bias is the list of gc content per bin

seqlast <- function (from, to, by) 
{
  vec <- do.call(what = seq, args = list(from, to, by))
  if ( tail(vec, 1) != to ) {
    return(c(vec, to))
  } else {
    return(vec)
  }
}

gc.correct <- function(coverage, bias) {
  
  #get min/max of gc - vector in increment of 0.001 
  i <- seqlast(min(bias, na.rm=TRUE), max(bias, na.rm=TRUE), by = 0.001)
  
  #obtain counts correlcted to GC content
  coverage.trend <- loess(coverage ~ bias)
  
  #create a model of counts correlated with GC content
  coverage.model <- loess(predict(coverage.trend, i) ~ i)
  
  #get predicted counts correlated with GC 
  coverage.pred <- predict(coverage.model, bias)
  
  #subtract predicted coverage explained by GC from model
  coverage.corrected <- coverage - coverage.pred + median(coverage)
}

#Create filepath to RDS files found in fragdir 
fragfiles <- list.files(fragdir, full.names=TRUE)
print(fragfiles)

#Get name of the sample
id <- basename(fragfiles)
id <- gsub("_frags_5mb.bed", "", id)

for(x in 1:length(fragfiles)) {
  frags <- read.table(fragfiles[x],header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
  #remove unwanted columns 

  #Rename columns and round gc content 
  frags <- frags %>% 
    select(-V5,-V6, -V7, -V8) %>%
    rename( #rename column headers
      "chr" = "V1",
      "start" = "V2",
      "end" = "V3",
      "arm" = "V4", 
      #"bin_gc" = "V5", 
      #"bin" = "V7",
      "frag" = "V9", 
      "gc" = "V10" #frag gc
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
    
  #create df with counts for each gc strata per nucleosome fraction 
  df_gc <- frags %>% 
    filter(frag >= 100 & frag <= 650) %>%  
    mutate(
      size = case_when(
        frag >= 100 & frag <=150 ~ "short1", 
        frag >= 151 & frag <=250 ~ "long1",
        frag >= 251 & frag <=300 ~ "short2",
        frag >= 301 & frag <=450 ~ "long2",
        frag >= 451  & frag <=500 ~ "short3",
        frag >= 501 & frag <=650 ~ "long3",
      ),
      frac = case_when(
        frag >= 100 & frag <=250 ~ "mono", 
        frag >= 251 & frag <=450 ~ "di", 
        frag >= 451 & frag <=650 ~ "tri")
        #frag >= 651 & frag <=1000 ~ "high_mw")
    ) %>% 
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
    ) %>% ungroup() 
  
  #select variables to merge   
  df_gc1 <- df_gc %>% select(frac, gc, gc.scale)
  
  #plot gc bias 
  gc.bias.overall<- df_gc %>%
    group_by(frac) %>% 
    group_map(~ggplot(.) + 
                aes(x = gc, y = count) +
                geom_line() + 
                ggtitle(.y[[1]])+ 
                geom_line(aes(x = gc, y = count.corrected, color = "red"))+
                labs(color = "corrected.counts"))
  #get end motifs for labeling 
  frac <- unique(df_gc$frac)
  
  #save GC plot
  for(i in 1:length(gc.bias.overall)){
    #create png file for output images 
    print(frac[i]) #THIS MIGHT NOT RUN 
    out.gc_plot <- file.path(imagedir, paste0(id[x], "/", id[[x]], frac[i], "_2bpmotif_gc_overallcorrect.png"))
    print(out.gc_plot)
    #save p lot 
    ggsave(filename = out.gc_plot, plot = gc.bias.overall[[i]], width = 7, height = 5)
  }
  
  #get frag count per bin and gc strata
  frags1 <- frags %>% 
    filter(frag >=100 & frag <= 650) %>% 
    group_by(bin, frag, gc) %>%
    summarise(
      counts = n(), 
      chr = first(chr), 
      start = first(start),
      end = first(end),
      arm = first(arm)
    ) %>%
    mutate(
      size = case_when(
        frag >= 100 & frag <=150 ~ "short1", 
        frag >= 151 & frag <=250 ~ "long1",
        frag >= 251 & frag <=300 ~ "short2",
        frag >= 301 & frag <=450 ~ "long2",
        frag >= 451  & frag <=500 ~ "short3",
        frag >= 501 & frag <=650 ~ "long3",
      ),
      frac = case_when(
        frag >= 100 & frag <=250 ~ "mono", 
        frag >= 251 & frag <=450 ~ "di", 
        frag >= 451 & frag <=650 ~ "tri")) %>%
    ungroup()
        #frag >= 651 & frag <=1000 ~ "high_mw"))
  
  #Correcting  fragment per bin
  frags2 <- frags1 %>% 
    inner_join(df_gc1, by=c("frac", "gc")) %>%
    mutate(
      counts.corrected = counts*gc.scale
    ) %>% ungroup() %>%
    group_by(bin, frag) %>% 
    summarise(
      counts = sum(counts), 
      counts.corrected = sum(counts.corrected), 
      chr = first(chr), 
      start = first(start),
      end = first(end),
      arm = first(arm)
    ) %>% ungroup()
  
  #For QC, get frag distribution for counts and counts.corrected
  frag.dist <- frags2 %>% 
    ungroup() %>%
    group_by(frag) %>%
    summarise(
      counts = sum(counts),
      counts.corrected = sum(counts.corrected)
    )%>%
    pivot_longer(cols = c(counts, counts.corrected),
                 names_to = "count_type",
                 values_to = "count_value")
  
  #plot frag.distribution
  p.distr <- ggplot(frag.dist, aes(x = frag, y = count_value, color = count_type)) + 
    geom_line() 
  
  #save plot 
  out.gc_plot <- file.path(imagedir, paste0(id[x], "/", id[x], "_frag_dist.png"))
  print(out.gc_plot) 
  ggsave(filename = out.gc_plot, plot = p.distr, width = 7, height = 5)
  
  #Reshap df using pivot wider
  frags3 <- frags2 %>% pivot_wider(names_from = frag, values_from = c(counts, counts.corrected),values_fill = 0)
    
  #create dataframe with short/long/ratio columns
  #merge with frags 3
  frags4 <- frags2 %>% 
    mutate(
      size = case_when(
        frag >= 100 & frag <=150 ~ "short1", 
        frag >= 151 & frag <=250 ~ "long1",
        frag >= 251 & frag <=300 ~ "short2",
        frag >= 301 & frag <=450 ~ "long2",
        frag >= 451  & frag <=500 ~ "short3",
        frag >= 501 & frag <=650 ~ "long3")
      ) %>%
    group_by(bin, size) %>%
    summarise(
      counts = sum(counts),
      counts.corrected = sum(counts.corrected),
      chr = first(chr), 
      start = first(start),
      end = first(end),
      arm = first(arm), 
      size = first(size), 
      #frac = first(frac)
    ) %>% ungroup()
  
  #Reshape with pivot wider
  frags5 <- frags4 %>% pivot_wider(names_from = size, 
                                 values_from = c(counts, counts.corrected),
                                 names_sep = "_",
                                 values_fill = 0)
  
  # Calculate ratio/ratio.corrected
  frags5 <- frags5 %>%
    mutate(
      ratio1 = counts_short1/counts_long1,
      ratio1.corrected = counts.corrected_short1/counts.corrected_long1,
      ratio2 = counts_short2/counts_long2, 
      ratio2.corrected = counts.corrected_short2/counts.corrected_long2,
      ratio3 = counts_short3/counts_long3,
      ratio3.corrected = counts.corrected_short3/counts.corrected_long3
    )
  
  # Merge 
  frags6 <- frags5 %>% 
    inner_join(frags3, by=c("bin", "chr", "start", "end", "arm")) 
  
  #save as RDS
  saveRDS(frags6, file.path(outdir, paste0(id[x], "_bin_5mb.rds")))
}
