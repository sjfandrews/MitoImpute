##======================================================================##
##        R Script for 1) counting number of MT snps on each platform   ##
##        and 2) Writing out a list of platforms with MT SNPs           ##
##======================================================================##

# Load Tidyverse 
suppressPackageStartupMessages(library(tidyverse))

## Read in Arguments 
# 1: Current working directory 
args = commandArgs(trailingOnly=TRUE)
rwd <- args[1]

# Read in file paths of individule platform strand files
b37 <- list.dirs(path = paste0(rwd, "/data/platforms"), full.names = TRUE, recursive = F)

## for each strand file
#   read in file and rename columns
#   if MT is present in CHROM column
#     count number of MT SNPs
#   else create empty tibble
#   write out
b37.ls <- lapply(b37, function(x){
  cat('IN PROGRESS: ', x)
  df <- read_tsv(paste0(x, '/platform.strand'), col_names = F, col_types = 'ccnncc') %>% 
    rename(SNP = X1, chr = X2, pos = X3, match = X4, strand = X5)
  if('MT' %in% df$chr){
    platform <- str_split_fixed(x, "platforms/", n =2)[,2]
    out <- count(df, chr) %>% filter(chr == 'MT') %>% mutate(platform = platform)
    cat('\t DONE \n'); out
  } else {
    platform <- str_split_fixed(x, "platforms/", n =2)[,2]
    out <- tibble(chr = 'MT', n = 0, platform = platform)
    cat('\t DONE \n'); out
  }
})

# Gather counts into tibble
platforms <- do.call(bind_rows, b37.ls)

# Write out tibble of number of MT snps on each platform
write_tsv(platforms %>% arrange(-n), paste0(rwd, '/data/platforms/Nsnps_Mt_platforms.txt'))
cat('Written: \t /data/platforms/Nsnps_Mt_platforms.txt \n')

# Write out list of platforms containing > 1 MT SNPs
write_tsv(platforms %>% filter(n >1) %>% select(platform), paste0(rwd, '/data/platforms/Mt_platforms.txt'), col_names = F)
cat('Written: \t /data/platforms/Mt_platforms.txt \n')

