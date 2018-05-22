library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

rwd <- args[1]

b37 <- list.dirs(path = paste0(rwd, "/data/platforms"), full.names = TRUE, recursive = F)

b37.ls <- lapply(b37, function(x){
  cat(x, 'IN PROGRESS \n')
  df <- read_tsv(paste0(x, '/platform.strand'), col_names = F, col_types = 'ccnncc') %>% 
        rename(SNP = X1, chr = X2, pos = X3, match = X4, strand = X5)
  if('MT' %in% df$chr){
    platform <- str_split_fixed(x, "platforms/", n =2)[,2]
    df %>% 
      filter(chr == 'MT') %>% 
      rename(CHROM = chr, POS = pos) %>% 
      select(CHROM, POS) %>%
      write_tsv(paste0(rwd, platform, "_MT_snps.txt"))
    out <- count(df, chr) %>% filter(chr == 'MT') %>% mutate(platform = platform)
    out
  } else {
    platform <- str_split_fixed(x, "platforms/", n =2)[,2]
    out <- tibble(chr = 'MT', n = 0, platform = platform)
    out
  }
})

platforms <- do.call(bind_rows, test)

write_tsv(platforms %>% arrange(-n), paste0(rwd, 'data/platforms/Nsnps_Mt_platforms.txt'))
write_tsv(platforms %>% filter(n >1) %>% select(platform), paste0(rwd, 'data/platforms/Mt_platforms.txt'), col_names = F)


