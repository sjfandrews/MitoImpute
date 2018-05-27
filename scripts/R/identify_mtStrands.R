library(tidyverse)

out.dir = "/Volumes/TimMcInerney/Other/Shea_Imputation/metadata/chip_strand"
if (!endsWith(out.dir, "/")) {
  out.dir = paste0(out.dir, "/")
}
strand.dir = "/Volumes/TimMcInerney/Other/Shea_Imputation/metadata/chip_strand/unzipped"
if (!endsWith(strand.dir, "/")) {
  strand.dir = paste0(strand.dir, "/")
}

b37 <- list.files(strand.dir, pattern = "*.strand")
b37 <- b37[-grep('.zip', b37)]

b37.ls <- lapply(b37, function(x){
  cat(x, 'IN PROGRESS \n')
  df <- read_tsv(paste0(strand.dir, x), col_names = F, col_types = 'ccnncc') %>% 
    rename(SNP = X1, chr = X2, pos = X3, match = X4, strand = X5)
  if('MT' %in% df$chr){
    df %>% filter(chr == 'MT') %>% write_tsv(paste0(strand.dir, x, "_MT_snps.txt"))
    out <- count(df, chr) %>% filter(chr == 'MT') %>% mutate(platform = x)
    out
  } else {
    tibble(chr = 'MT', n = 0, platform = x)
  }
})

platforms <- do.call(bind_rows, test)
platforms %>% arrange(-n) %>% print(n = Inf)

write_tsv(platforms, paste0(out.dir, "MT_platforms.txt"))