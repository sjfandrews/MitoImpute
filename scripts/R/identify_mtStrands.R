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
#b37 <- b37[-grep('.zip', b37)]

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

platforms <- do.call(bind_rows, b37.ls)
platforms %>% arrange(-n) %>% print(n = Inf)

write_tsv(platforms, paste0(out.dir, "MT_platforms.txt"))

strand.files = list.files(strand.dir, pattern = "_MT_snps.txt")
tmp = read.table(paste0(strand.dir, strand.files[2]), header = T)
m = matrix(0, nrow = length(strand.files), ncol = 16569)

for (i in 1:length(strand.files)) {
  tmp = read.table(paste0(strand.dir, strand.files[i]), header = T)
  td = c()
  for (j in sort(tmp$pos)) {
    m[i,j] = 1
    td = c(td, paste0("MT\t", j))
  }
  file.out = file(paste0(out.dir, "chip_regions/", strand.files[i]))
  writeLines(td, file.out)
  close(file.out)
}

row.names(m) = strand.files
mtn = c(1:16569)
mtn = unlist(lapply(mtn, function(x) paste0("MT.", x)))
colnames(m) = mtn

M = m[, colSums(m != 0) > 0]
write.csv(M, paste0(out.dir, "SNP_inclusion_matrix.csv"), quote = F)

for (i in 1:nrow(M)) {
  print(row.names(M)[i])
}

png(paste0(out.dir, "heat_map.png"), width = 30, height = 20, units = "in", res = 300)
heatmap(M, Rowv=NA, Colv=NA, scale="column")
dev.off()
