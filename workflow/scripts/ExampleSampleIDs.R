## Extract Thousand Genome Samples for use in example Sample panle

library(tidyverse)

infile = snakemake@input[['info']]
outfile = snakemake@output[['out']]

eur <- readxl::read_xlsx(infile) %>%
  select(Sample, FID = `Family ID`, Population) %>% 
  filter(Population %in% c('CEU', 'GBR', 'TSI', 'FIN', 'IBS'))

eur %>% select(Sample) %>% write_tsv(outfile, col_names = F)
