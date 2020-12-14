## Extract Samples for use in example reference panle

library(tidyverse)

infile = snakemake@input[['ref']]
outfile = snakemake@output[['out']]

set.seed(333)
mtref <- read_tsv(infile) %>% 
  rename(haplogroup = Haplogroup) %>% 
  mutate(macro = ifelse(str_detect(haplogroup, "^L|^HV"),
                        substr(haplogroup, start = 1, stop = 2),
                        substr(haplogroup, start = 1, stop = 1))) %>% 
  filter(macro %in% c('I', "W", "X", "J", "T", "U", "K", "HV", "H", "V")) %>% 
  sample_frac(., 0.1, replace = FALSE)

mtref %>% select(SampleID) %>% write_tsv(outfile)