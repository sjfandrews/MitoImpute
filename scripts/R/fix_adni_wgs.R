#!/usr/bin/Rscript
## Fix ADNI mitochondrial file to remove the last column, Mitochondria Information
require(tidyverse)
args = commandArgs(trailingOnly=TRUE)

filename <- args[1]
outfile <- args[2]

vcf <- read_tsv(filename, comment = '##')
vcf2 <- select(vcf, -Mitochondria_Information)
write_tsv(vcf2, append = T, outfile, col_names = T)



