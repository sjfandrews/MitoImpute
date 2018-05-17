#!/usr/bin/Rscript
##============================================================================##
##      R script to construct a recombiation map file for Impute2
##      Takes a VCF as input and outputs a .txt MAP file
##============================================================================##

suppressPackageStartupMessages(library(tidyverse))

args = commandArgs(trailingOnly = TRUE) # Set arguments from the command line
vcf.file = args[1] # SPECIFY IN THE VCF FILE
out.file = args[2] # SPECIFY THE OUT MAP FILE

message("")
if (is.na(out.file) == TRUE | is.null(out.file) == TRUE) {
  message('OUTPUT FILE NOT SPECIFIED')
} else {
  message(paste('OUTPUT FILE ASSIGNED TO ', out.file))
}

MT_RefPanel.raw <- read_tsv(vcf.file, comment = '##')
message(vcf.file, " SUCCESSFULLY PARSED")
message("")

MT_RefPanel <- MT_RefPanel.raw %>%
  separate(INFO, c('AC', 'AN', 'NS', 'AF', 'MAF',  'AC_Het', 'AC_Hom', 'AC_Hemi','HWE'), sep = ';') %>%
  select(`#CHROM`, POS, ID, REF, ALT, QUAL, FILTER, NS, AN, AF, MAF, AC) %>%
  rename(CHROM = `#CHROM`) %>%
  mutate(AF = as.numeric(gsub('AF=', '', AF))) %>%
	mutate(MAF = as.numeric(gsub('MAF=', '', MAF))) %>%
	mutate(NS = as.numeric(gsub('NS=', '', NS))) %>%
	mutate(AC = as.numeric(gsub('AC=', '', AC))) %>%
	mutate(AN = as.numeric(gsub('AN=', '', AN)))

print(MT_RefPanel)
cat('\n', length(unique(MT_RefPanel$POS)), 'Unique SNPs', '\n')

##  recombination map
# Fine-scale recombination map for the region to be analyzed. This file should have three columns:
# 1 = physical position (in base pairs)
# 2 = recombination rate between current position and next position in map (in cM/Mb)
# 3 = genetic map position (in cM).
# The file should also have a header line with an unbroken character string for each column (e.g., "position COMBINED_rate(cM/Mb) Genetic_Map(cM)").
# For the mitochondria, just set it to 0 and tbe bp for each snp

map.file <- tibble(position = MT_RefPanel$POS, COMBINED_rate.cM.Mb. = 0, Genetic_Map.cM. = MT_RefPanel$POS)
write_delim(map.file, out.file, delim = " ")

message(out.file, " SUCCESSFULLY WRITTEN")
