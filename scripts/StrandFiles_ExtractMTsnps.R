##########################################################################
##        R Script for extracting a list of MT SNPs on each platform    ##
##########################################################################

## Load Tidyverse
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))

## Arguments
# 1: Stand file
# 2: Output file
args <- commandArgs(trailingOnly = T)
strand.file <- args[1]
out.file <- args[2]

## Extract MT snps from strand file
#     Read in strand file
#     Rename columns
#     Filter for MT in CHROM
#     Select on CHROM and POS columns
#     Arrange tibble by POS
#     write out file
message(sprintf("IN PROGRESS: %s", strand.file))
read_tsv(strand.file, col_names = F, col_types = "ccnncc") %>%
  rename(SNP = X1, CHROM = X2, POS = X3, match = X4, strand = X5) %>%
  filter(CHROM == "MT") %>%
  select(CHROM, POS) %>%
  arrange(POS) %>%
  write_tsv(out.file, col_names = F)
