#!/usr/bin/Rscript
##============================================================================##
##      R script to construct a recombiation map file for Impute2
##      and a strand file for Impute2
##      Takes a VCF as input and outputs a .txt MAP and Strand File
##============================================================================##

suppressPackageStartupMessages(library(tidyverse))

vcf.file = snakemake@input[['vcf']] # SPECIFY IN THE VCF FILE
outfile.map = snakemake@output[['map']] # SPECIFY THE OUT MAP FILE
outfile.strand = snakemake@output[['strand']] # SPECIFY THE OUT STRAND FILE
no_genos = T # Don't bother reading genotypes

message("")
if (is.na(outfile.map) == TRUE | is.null(outfile.map) == TRUE) {
  message("OUTPUT FILE NOT SPECIFIED")
} else {
  message(paste("OUTPUT FILE ASSIGNED TO ", outfile.map))
}

message("")
if (is.na(outfile.strand) == TRUE | is.null(outfile.strand) == TRUE) {
  message("OUTPUT FILE NOT SPECIFIED")
} else {
  message(paste("OUTPUT FILE ASSIGNED TO ", outfile.strand))
}

numstr <- . %>%
  na.omit %>% # remove missing
  {suppressWarnings(as.numeric(.))} %>% # convert to numeric
  is.na %>% # see which failed
  any %>% # see if any failed
  `!` # reverse to success

sep_infolist <- function(df) {
 info <- df %>%
   select(CHROM, POS, REF, ALT, INFO) %>% #don't fill up memory
   mutate(INFO = strsplit(INFO, ";")) %>% #info to list of k=v pairs
   unnest(INFO) %>% # each info field gets its own row
   separate(INFO, into = c("key", "val"), sep = "=") %>% # separate k=v pairs
   spread(key, val) %>% # spread back to one row per record
   mutate_if(numstr, as.numeric) #correct column types
  df %>%
    left_join(info, by = c("CHROM", "POS", "REF", "ALT")) %>%
    select(names(info), everything())
}

if (no_genos) {
  vcfcols <- cols_only(
    POS = "i", REF = "c", "ALT" = "c", FILTER = "c",
    FORMAT = "c", INFO = "c", "#CHROM" = "c"
  )
} else {
  vcfcols <- cols(
    .default = "i", POS = "i", REF = "c", "ALT" = "c",
    FILTER = "c", FORMAT = "c", INFO = "c", "#CHROM" = "c"
  )
}

vcf.raw <- read_tsv(vcf.file, comment = "##",
                    col_types = vcfcols, na = ".", trim_ws = T)
message(vcf.file, " SUCCESSFULLY PARSED")
message("")

message("PARSING INFO")
mt_panel <- vcf.raw %>%
  rename(CHROM = "#CHROM") %>%
  sep_infolist

print(mt_panel)
cat("\n", length(unique(mt_panel$POS)), "Unique SNPs", "\n")

##  recombination map
# Fine-scale recombination map for the region to be analyzed. This file should have three columns:
# 1 = physical position (in base pairs)
# 2 = recombination rate between current position and next position in map (in cM/Mb)
# 3 = genetic map position (in cM).
# The file should also have a header line with an unbroken character string for each column (e.g., "position COMBINED_rate(cM/Mb) Genetic_Map(cM)").
# For the mitochondria, just set it to 0 and tbe bp for each snp

tibble(position = mt_panel$POS,
  COMBINED_rate.cM.Mb. = 0,
  Genetic_Map.cM. = mt_panel$POS) %>%
  write_delim(outfile.map, delim = " ")

message(outfile.map, " SUCCESSFULLY WRITTEN")

## Strand File
# File showing the strand orientation of the SNP allele codings in the -g file, relative to a fixed reference point.
# Each SNP occupies one line, and the file should have two columns:
# (i) the base pair position of the SNP and
# (ii) the strand orientation ('+' or '-') of the alleles in the genotype file
# the columns should be separated by a single space.
# For the mitochondria, everything is + as it is alligned to the rCRS

tibble(position = mt_panel$POS, strand = "+") %>%
  write_delim(outfile.strand, delim = " ", col_names = F)

message(outfile.strand, " SUCCESSFULLY WRITTEN")
