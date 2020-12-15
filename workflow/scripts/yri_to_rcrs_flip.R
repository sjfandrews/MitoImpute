#!/usr/bin/env Rscript

## Snakemake INput
reference <- snakemake@input[['referenceSnps']]
bimfile <- snakemake@input[['bim']]
outfile <- snakemake@output[['bim']]
log_path = snakemake@log[[1]]

## Logging messages
con <- file(log_path, open = "wt")
sink(con, type = "message")

suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))

#-----------------------------------------------------##
## 	Converting Yoruba (YRI - AF347015) SNP positions
##	to the revised cambrige reference sequence (rCRS)
##	https://www.mitomap.org/MITOMAP/YorubanConversion
##  http://haplogrep.uibk.ac.at/blog/rcrs-vs-rsrs-vs-hg19/
#	Yoruba		Cambridge
#	1-309		Same
#	310			Deleted
#	311-316		Subtract One
#	317			Deleted
#	318-3108	Subtract Two
#	Deleted		Placeholder at 3107
#	3109-16194	Subtract One
#	16195		Deleted
#	16196-16571	Subtract Two
#-----------------------------------------------------##
yri_rcrs <- Vectorize(function(x){
  if (x < 310) {
    return(x)
  } else if (x == 310) {
    return(NA)
  } else if (x < 317) {
    return(x - 1)
  } else if (x == 317) {
    return(NA)
  } else if (x < 3109) {
    return(x - 2)
  } else if (x < 16195) {
    return(x - 1)
  } else if (x == 16195) {
    return(NA)
  } else if (x < 16571) {
    return(x - 2)
  } else {
    return(x)
  }
})

## Funciton to flip nuclotides
nucleotide_flip <- function(x){
  dplyr::recode(x, "A" = "T", "T" = "A", "G" = "C", "C" = "G", .default = x)
}

##  Read in reference
ref <- read_tsv(reference, col_types = "cicciid",
  col_names = c("CHROM", "POS", "Ref", "Alt", "AC", "AN", "AF"))

## Read in sample bim file
bim <- read_table2(bimfile, col_types = "iciicc",
  col_names = c("CHROM", "SNP", "cm", "POS", "A1", "A2")) %>%
  filter(CHROM == 26) %>%
  arrange(POS)

## Convert from YRI to rCRS
bim <- bim %>%
  mutate(POSrcrs = yri_rcrs(POS)) %>%
  filter(!is.na(POSrcrs))

## Join bim to referenece SNPs
pos_join <- left_join(bim, select(ref, -CHROM), by = "POS")
rcrs_join <- left_join(bim, select(ref, -CHROM), by = c("POSrcrs" = "POS"))

## Percentage of SNPs alligning to reference panel:
#   after YRI to rCRS conversion
perc.rcrs <- sum(!is.na(rcrs_join$AF)) / nrow(rcrs_join)
#   prior YRI to rCRS conversion
perc.pos <- sum(!is.na(pos_join$AF)) / nrow(pos_join)

check_aln <- Vectorize(function(ret_allele, oth_allele, ref) {
  if (is.na(ref)) {
  	return(stringr::str_to_lower(ret_allele))
  }
  else if (ret_allele == ref | oth_allele == ref) {
    return(ret_allele)
  } else {
    return(NA)
  }
})

flip <- . %>%
  mutate(noflip = !is.na(Ref) & (A1 == Ref | A2 == Ref)) %>%
  mutate(A1_flip = ifelse(noflip, A1, nucleotide_flip(A1))) %>%
  mutate(A2_flip = ifelse(noflip, A2, nucleotide_flip(A2))) %>%
  mutate(chk_a1 = check_aln(A1_flip, A2_flip, Ref)) %>%
  mutate(chk_a2 = check_aln(A2_flip, A1_flip, Ref))

strng <- "Percentage of SNPs aligining to reference panel %s to rCRS conversion: %s"
message(sprintf(strng, "prior", perc.pos))
message(sprintf(strng, "after", perc.rcrs))
if (perc.rcrs > perc.pos) {
  message("Sample is likely aligned to YRI - lifting over and aligning")
  cat("\n")
  out <- flip(rcrs_join) %>%
    select(-POS) %>%
    rename(POS = POSrcrs)
} else {
  message("Sample is likely aligned to rCRS - aligning original")
  cat("\n")
  out <- flip(pos_join)
}

out %>% filter(!is.na(Ref)) %>%
	summarise(nmiss_a1 = sum(is.na(chk_a1)),
            nmiss_a2 = sum(is.na(chk_a2))) %>%
	mutate(perc_a1 = nmiss_a1 / nrow(filter(out, !is.na(Ref))),
         perc_a2 = nmiss_a2 / nrow(filter(out, !is.na(Ref))))

strng <- "flipped strand due to allele mismatch at %i out of %i SNPs."
message(sprintf(strng, sum(out$flip, na.rm = T), nrow(out)))

out %>%
  select(CHROM, SNP, cm, POS, A1_flip, A2_flip)  %>%
  write_tsv(outfile, col_names = F)

sink()
