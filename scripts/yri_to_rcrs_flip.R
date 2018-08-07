#!/usr/bin/Rscript

suppressPackageStartupMessages(library(tidyverse))

#-----------------------------------------------------##
## 	Converting Yoruba (YRI - AF347015) SNP positions
##	to the revised cambrige reference sequence (rCRS)
##	https://www.mitomap.org/MITOMAP/YorubanConversion
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
yri_rcrs <- function(x){
  if(x < 311){
    x
  }else{
    if(x >= 311 & x < 318){
      x - 1
    }else {
      if(x >= 318 & x < 3109){
        x - 2
      }else {
        if(x >= 3109 & x < 16195){
          x - 1
        }else {
          if(x >= 16196 & x < 16571){
            x - 2
          }else {
            x
          }}}}}
}

## Funciton to flip nuclotides
nucleotide_flip <- function(x){dplyr::recode(x, 'A' = 'T', 'T' = 'A', 'G' = 'C', 'C' = 'G')}

## Arg parse
args = commandArgs(trailingOnly=TRUE)
reference <- args[1]
#reference <- '/Users/sheaandrews/Dropbox/src/MitoImputePrep/DerivedData/ReferencePanel/ReferenceSNPs.txt'
bimfile <- args[2]
#bimfile <- '/Users/sheaandrews/LOAD_minerva/dummy/Data/ADGC/ADGC_2018/ADGC_Phase1+Phase2/1000Genomes/ROSMAP1/RawGenotypes/ROSMAP1.bim'
#bimfile <- '/Users/sheaandrews/LOAD_minerva/dummy/Data/ADRC/GSA_2018/GSA/GSA/MSMD/ADRC/Cleaning/raj_ADRC.bim'
outfile <- args[3]

##  Read in reference
ref <- read_tsv(reference, col_names = F) %>% 
  rename(CHROM = X1, POS = X2, Ref = X3, Alt = X4, AC = X5, AN = X6, AF = X7)

## Read in sample bim file
bim <- read_table2(bimfile, col_names = F) %>% 
  rename(CHROM = X1, SNP = X2, cm = X3, POS = X4, A1 = X5, A2 = X6) %>% 
  filter(CHROM == 26) %>% arrange(POS)

## Convert from YRI to rCRS
bim$POSrcrs <- sapply(bim$POS, yri_rcrs)

## Join bim to referenece SNPs
pos_join <- left_join(bim, select(ref, -CHROM), by = 'POS') 
rcrs_join <- left_join(bim, select(ref, -CHROM), by = c('POSrcrs' = 'POS'))

## Percentage of SNPs alligning to reference panel: 
#   after YRI to rCRS conversion
perc.rcrs <- sum(!is.na(rcrs_join$AF)) / nrow(rcrs_join)
#   prior YRI to rCRS conversion
perc.pos <- sum(!is.na(pos_join$AF)) / nrow(pos_join)

if(perc.rcrs > perc.pos){
  cat('\nPercentage of SNPs aligining to reference panel prior to rCRS conversion:', perc.pos,
      '\nPercentage of SNPs aligining to reference panel after rCRS conversion:', perc.rcrs,
      '\nSample is likely aligned to YRI - lifting over \n')
}else{
  cat('\nPercentage of SNPs aligining to reference panel prior to rCRS conversion:', perc.pos,
      '\nPercentage of SNPs aligining to reference panel after rCRS conversion:', perc.rcrs,
      '\nSample is likely aligned to rCRS \n')
}

if(perc.rcrs > perc.pos){
  rcrs_join %>% 
    mutate(A1_flip = ifelse(A1 == Ref | A2 == Ref, A1, nucleotide_flip(A1))) %>% 
    mutate(A2_flip = ifelse(A1 == Ref | A2 == Ref, A2, nucleotide_flip(A2))) %>% 
    mutate(A1_flip = ifelse(is.na(A1_flip), A1, A1_flip)) %>% 
    mutate(A2_flip = ifelse(is.na(A2_flip), A2, A2_flip)) %>%
    select(CHROM, SNP, cm, POSrcrs, A1_flip, A2_flip)  %>%
    write_tsv(outfile, col_names = F)
} else {
  pos_join %>% 
    mutate(A1_flip = ifelse(A1 == Ref | A2 == Ref, A1, nucleotide_flip(A1))) %>% 
    mutate(A2_flip = ifelse(A1 == Ref | A2 == Ref, A2, nucleotide_flip(A2))) %>% 
    mutate(A1_flip = ifelse(is.na(A1_flip), A1, A1_flip)) %>% 
    mutate(A2_flip = ifelse(is.na(A2_flip), A2, A2_flip)) %>%
    select(CHROM, SNP, cm, POS, A1_flip, A2_flip) %>%
    write_tsv(outfile, col_names = F)
}


