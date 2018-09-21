#!/usr/bin/Rscript

args = commandArgs(trailingOnly = TRUE) # Set arguments from the command line
wd = args[1]
sample = args[2]
info_cut = args[3]

## ===============================================## 
## Functions, Librays etc 
## ===============================================## 

library(tidyverse)
library(ggplot2)
library(ggforce)
library(HiMC); data(nodes)

## Function for reading in Map and Ped files
generate_snp_data_fixed <- function (map_file, ped_file) 
{
  map <- read.csv(map_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  header_row <- c("Family", "Individual", "Father", "Mother", 
                  "Sex", "Phenotype")
  snps = map[, 4]
  new_header = c(header_row, snps)
  ped <- read.csv(ped_file, sep = " ", header = FALSE, stringsAsFactors = FALSE, colClasses = 'character')
  range1 = seq(1, 6, by = 1)
  snp_data = data.frame(seq(1, nrow(ped), by = 1))
  for (i in range1) {
    snp_data[, i] = ped[, i]
  }
  range2 = seq(7, ncol(ped), by = 2)
  for (i in range2) {
    index = ((i - 7)/2) + 1
    snp_data[, index + 6] = paste(ped[, i], ped[, i + 1])
  }
  names(snp_data) <- new_header
  return(snp_data)
}


## Tibble of SNPs used in Hi-MC and the assocaited haplogroup
HiMC <- tibble(
  himc = 'yes',
  pos = as.numeric(c('10115', '1018', '10398', '10400', '10550', '11177', '11251', '11719', '11947', '12007', '12308', '12414', '12705', '13263', '13368', '13506', '13708', '13789', '14178', '14318', '1438', '14470', '14560', '14668', '14766', '14905', '15043', '15326', '15452', '15535', '16111', '16189', '16271', '16362', '16390', '16391', '16391', '1719', '1736', '2092', '3505', '3552', '3594', '4580', '4769', '4883', '4917', '4977', '5178', '5442', '6371', '7028', '825', '8251', '8414', '8468', '8703', '9042', '9055', '9347', '9950')),
  Haplogroup = c("L2", "L3", "K1", "M", "K", "B2", "JT", "R0", "W", "A2", "U", "N2", "R", "C", "T", "L2'3'4'6", "J", "L1", "L1", "C", "H2", "K", ".", "D4", "HV", "T", "N1a1b", "H2a2a", "JT", "B4b'd'e", "A2", "T1", "JT", "L4", "L2", "I", "I", "X2", "A", "D1", "W", "C", "L3'4", "V", "H2a", "M80'D", "T", "B2", "D ", "L0", "X ", "H", "L2'3'4'6", "N1a1b", "D4", "L2'3'4'6", "D2", "L0", "U8b", "L0", "B2"))


## ===============================================## 
## Specificy files  
## ===============================================## 
wd <- '~/LOAD_minerva/dummy/shea/Projects/3_mitoWAX/2_DerivedData/ADGC_180920_iter2'
setwd(wd)

sample <- 'TARC1_080213.ATGC'
typ_map <- paste0(sample, '/', sample, '_typedOnly.map') ## Map file for pre-imputed mtSNPs - sample_typedOnly.map
typ_ped <- paste0(sample, '/', sample, '_typedOnly.ped') ## Ped file for pre-imputed mtSNPs - sample_typedOnly.ped
imp_map <- paste0(sample, '/', 'Imputed_', sample, '.map') ## Map file for post-imputed mtSNPs - Imputed_{sample}.map
imp_ped <- paste0(sample, '/', 'Imputed_', sample, '.ped') ## Ped file for post-imputed mtSNPs - Imputed_{sample}.ped
imp_info <- paste0(sample, '/', sample, '_imputed_info')
info_cut <- 0.3 ## info score cut off 
#info_cut <- 0 ## info score cut off 

## ===============================================## 
## Read in / munge data
## ===============================================## 

##  Read in - Typed only
#     Plink .ped files 
samp_typ.ped <- generate_snp_data_fixed(typ_map, typ_ped)

##  Read in - Imputed 
#     Plink .ped files 
samp_imp.ped <- generate_snp_data_fixed(imp_map, imp_ped) 
samp_imp.ped <- samp_imp.ped[,-c(grep("\\<189\\>", colnames(samp_imp.ped)), grep("\\<16183\\>", colnames(samp_imp.ped)))]
samp_imp.info <- read_delim(imp_info, delim = " ")

## Combine info scores
samp_imp.info <- samp_imp.info %>% 
  mutate(info_comb = ifelse(info_type0 == -1, info,info_type0 )) %>% 
  left_join(HiMC, by = c('position' = 'pos')) %>%
  mutate(himc = ifelse(is.na(himc), 'no', himc))

## ===============================================## 
## INFO Score Plots
## ===============================================## 

# Info Score of imputed mtSNPs across the mitochondrial genome
ggplot(samp_imp.info, aes(x = position, y = info_comb, colour = as.factor(type))) + 
  geom_point() + theme_bw() + 
  labs(x = 'mtDNA position', y = 'Info Score') + geom_hline(yintercept = as.numeric(info_cut), colour = 'red', linetype = 2) +  
  guides(colour=guide_legend(title="SNP Type")) 
ggsave(filename = paste0(sample, '/', sample, '_Info.pdf'), device = 'pdf', width = 19.2, height = 6.87, units = 'cm')

# Info Score of imputed mtSNPs across the mitochondrial genome highlighting Hi-MC mtSNPs used in haplgroup assignment
ggplot(samp_imp.info, aes(x = position, y = info_comb, colour = as.factor(type), alpha = himc)) + geom_point() + theme_bw() + 
  labs(x = 'mtDNA position', y = 'Info Score') + geom_hline(yintercept = as.numeric(info_cut), colour = 'red', linetype = 2) + 
  guides(colour=guide_legend(title="SNP Type")) 
ggsave(filename = paste0(sample, '/', sample, '_InfoHiMC.pdf'), device = 'pdf', width = 19.2, height = 6.87, units = 'cm')

## ===============================================## 
## Haplogroup Classification
# [Hi-MC (Smieszek et al PeerJ 2018)](https://peerj.com/articles/5149/), a broad haplogroup classifier for European, African, and Native American mitochondrial haplogroups, is used to assign samples to macro haplogroups. Hi-MC was developed using a reduced mitochondrial phylogenetic tree targeting specific single nucleotide polymorphisms (SNPs) for genotyping.
## ===============================================## 

#   Typed Only
MTtyp.classifications <- HiMC::getClassifications(samp_typ.ped)

#   Imputed 
## remove sites with poor info scores (< 0.3)
rm.info <- filter(samp_imp.info, info > info_cut)
samp_imp.filtered <- samp_imp.ped[,colnames(samp_imp.ped) %in% c('Individual', rm.info$position)]

MTimp.classifications <- HiMC::getClassifications(samp_imp.filtered)

## Merge Typed and Imputed datasets
MT_haps <- MTtyp.classifications %>% 
  left_join(MTimp.classifications, by = 'Individual', suffix = c("_typ", "_imp")) %>% 
  as.tibble()

## Write out Haplogroups
out <- full_join(MTimp.classifications, samp_imp.filtered, by = 'Individual')
write_tsv(out, paste0(sample, '/', sample, '_Haplogroups.txt'))

## ===============================================## 
## Alluvial Diagram
## ===============================================## 

## Count pairs of haplogroups of imputed and WGS assignments
hap.match <- MT_haps %>%
  count(haplogroup_typ, haplogroup_imp) %>% 
  mutate(match = haplogroup_typ == haplogroup_imp) 

## Use ggforce to tidy data for geom_parallel Sets 
# Requires developmental version of ggforce
dat_ggforce <- hap.match  %>%
  gather_set_data(1:2) %>%       
  arrange(x,haplogroup_imp,desc(haplogroup_typ)) %>% 
  mutate(x = recode(x, 'haplogroup_imp' = '2', 'haplogroup_typ' = '0'))

#Haplogroup assignments pre- and post-imputation
ggplot(dat_ggforce, aes(x = x, id = id, split = y, value = n)) +
  geom_parallel_sets(aes(fill = match), alpha = 0.5, axis.width = 0.2) +
  geom_parallel_sets_labels(colour = 'black', angle = 0, size = 3) + 
  theme_classic() + theme(legend.position = 'bottom') + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text = element_text(size=12)) + 
  scale_x_discrete(labels=c("Typed", "Imputed")) + 
  labs(x = "Mitochondrial Haplogroups") + scale_fill_brewer(palette = 'Paired')
ggsave(filename = paste0(sample, '/', sample, '_HaplogroupMatch.pdf'), device = 'pdf', width = 19.2, height = 19, units = 'cm')
























