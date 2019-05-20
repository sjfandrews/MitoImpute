#!/usr/bin/Rscript

args = commandArgs(trailingOnly = TRUE) # Set arguments from the command line
typ_map.path = args[1]
typ_ped.path = args[2]
imp_map.path = args[3]
imp_ped.path = args[4]
imp_info.path = args[5]
outdir = args[6]
sample = args[7]
info_cut = args[8]

## ===============================================## 
## Functions, Librays etc 
## ===============================================## 

library(tidyverse)
library(ggplot2)
library(ggfittext)
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
## Read in / munge data
## ===============================================## 

##  Read in - Typed only
#     Plink .ped files 
samp_typ.ped <- generate_snp_data_fixed(typ_map.path, typ_ped.path)

##  Read in - Imputed 
#     Plink .ped files 
samp_imp.ped <- generate_snp_data_fixed(imp_map.path, imp_ped.path) 
## Remove columns with duplicate names
samp_imp.ped <- samp_imp.ped[,-c(grep("\\<189\\>", colnames(samp_imp.ped)), 
                                 grep("\\<16183\\>", colnames(samp_imp.ped)))]
samp_imp.info <- read_delim(imp_info.path, delim = " ")

## Combine info scores
samp_imp.info <- samp_imp.info %>% 
  mutate(info_comb = ifelse(info_type0 == -1, info,info_type0 )) %>% 
  left_join(HiMC, by = c('position' = 'pos')) %>%
  mutate(himc = ifelse(is.na(himc), 'no', himc))

## ===============================================## 
## INFO Score Plots
## ===============================================## 

# Info Score of imputed mtSNPs across the mitochondrial genome
p1 <- ggplot(samp_imp.info, aes(x = position, y = info_comb, colour = as.factor(type))) + 
  geom_point() + theme_bw() + 
  labs(x = 'mtDNA position', y = 'Info Score', title = 'mtSNP INFO scores') + 
  geom_hline(yintercept = as.numeric(info_cut), colour = 'red', linetype = 2) +  
  guides(colour=guide_legend(title="SNP Type")) 
ggsave(plot = p1, filename = paste0(outdir, '/', sample, '_Info.png'), device = 'png', width = 19.2, height = 6.87, units = 'cm')

# Info Score of imputed mtSNPs across the mitochondrial genome highlighting Hi-MC mtSNPs used in haplgroup assignment
p2 <- ggplot(samp_imp.info, aes(x = position, y = info_comb, colour = as.factor(type), alpha = himc)) + geom_point() + theme_bw() + 
  labs(x = 'mtDNA position', y = 'Info Score', title = 'mtSNP INFO scores with Haplogroup defining SNPs highlighted') + 
  geom_hline(yintercept = as.numeric(info_cut), colour = 'red', linetype = 2) + 
  guides(colour=guide_legend(title="SNP Type")) 
ggsave(plot = p2, filename = paste0(outdir, '/', sample, '_InfoHiMC.png'), device = 'png', width = 19.2, height = 6.87, units = 'cm')

## ===============================================## 
## Haplogroup Classification
# Hi-MC (Smieszek et al PeerJ 2018), a broad haplogroup classifier for European, African, and Native American
# mitochondrial haplogroups, is used to assign samples to macro haplogroups. Hi-MC was developed using a reduced mitochondrial phylogenetic tree
# targeting specific single nucleotide polymorphisms (SNPs) for genotyping.
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
  as_tibble()

## Write out Haplogroups
out <- full_join(MTimp.classifications, samp_imp.filtered, by = 'Individual')
write_tsv(out, paste0(outdir, '/', sample, '_Haplogroups.txt'))

## ===============================================## 
## Confusion Matrix
## ===============================================## 

## Count pairs of haplogroups of imputed and WGS assignments
hap.match <- MT_haps %>%
  count(haplogroup_typ, haplogroup_imp) %>% 
  mutate(match = haplogroup_typ == haplogroup_imp) 

#Haplogroup assignments pre- and post-imputation
dat.p <- expand.grid(c(hap.match$haplogroup_imp, hap.match$haplogroup_typ), 
                     c(hap.match$haplogroup_imp, hap.match$haplogroup_typ)) %>% 
  rename(haplogroup_typ = Var1, haplogroup_imp = Var2) %>%
  left_join(hap.match) %>% 
  as_tibble() %>% 
  mutate(haplogroup_imp = as_factor(haplogroup_imp), haplogroup_typ = as_factor(haplogroup_typ)) 

p3 <- ggplot(dat.p, aes(x = haplogroup_imp, y = haplogroup_typ, fill = n, label = n)) + 
  geom_tile() + 
  geom_fit_text() +
  geom_vline(xintercept=seq(0.5, 50.5, 1),color="white") +
  geom_hline(yintercept=seq(0.5, 50.5, 1),color="white") +
  scale_fill_gradient(low = 'white', high = 'steelblue', na.value = 'grey90') + 
  theme_bw() + coord_equal() + 
  labs(x = 'Haplgroups assigned using \nimputed data', y = 'Haplgroups assignedusing \nun-imputed data')

ggsave(plot = p3, filename = paste0(outdir, '/', sample, '_HaplogroupMatch.png'), device = 'png', width = 20, height = 20, units = 'cm')
























