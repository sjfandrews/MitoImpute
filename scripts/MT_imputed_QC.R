library(tidyverse)
library(readxl)
library(HiMC); data(nodes)

##  Function
#source('/Users/u5015730/GitCode/timmci/scripts/Shea_Imputation/R/import_snps.R', chdir = TRUE)
source('~/GitCode/MitoImpute/scripts/import_snps.R', chdir = TRUE)

##  Function for calculating mathews correlation coefficent (MCC)
mccr <- function (act, pred) 
{
  TP <- sum(act %in% 1 & pred %in% 1)
  TN <- sum(act %in% 0 & pred %in% 0)
  FP <- sum(act %in% 0 & pred %in% 1)
  FN <- sum(act %in% 1 & pred %in% 0)
  denom <- as.double(TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)
  if (any((TP + FP) == 0, (TP + FN) == 0, (TN + FP) == 0, (TN + FN) == 0)) 
    denom <- 1
  mcc <- ((TP * TN) - (FP * FN))/sqrt(denom)
  return(mcc)
}

###---------------------------------------------------##
###                 Read in Data                      ##
###---------------------------------------------------##
## Read in info scores from impute2
info.score <- read_delim('/Volumes/TimMcInerney/Other/Shea_Imputation/data_files/IMPUTE2/mito_snps_rcrs_ed_nameMod_info', delim = " ")

##  Readin .ped files 
MTwgs <- generate_snp_data("/Volumes/TimMcInerney/Other/Shea_Imputation/data_files/plink_files/ADNI_reseq/adni_mito_genomes_180214.map",
                               "/Volumes/TimMcInerney/Other/Shea_Imputation/data_files/plink_files/ADNI_reseq/adni_mito_genomes_180214.ped")
MTadni <- generate_snp_data("/Volumes/TimMcInerney/Other/Shea_Imputation/data_files/IMPUTE2/imputed_plink/mito_snps_rcrs_ed_nameMod.map",
                            "/Volumes/TimMcInerney/Other/Shea_Imputation/data_files/IMPUTE2/imputed_plink/mito_snps_rcrs_ed_nameMod.ped")

# HiMC Haplotypes
#snp_df = HiMC::generate_snp_data("/Volumes/TimMcInerney/Other/Shea_Imputation/data_files/plink_files/ADNI_reseq/adni_mito_genomes_180214.map", "/Volumes/TimMcInerney/Other/Shea_Imputation/data_files/plink_files/ADNI_reseq/adni_mito_genomes_180214.ped")
#dat.hap = HiMC::getClassifications(snp_df)
#names(dat.hap) = c("SAMPLE_NUMBER", names(dat.hap[2:3]))

##  Read in vcf files
adni_imputed.raw <- read_tsv('/Volumes/TimMcInerney/Other/Shea_Imputation/data_files/IMPUTE2/imputed_plink/mito_snps_rcrs_ed_nameMod.vcf', comment = '##')
adni_wgs.raw <- read_tsv('/Volumes/TimMcInerney/Other/Shea_Imputation/data_files/VCF_files/ADNI_reseq/adni_mito_genomes_180214.vcf.gz', comment = '##')

##  Sample names for ADNI wgs samples and Phy-mer haplotypes
dat.hap <- read_excel('/Volumes/TimMcInerney/Other/Shea_Imputation/metadata/Haplo_Grouping/ADNI_Mitochondrial_Haplotypes_PhyMer.xlsx')

###---------------------------------------------------##
###             Munge Dataframes                      ##
###---------------------------------------------------##

##  subset imputed SNPs that are present in the MT WGS (mt wgs does not contain control region)
imputed.snps <- info.score %>% 
  filter(type == 0) %>% 
  filter(position <= 15693) %>% 
  select(position) %>% mutate(position = paste0('mt', position))

## Subset Sequenced SNPs 
wgs.snps <- adni_wgs.raw %>% 
  select(POS) %>% mutate(POS = paste0('mt', POS))

snps <- semi_join(imputed.snps, wgs.snps, by = c('position' = 'POS'))

##  Munge imputed dataframe
# Drop all columns but position and genotypes
adni_imputed <- adni_imputed.raw %>% 
  select(-`#CHROM`, -ID, -REF, -ALT, -QUAL, -FILTER, -INFO, -FORMAT) %>% 
  distinct(adni_imputed, POS, .keep_all = T) %>% 
  mutate(POS = paste0('mt', POS))
# Transpose dataframe
adni_imputed <- adni_imputed %>%
  gather(key = var_name, value = value, 2:ncol(adni_imputed)) %>% 
  spread_(key = names(adni_imputed)[1],value = 'value')
# recode sample anmes
adni_imputed <- adni_imputed %>% separate(var_name,c('SID', 'var_name') ,'_S_', ) %>% 
  mutate(var_name = as.integer(var_name)) %>% 
  select(-SID) %>% 
  arrange(var_name)
# substitute allels for NA, 0, 1 calls, change formate to interger
adni_imputed[adni_imputed == './.'] <- NA
adni_imputed[adni_imputed == '0/0'] <- 0
adni_imputed[adni_imputed == '1/1'] <- 1
adni_imputed <- adni_imputed %>% mutate_if(is.character, as.integer)

##  Munge ADNI wgs dataframe
# Drop all columns but position and genotypes
adni_wgs <- adni_wgs.raw %>% 
  select(-`#CHROM`, -ID, -REF, -ALT, -QUAL, -FILTER, -INFO, -FORMAT) %>% 
  distinct(adni_imputed, POS, .keep_all = T) %>% 
  mutate(POS = paste0('mt', POS))
# Transpose dataframe
adni_wgs <- adni_wgs %>%
  gather(key = var_name, value = value, 2:ncol(adni_wgs)) %>% 
  spread_(key = names(adni_wgs)[1],value = 'value')
##  correct sample names
adni_wgs <- left_join(adni_wgs, dat.hap, by = c('var_name' = 'SAMPLE_NUMBER')) %>% 
  separate(PATNO,c('SID', 'PATNO') ,'_S_', )  %>% 
  mutate(var_name = as.integer(PATNO)) %>%
  select(-SID) %>% select(-PATNO) %>%
  arrange(var_name)
# substitute allels for NA, 0, 1 calls, change formate to interger
adni_wgs[adni_wgs == './.'] <- NA
adni_wgs[adni_wgs == '0/0'] <- 0
adni_wgs[adni_wgs == '1/1'] <- 1
adni_wgs <- adni_wgs %>% mutate_if(is.character, as.integer)

##  obtain intersect of samples from ADNI1 and WGS; select imputed SNPs only
adni1_sample <- semi_join(adni_imputed, adni_wgs, by = 'var_name') %>% 
  select(var_name, snps$position)

wgs_sample <- adni_wgs %>% 
  semi_join(adni_imputed, by = 'var_name') %>% 
  select(var_name, snps$position)

###---------------------------------------------------##
###             Concordance statistics - SNPs         ##
###---------------------------------------------------##

## calculate MCC between imputed and typed SNPs
mccr.geno <- unlist(map2(wgs_sample[,2:ncol(wgs_sample)], adni1_sample[,2:ncol(adni1_sample)], function(a,b) mccr(a,b)))
##  Calculate concordance between imputed and typed SNPs
concodance.geno <- unlist(map2(adni1_sample[,2:ncol(adni1_sample)], wgs_sample[,2:ncol(wgs_sample)], function(a, b) sum(a == b, na.rm = T)/length(a)))
##  Calculate allele frequency of typed SNPs from WGS data
af <- unlist(map(wgs_sample[,2:ncol(wgs_sample)], function(a) sum(a, na.rm = T)/length(a)))

##  Dataframe for summary stats
summary.stats <- tibble(mtSNP = colnames(wgs_sample[,2:ncol(wgs_sample)]), 
                        pos = as.integer(gsub('mt', '', colnames(wgs_sample[,2:ncol(wgs_sample)]))),
                        af = af,
                        mcc = mccr.geno,
                        concodance = concodance.geno)
##  merge on info.score file
summary.stats <- left_join(summary.stats, select(info.score, c(-snp_id, -rs_id)), by = c('pos' = 'position'))

##  basic summary stats
summary.stats %>% count(af > 0.01)
summary.stats %>% count(mcc > 0.4)
summary.stats %>% count(concodance < 0.9)
summary.stats %>% count(info > 0.3); summary.stats %>% count(info > 0.5)

#---------------------------------------------------#
#             Visualization                         #

##  Plots by bp position - all SNPs
mcc = ggplot(summary.stats, aes(x = pos, y = mcc, colour = info > 0.4, size = af)) + geom_point() + theme_bw() + 
  labs(x = 'mtDNA position', y = 'MCC', title = 'Imputed mtSNP allele vs genotyped mtSNP allele') + 
  guides(colour=guide_legend(title="Impute2 Info \n Score > 0.4")) + 
  guides(size=guide_legend(title="Allele Frequency"))
ggsave('/Volumes/TimMcInerney/Other/Shea_Imputation/metadata/QC/tim_summary/MCC_imputed_v_genotyped.png', mcc)
ss = ggplot(summary.stats, aes(x = pos, y = concodance, colour = info > 0.3, size = af)) + geom_point() + theme_bw()
ggsave('/Volumes/TimMcInerney/Other/Shea_Imputation/metadata/QC/tim_summary/summary_stats_concordance.png', ss)

##  Plots by bp position - remove SNPs af < 0.01 & info < 0.3
ggplot(filter(summary.stats, af > 0.01, info > 0.3), aes(x = pos, y = mcc, colour = info > 0.4, size = af)) + geom_point() + theme_bw()
  ggplot(filter(summary.stats, af > 0.01, info > 0.3), aes(x = pos, y = concodance)) + geom_point()

##  Plots for typed SNPs
info_typed <- filter(info.score, type == 2)
ggplot(info_typed, aes(x = position, y = concord_type0, colour = info > 0.4, size = exp_freq_a1)) + geom_point() + theme_bw()
ggplot(info_typed, aes(x = position, y = r2_type0, colour = info > 0.4, size = exp_freq_a1)) + geom_point() + theme_bw()


###---------------------------------------------------##
###     Concordance statistics - MT Haplogroyps       ##
###---------------------------------------------------##

##  Filter SNPs that have an info score > 0.4
inc.snps <- info.score %>% filter(info > 0.4) %>% select(position)
keep.columns <- c('Family', 'Individual', 'Father', 'Mother', 'Sex', 'Phenotype', inc.snps$position)
MTadni <- MTadni[,keep.columns]

##  Hi-MC MT haplogroups classifications
MTadni.classifications <- HiMC::getClassifications(MTadni)
MTwgs.classifications <- HiMC::getClassifications(MTwgs)

adni_mt_haps <- MTwgs.classifications %>% 
  left_join(dat.hap, by = c('Individual' = 'SAMPLE_NUMBER')) %>% 
  rename(full_path_wgs = full_path, haplogroup_wgs = haplogroup) %>% 
  inner_join(MTadni.classifications, by = c('PATNO' = 'Individual')) %>% 
  rename(full_path_impute = full_path, haplogroup_impute = haplogroup) %>%
  select(PATNO, haplogroup_impute, haplogroup_wgs, HAPLOTYPE) %>% 
  as.tibble()

out <- adni_mt_haps %>% filter(haplogroup_impute != haplogroup_wgs) %>% print(n = Inf)
write_csv(out, '/Volumes/TimMcInerney/Other/Shea_Imputation/metadata/imputed/hap_comparison.csv')

adni_haps <- MTadni %>% left_join(MTwgs.classifications) %>% 
  as.tibble()

write_csv(adni_haps, '/Volumes/TimMcInerney/Other/Shea_Imputation/metadata/imputed/adni_haps.csv')











































