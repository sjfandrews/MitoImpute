library(tidyverse)

##  Identifiy Samples with excisve missing data in the master reference alignment
##  and excessive non-ACTG annotations

# IUPAC notation
# https://en.wikipedia.org/wiki/Nucleic_acid_notation
# Symbol = Description
# A =	A
# C =	C		
# G =	G	
# T =	T
# U =	U
# W =	A, T
# S =	C, G	
# M = A, C		
# K =	G, T
# R =	A, G	
# Y =	C, T
# B =	C, G,	T	
# D =	A, G,	T
# H =	A, C, T
# V =	A, C, G	
# N =	A, C, G, T	
# Z	Zero

##  Import Master alignment
dat.fasta <- read_table("/Users/sheaandrews/Dropbox/Research/PostDoc/MitoWax/1_Raw_Data/Mito_reference/McInerney_Master_Alignment_Nov30_2017.fasta", col_names = F)

##  Spread sequence names and sequecents 
dat <- filter(dat.fasta, str_detect(X1, ">")); dat$X1 <- gsub('>', '', dat$X1)
dat$seq <- filter(dat.fasta, !str_detect(X1, ">"))$X1

## Count number of nucleiac acid notations in each sequence
dat <- dat %>% 
  mutate('A' = str_count(seq, "A")) %>% 
  mutate('C' = str_count(seq, "C")) %>%
  mutate('G' = str_count(seq, "G")) %>%
  mutate('T' = str_count(seq, "T")) %>%
  mutate('U' = str_count(seq, "U")) %>%
  mutate('W' = str_count(seq, "W")) %>%
  mutate('S' = str_count(seq, "S")) %>%
  mutate('M' = str_count(seq, "M")) %>%
  mutate('K' = str_count(seq, "K")) %>%
  mutate('R' = str_count(seq, "R")) %>%
  mutate('Y' = str_count(seq, "Y")) %>%
  mutate('B' = str_count(seq, "B")) %>%
  mutate('D' = str_count(seq, "D")) %>%
  mutate('H' = str_count(seq, "H")) %>%
  mutate('V' = str_count(seq, "V")) %>%
  mutate('N' = str_count(seq, "N")) %>%
  mutate('Z' = str_count(seq, "Z")) %>%
  mutate('gap' = str_count(seq, "-"))

dat$non_actg <- rowSums(select(dat, U, W, S, M, K, R, Y, B, D, H, V, Z))
dat$n_gap <- rowSums(select(dat, N, gap))
dat$total <- rowSums(select(dat, U, W, S, M, K, R, Y, B, D, H, V, N, Z, gap))
dat <- select(dat, X1, A, C, G, 'T', U, W, S, M, K, R, Y, B, D, H, V, N, Z, gap, non_actg, n_gap, total, seq)


## Cross-tab of samples with non-ACTG
count(dat, U); count(dat, W); count(dat, S); count(dat, M); count(dat, K); count(dat, R); count(dat, Y); count(dat, B); count(dat, D); count(dat, H); count(dat, V); count(dat, N); count(dat, Z)

##  Plots
count(dat, non_actg); ggplot(dat, aes(x = non_actg)) + geom_histogram() + scale_x_log10() + theme_bw()
count(dat, gap); ggplot(dat, aes(x = gap)) + geom_histogram() + scale_x_log10() + theme_bw()
count(dat, total); ggplot(dat, aes(x = total)) + geom_histogram() + scale_x_log10() + theme_bw()

##  Write out a list of samples with more than 5 gaps
reject.gap <- filter(dat, gap > 5)
write_tsv(select(reject.gap, X1), '~/Dropbox/Research/PostDoc/MitoWax/1_Raw_Data/Mito_reference/Exclusion_lists/sample_missing.txt', col_names = F)

##  Write out a list of samples with more than 0 non-ACTG
reject.non_actg <- filter(dat, non_actg > 0)
write_tsv(select(reject.non_actg, X1), '~/Dropbox/Research/PostDoc/MitoWax/1_Raw_Data/Mito_reference/Exclusion_lists/non_actg.txt', col_names = F)

##  Write out a list of samples with more than 1 non-ACTG
reject.non_actg <- filter(dat, N > 2)
reject.non_actg <- filter(dat, n_gap > 2)
write_tsv(select(reject.non_actg, X1), '~/Dropbox/Research/PostDoc/MitoWax/1_Raw_Data/Mito_reference/Exclusion_lists/non_actg.txt', col_names = F)

##  Write out a list of samples with more than 5 non-ACTG or gaps
reject.total <- filter(dat, total > 5)
##  Write out a list of samples with more than 5 gaps or more then 1 non-ACTG
reject.total <- filter(dat, n_gap > 5 | non_actg > 0) 
write_tsv(select(reject.total, X1), '~/Dropbox/Research/PostDoc/MitoWax/1_Raw_Data/Mito_reference/Exclusion_lists/non_actg_gaps.txt', col_names = F)






#############################################################
##  OLD CODE - use FASTA rather then VCF
#############################################################

##  From BCF tools, normalize master VCF to split multiallelic sities 
##  using SNP-sites to convert fasta to VCF
##  https://github.com/sanger-pathogens/snp-sites
#snp-sites -v -o \
#/Users/sheaandrews/Dropbox/Research/PostDoc/MitoWax/1_Raw_Data/Mito_reference/McInerney_Master_Alignment_Nov30_2017.vcf \
#/Users/sheaandrews/Dropbox/Research/PostDoc/MitoWax/1_Raw_Data/Mito_reference/McInerney_Master_Alignment_Nov30_2017.fasta

#bcftools annotate --rename-chrs /Users/sheaandrews/Dropbox/Research/PostDoc/MitoWax/1_Raw_Data/Mito_reference/rename_chr.txt \
#/Users/sheaandrews/Dropbox/Research/PostDoc/MitoWax/1_Raw_Data/Mito_reference/McInerney_Master_Alignment_Nov30_2017.vcf \
#--output /Users/sheaandrews/Dropbox/Research/PostDoc/MitoWax/1_Raw_Data/Mito_reference/MTReference_180215.vcf

## Split multiallelic sites
#bcftools norm -m -any -Ov -o \
#/Users/sheaandrews/Dropbox/Research/PostDoc/MitoWax/1_Raw_Data/Mito_reference/MTReference_multiallelic_180215.vcf \
#/Users/sheaandrews/Dropbox/Research/PostDoc/MitoWax/1_Raw_Data/Mito_reference/MTReference_180215.vcf

multiallelic.dat <- read_tsv('/Users/sheaandrews/Dropbox/Research/PostDoc/MitoWax/1_Raw_Data/Mito_reference/MTReference_multiallelic_180215.vcf', comment = '##')

##  Select SNPs with deletions
del.dat <- filter(multiallelic.dat, ALT == '*')

##  Number of deletions per sample
sample_missing <- colSums(select(del.dat, -`#CHROM`, -POS, -ID, -REF, -ALT, -QUAL, -FILTER, -INFO, -FORMAT))
sample_missing <- tibble(names(sample_missing), sample_missing)
sample_missing <- mutate(sample_missing, perc = (sample_missing/16569)*100)

arrange(sample_missing, -sample_missing)
count(sample_missing, sample_missing) %>% print(n = 30)
ggplot(count(test, sample_missing), aes(x = sample_missing, y = n))  + geom_point() + theme_bw() + scale_y_log10() 

ggplot(sample_missing, aes(x = perc)) + geom_histogram() + scale_y_log10() 

##  Write out a list of samples with more than 5 
reject <- filter(sample_missing, sample_missing > 5)
write_tsv(select(reject,`names(sample_missing)`), '~/Dropbox/Research/PostDoc/MitoWax/1_Raw_Data/Mito_reference/Exclusion_lists/sample_missing.txt', col_names = F)

## Missing at MT BP
pos_missing <- rowSums(select(del.dat, -`#CHROM`, -POS, -ID, -REF, -ALT, -QUAL, -FILTER, -INFO, -FORMAT))
pos_missing <- tibble(del.dat$POS, pos_missing)

pos_missing %>% filter(pos_missing > 1000) %>% print(n = 32)
ggplot(pos_missing, aes(x = `del.dat$POS`, y = pos_missing)) + geom_point() 
