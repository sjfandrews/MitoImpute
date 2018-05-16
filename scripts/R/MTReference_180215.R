require(tidyverse)
## Creating list of sample names to be retained in MT references file

#reference.haps <- read_csv('~/Dropbox/Research/Collaboration/MitoWAX/Data/Genetics/ReferenceGT/Genious/McInerney_Master_Alignment_Nov30_2017_HAPLOGROUPS_HAPLOFIND2.csv', col_names = T)
#species.exclude <- read_tsv('~/Dropbox/Research/PostDoc/MitoWax/1_Raw_Data/Mito_reference/Exclusion_lists/species.exclude_MOD.txt', col_names = F)
#partial.exclude <- read_tsv('~/Dropbox/Research/PostDoc/MitoWax/1_Raw_Data/Mito_reference/Exclusion_lists/partial.exclude_MOD.txt', col_names = F)
#ancient.exclude <- read_tsv('~/Dropbox/Research/PostDoc/MitoWax/1_Raw_Data/Mito_reference/Exclusion_lists/ancient.exclude_MOD.txt', col_names = F)
#gaps_actg.exclude <- read_tsv('~/Dropbox/Research/PostDoc/MitoWax/1_Raw_Data/Mito_reference/Exclusion_lists/non_actg_gaps.txt', col_names = F)
##############################
reference.haps <- read_csv('/Volumes/TimMcInerney/Master_Alignment/HAPLOGROUPINGS/McInerney_Master_Alignment_Nov30_2017_HAPLOGROUPS_HAPLOFIND2.csv', col_names = T)
species.exclude <- read_tsv('/Volumes/TimMcInerney/Other/Shea_Imputation/metadata/exclusion_lists/species.exclude_MOD.txt', col_names = F)
partial.exclude <- read_tsv('/Volumes/TimMcInerney/Other/Shea_Imputation/metadata/exclusion_lists/partial.exclude_MOD.txt', col_names = F)
ancient.exclude <- read_tsv('/Volumes/TimMcInerney/Other/Shea_Imputation/metadata/exclusion_lists/ancient.exclude_MOD.txt', col_names = F)
gaps_actg.exclude <- read_tsv('/Volumes/TimMcInerney/Other/Shea_Imputation/metadata/exclusion_lists/gap7_missing5_dupsRemoved.txt', col_names = F)


#subsequently filtered to match the Major European haplogroups (H, V, J, T, U, K, W, X, I, R and N)
reference.haps$macro <- substr(reference.haps$Haplogroup, start = 1, stop = 1)
reference.filtered <- reference.haps %>% 
  filter(macro %in% c('H', 'V', 'J', 'T', 'U', 'K', 'W', 'X', 'I', 'R', 'N', 'rCRS')) %>% 
  select(Sequence.id, Haplogroup, macro) %>% 
  anti_join(species.exclude, by = c('Sequence.id' = 'X1')) %>% 
  anti_join(partial.exclude, by = c('Sequence.id' = 'X1')) %>% 
  anti_join(ancient.exclude, by = c('Sequence.id' = 'X1')) %>%
  anti_join(gaps_actg.exclude, by = c('Sequence.id' = 'X1'))  

write_tsv(select(reference.filtered,Sequence.id), '~/Dropbox/Research/PostDoc/MitoWax/1_Raw_Data/Mito_reference/Exclusion_lists/inclusion_list.txt', col_names = F)

sex <- select(reference.filtered,Sequence.id); sex$sex <- 'M'

write_tsv(sex, '~/Dropbox/Research/PostDoc/MitoWax/1_Raw_Data/Mito_reference/Exclusion_lists/h_samples.txt', col_names = F)

reference.haps$macro <- substr(reference.haps$Haplogroup, start = 1, stop = 1)
reference.filtered.all <- reference.haps %>% 
  #filter(macro %in% c('H', 'V', 'J', 'T', 'U', 'K', 'W', 'X', 'I', 'R', 'N', 'rCRS')) %>% 
  select(Sequence.id, Haplogroup, macro) %>% 
  anti_join(species.exclude, by = c('Sequence.id' = 'X1')) %>% 
  anti_join(partial.exclude, by = c('Sequence.id' = 'X1')) %>% 
  anti_join(ancient.exclude, by = c('Sequence.id' = 'X1')) %>%
  anti_join(gaps_actg.exclude, by = c('Sequence.id' = 'X1'))

write_tsv(select(reference.filtered.all,Sequence.id), '/Volumes/TimMcInerney/Other/Shea_Imputation/metadata/exclusion_lists/BIG_INCLUSION_EuroNonEuro.txt', col_names = F)


gm = read.table('/Volumes/TimMcInerney/Other/Shea_Imputation/metadata/exclusion_lists/gap7_missing5.txt', header = F, sep = '\t')
gm.u = unique(gm)
write.table(gm.u, '/Volumes/TimMcInerney/Other/Shea_Imputation/metadata/exclusion_lists/gap7_missing5_dupsRemoved.txt', col.names = F, row.names = F, quote = F)
