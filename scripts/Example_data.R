library(HiMC); data(nodes)
library(caret)
library(tidyverse)
set.seed(333)

##===========================##
##  Example Reference Panel
##===========================##

## Readin Reference Panel .plink files
ref <- generate_snp_data('/Users/sheaandrews/Dropbox/src/MitoImpute/DerivedData/ReferencePanel/ReferencePanel.map', 
                         '/Users/sheaandrews/Dropbox/src/MitoImpute/DerivedData/ReferencePanel/ReferencePanel.ped')

## Assign haplogroups to Reference samples using Hi-MC
ref.classifications <- HiMC::getClassifications(ref)
ref.classifications %>% count(haplogroup) %>% print(n = Inf)

## Extract European only haplogroup samples
eur <- ref.classifications %>% 
  filter(haplogroup %in% c('H', 'H2a', 'H2a2a', 'HV', 'I', 'J', 'JT', 'K', 'K1', 'T', 'T1', 'U', 'U8b', 'V', 'W', 'X', 'X2')) %>% 
  as.tibble()

## Extract 10% of Eurorpean haplogroup samples, retaining the same proportion of haplgroups
red <- as.tibble(eur[createDataPartition(eur$haplogroup, p = .1, list = FALSE, times = 1),])

eur %>% count(haplogroup) %>% print(n = Inf)
red %>% count(haplogroup) %>% print(n = Inf)

## Write out reduced samples
write_tsv(select(red, Individual), '/Users/sheaandrews/Dropbox/src/MitoImpute/example/ExampleReferenceSamples.txt', col_names = F)

##===========================##
##  Example Sample Panel
##===========================##

EUR <- read_tsv('/Users/sheaandrews/Dropbox/src/MitoImpute/scripts/INFORMATION_LISTS/1000genomes_pops.txt') %>% 
  filter(Population %in% c('CEU', 'TSI', 'FIN', 'GBR', 'IBS')) 
  
## Write out EUR only samples
write_tsv(select(EUR, IID), '/Users/sheaandrews/Dropbox/src/MitoImpute/example/ExampleSamplePanel.txt', col_names = F)

  
