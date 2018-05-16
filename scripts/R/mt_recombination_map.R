MT_RefPanel <- read_tsv('/Users/sheaandrews/Dropbox/Research/PostDoc/MitoWax/1_Raw_Data/Mito_reference/MT_RefPanel_180221.vcf', comment = '##')
MT_RefPanel <- MT_RefPanel %>% separate(INFO, c('NS', 'AN', 'AF', 'MAF', 'AC', 'AC_Het', 'AC_Hom', 'AC_Hemi','HWE'), sep = ';')
MT_RefPanel <- MT_RefPanel %>% select(`#CHROM`, POS, ID, REF, ALT, QUAL ,FILTER, NS, AN, AF, MAF, AC)
MT_RefPanel$AF <- as.numeric(gsub('AF=', '', MT_RefPanel$AF))
MT_RefPanel$MAF <- as.numeric(gsub('MAF=', '', MT_RefPanel$MAF))

print(MT_RefPanel, n = Inf)
length(unique(MT_RefPanel$POS))

##  recombination map
# Fine-scale recombination map for the region to be analyzed. This file should have three columns: 
# 1 = physical position (in base pairs)
# 2 = recombination rate between current position and next position in map (in cM/Mb)
# 3 = genetic map position (in cM). 
# The file should also have a header line with an unbroken character string for each column (e.g., "position COMBINED_rate(cM/Mb) Genetic_Map(cM)"). 
# For the mitochondria, just set it to 0 and tbe bp for each snp

map.file <- tibble(position = MT_RefPanel$POS, COMBINED_rate.cM.Mb. = 0, Genetic_Map.cM. = MT_RefPanel$POS)
write_delim(map.file, '/Users/sheaandrews/Dropbox/Research/PostDoc/MitoWax/1_Raw_Data/Mito_reference/mt.map', delim = " ")
