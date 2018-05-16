require(phangorn)

infile = '/Volumes/MHS/Master_Alignment/FASTA/McInerney_Master_Alignment_Nov30_2017_ambig2missing.fasta'

aln = read.dna(infile, format = 'fasta')
base.freq(aln, T, T)

BC = read.csv('/Volumes/MHS/Other/Shea_Imputation/analyses/base_composition/McInerney_Master_Alignment_Nov30_2017_baseComposition.csv', header = T)

infile = '/Volumes/MHS/Other/Shea_Imputation/metadata/Haplo_Grouping/McInerney_Master_Alignment_Nov30_2017_HAPLOGROUPS_HAPLOFIND2.csv'

m = read.csv(infile, header = T)

hgs = m$Haplogroup
m
macro.hgs = substr(m$Haplogroup, 1, 1)

all.macro.hgs = unique(macro.hgs)
hgs.wanted = c('H','I','J','K', 'N','R', 'T','U','V','W','X')
nonEuro.macro = subset(all.macro.hgs, all.macro.hgs != hgs.wanted)

Euro.Only = subset(m, m$macrogroup %in% c('H', 'V', 'J', 'T', 'U', 'K', 'W', 'X', 'I', 'R', 'N', 'rCRS'))
Non.Euro = subset(m, !(m$Sequence.id %in% Euro.Only$Sequence.id))

write.table(Euro.Only$Sequence.id, '/Volumes/MHS/Other/Shea_Imputation/metadata/exclusion_lists/EuropeanOnlySeqs.txt', sep = '\t', quote = F, col.names = F, row.names = F)   
write.table(Non.Euro$Sequence.id, '/Volumes/MHS/Other/Shea_Imputation/metadata/exclusion_lists/NonEuropeanOnlySeqs.txt', sep = '\t', quote = F, col.names = F, row.names = F)

B = read.table('/Volumes/MHS/Other/Shea_Imputation/metadata/exclusion_lists/BIG_EXCLUSION_dups.txt', header = F, sep = '\t')

BdR = data.frame(unique(B))
write.table(BdR,'/Volumes/MHS/Other/Shea_Imputation/metadata/exclusion_lists/BIG_EXCLUSION_dupsRemoved.txt', row.names = F, col.names = F, quote = F, sep = '\t')

