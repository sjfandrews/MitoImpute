require(phangorn)

harmonic.mean = function (x, na.rm = TRUE, zero = TRUE){
  # From package 'psych' available at: https://cran.r-project.org/web/packages/psych/index.html
  if (!zero) {
    x[x == 0] <- NA
  }
  if (is.null(nrow(x))) {
    1/mean(1/x, na.rm = na.rm)
  }
  else {
    1/(apply(1/x, 2, mean, na.rm = na.rm))
  }
}

infile = '/Volumes/TimMcInerney/Master_Alignment/McInerney_Master_Alignment_Nov30_2017.fasta'
hgfile = '/Volumes/TimMcInerney/Master_Alignment/McInerney_Master_Alignment_Nov30_2017_HAPLOGROUPS_HAPLOFIND2.csv'
aln = read.dna(infile, format = 'fasta')
hg = read.csv(hgfile, header = T)
hg = hg[,1:2]

Index = 1:nrow(aln)
df = data.frame(Index)
df$seq = row.names(aln)

for (i in 1:nrow(aln)) {
  if ((i %% 1000) == 0) {
    print(i)
  }
  #df$haplogroup[i] = as.character(hg$Haplogroup[match(df$seq[i], hg$Sequence.id)])
  #if (df$haplogroup[i] == 'rCRS') {
  #  df$macro[i] = 'rCRS'
  #}
  #if (df$haplogroup[i] != 'rCRS') {
  #  df$macro[i] = substr(df$haplogroup[i], 1, 1)
  #}
  tb = base.freq(aln[i,], freq = T, all = T)
  df$a[i] = as.numeric(tb[1])
  df$c[i] = as.numeric(tb[2])
  df$g[i] = as.numeric(tb[3])
  df$t[i] = as.numeric(tb[4])
  df$r[i] = as.numeric(tb[5])
  df$m[i] = as.numeric(tb[6])
  df$w[i] = as.numeric(tb[7])
  df$s[i] = as.numeric(tb[8])
  df$k[i] = as.numeric(tb[9])
  df$y[i] = as.numeric(tb[10])
  df$v[i] = as.numeric(tb[11])
  df$h[i] = as.numeric(tb[12])
  df$d[i] = as.numeric(tb[13])
  df$b[i] = as.numeric(tb[14])
  df$n[i] = as.numeric(tb[15])
  df$gap[i] = as.numeric(tb[16])
  df$missing[i] = as.numeric(tb[17])
  #df$sum[i] = sum(as.numeric(df[i,5:21]))
}
head(df)

df.missing = subset(df, df$n > 5)
df.gap = subset(df, df$gap > 7)

write.table(df.missing$seq, '/Volumes/MHS/Other/Shea_Imputation/metadata/exclusion_lists/missing5.txt', sep = '\t', quote = F, col.names = F, row.names = F)   
write.table(df.gap$seq, '/Volumes/MHS/Other/Shea_Imputation/metadata/exclusion_lists/gap7.txt', sep = '\t', quote = F, col.names = F, row.names = F)

write.csv(df, '/Volumes/MHS/Other/Shea_Imputation/analyses/base_composition/McInerney_Master_Alignment_Nov30_2017_ambig2missing_baseComposition.csv', quote = F, row.names = F)

df$m = as.numeric(df$m)
df$h = as.numeric(df$h)
dfSum = summary(df[,5:22])
dfSum

write.csv(dfSum, '/Volumes/TimMcInerney/Other/Shea_Imputation/analyses/base_composition/McInerney_Master_Alignment_Nov30_2017_baseCompositionSummary.csv', quote = F, row.names = F)

macroHGs = unique(df$macro)
macroHGs.df = data.frame(macroHGs)
for (i in 1:nrow(macroHGs.df)) {
  tmp.df = subset(df, df$macro == macroHGs.df$macroHGs[i])
  macroHGs.df$obs[i] = nrow(tmp.df)
  macroHGs.df$gap.min[i] = min(tmp.df$gap)
  macroHGs.df$gap.1stQ[i] = quantile(tmp.df$gap, 0.25)
  macroHGs.df$gap.mean[i] = mean(tmp.df$gap)
  macroHGs.df$gap.median[i] = median(tmp.df$gap)
  macroHGs.df$gap.Harm.mean[i] = harmonic.mean(tmp.df$gap)
  macroHGs.df$gap.3rdQ[i] = quantile(tmp.df$gap, 0.75)
  macroHGs.df$gap.max[i] = max(tmp.df$gap)
  #
  macroHGs.df$n.min[i] = min(tmp.df$n)
  macroHGs.df$n.1stQ[i] = quantile(tmp.df$n, 0.25)
  macroHGs.df$n.mean[i] = mean(tmp.df$n)
  macroHGs.df$n.median[i] = median(tmp.df$n)
  macroHGs.df$n.Harm.mean[i] = harmonic.mean(tmp.df$n)
  macroHGs.df$n.3rdQ[i] = quantile(tmp.df$n, 0.75)
  macroHGs.df$n.max[i] = max(tmp.df$n)
}
macroHGs.df

write.csv(macroHGs.df, '/Volumes/TimMcInerney/Other/Shea_Imputation/analyses/base_composition/McInerney_Master_Alignment_Nov30_2017_baseComposition_HaplogroupSummary.csv', quote = F, row.names = F)

CutOff = 16569 * 0.05

df.5pcCut = subset(df, df$n < CutOff)
df.5pcCut.ex = subset(df, df$n > CutOff)

macroHGs.5pcCut = unique(df.5pcCut$macro)
macroHGs.df.5pcCut = data.frame(macroHGs.5pcCut)
for (i in 1:nrow(macroHGs.df.5pcCut)) {
  tmp.df = subset(df.5pcCut, df.5pcCut$macro == macroHGs.df.5pcCut$macroHGs[i])
  macroHGs.df.5pcCut$obs[i] = nrow(tmp.df)
  macroHGs.df.5pcCut$gap.min[i] = min(tmp.df$gap)
  macroHGs.df.5pcCut$gap.1stQ[i] = quantile(tmp.df$gap, 0.25)
  macroHGs.df.5pcCut$gap.mean[i] = mean(tmp.df$gap)
  macroHGs.df.5pcCut$gap.median[i] = median(tmp.df$gap)
  macroHGs.df.5pcCut$gap.Harm.mean[i] = harmonic.mean(tmp.df$gap)
  macroHGs.df.5pcCut$gap.3rdQ[i] = quantile(tmp.df$gap, 0.75)
  macroHGs.df.5pcCut$gap.max[i] = max(tmp.df$gap)
  #
  macroHGs.df.5pcCut$n.min[i] = min(tmp.df$n)
  macroHGs.df.5pcCut$n.1stQ[i] = quantile(tmp.df$n, 0.25)
  macroHGs.df.5pcCut$n.mean[i] = mean(tmp.df$n)
  macroHGs.df.5pcCut$n.median[i] = median(tmp.df$n)
  macroHGs.df.5pcCut$n.Harm.mean[i] = harmonic.mean(tmp.df$n)
  macroHGs.df.5pcCut$n.3rdQ[i] = quantile(tmp.df$n, 0.75)
  macroHGs.df.5pcCut$n.max[i] = max(tmp.df$n)
}
macroHGs.df.5pcCut

write.csv(macroHGs.df.5pcCut, '/Volumes/TimMcInerney/Other/Shea_Imputation/analyses/base_composition/McInerney_Master_Alignment_Nov30_2017_baseComposition_HaplogroupSummary_5pcGapMissingCO.csv', quote = F, row.names = F)


