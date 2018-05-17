Alignment procedure:

41,446 human mtDNA sequences were downloaded via Geneious v10.2.3 (Kearse et al., 2012) from NCBI's GenBank repository using the following search criteria:
	Homo [Organism] AND gene_in_mitochondrion[PROP] AND 16000:17000[SLEN] NOT pseudogene[All Fields]

The raw sequences were placed into a separate folder within Geneious: Shared Databases > geneious > TimMcInerney > GenBank_Seqs_Nov23_2017
Sequences were duplicated in batches of 1,000
Duplicates of those sequences were placed into Shared Databases > geneious > TimMcInerney > Shea_Imputation > Complete_Align
Sequences were de-circularised/changed to linear topology
Sequences were aligned to the Easteal Reference Alignment in batches of 1,000 using the MAFFT algorithm with default settings within Geneious
Once each batch alignment was completed, the resulting alignment was checked for so that orthologous positions matched up, and the 1,000 sequences were moved from the original folder to Shared Databases > geneious > TimMcInerney > Shea_Imputation > Complete_Align > DONE
If sequences were found to have starting positions inconsistent with the rCRS, then those sequences were re-circularised, and the starting position was moved to match the motif "GATCACAGGT" using the 'Change Residue Numbering' option.
Sequences that had to be re-numbered were:
	AM711903
	AM711904
	HUMMTA
If sequences could be aligned properly, they were removed, and the remaining sequences realigned
Removed sequences:
	AP008866*
	KC911619
	KF055290
	KF055291
	KF055292
	KF055293
	KF055294
	KF055295
	KF055296
	KF055297
	KF055298
	KF055299
	KF055300
	KF055301
	KF055302
	KF055303
	KF055304
	KF055305
	KF055306
	KF055307
	KF055308
	KF055309
	KF055310
	KF055311
	KF055312
	KF055313
	KF055314
	KF055315
	KF055316
	KF055317
	KF055318
	KF055319
	KF055320
	KF055321
	KF055322
	KF055323
	KF055324
	KF055325
	KF055326
	KF055327
	KF055328
	KF055329
	KF055330
	KF055331
	KF055332
	KF055333
The following sequences were removed due to being partial genomes:
	JX893366
	JX893367
The following sequences were removed due to being entirely made up of '?/ characters:
	CM001971
	CM003113
	CM003244
	CM003245
	CM003277
	CM003308
	CM003309
	CM003500
	CM003501
	CM003502
	CM003571
	CM003572
	CM003573
	CM003574
	CM003747
	CM003766

Once the alignment process was completed, the Easteal Reference Alignment sequences were removed, with the exception of the rCRS (NC_012920). The rCRS that was aligned during this novel procedure was removed, as the rCRS sequence in the Easteal Reference Alignment was given precedence.
The final alignment had a length of 41,400 sequences.

The alignment was then exported from Geneious 10.2.3 to fasta and fa.gz. The absolute path of the exported file is: /Volumes/MHS/Master_Alignment/McInerney_Master_Alignment_Nov30_2017.fasta

SNP-sites (Page et al., 2016) was used to convert the resulting fasta file into a haploid VCF file, using the following command:
	$ snp-sites -v -o /Volumes/MHS/Master_Alignment/McInerney_Master_Alignment_Nov30_2017.vcf.gz  /Volumes/MHS/Master_Alignment/McInerney_Master_Alignment_Nov30_2017.fasta
Plink v1.9 was then used to convert the haploid VCF file into a diploid  VCF file, consistent with the ADNI VCF, then into plink formatted .ped and .map files using the following command:
	$ plink1.9 --vcf /Volumes/MHS/Master_Alignment/McInerney_Master_Alignment_Nov30_2017.vcf.gz --recode vcf-iid --out /Volumes/MHS/Master_Alignment/McInerney_Master_Alignment_Nov30_2017_plinkVCF
	$ bcftools view -Oz -o /Volumes/MHS/Master_Alignment/McInerney_Master_Alignment_Nov30_2017_plinkVCF.vcf.gz /Volumes/MHS/Master_Alignment/McInerney_Master_Alignment_Nov30_2017_plinkVCF.vcf
	$ bcftools index /Volumes/MHS/Master_Alignment/McInerney_Master_Alignment_Nov30_2017_plinkVCF.vcf.gz
	$ plink1.9 --vcf /Volumes/MHS/Master_Alignment/McInerney_Master_Alignment_Nov30_2017_plinkVCF.vcf.gz --recode --out /Volumes/MHS/Master_Alignment/McInerney_Master_Alignment_Nov30_2017_plinkVCF

To assign haplogroups based off the PhyloTree nomenclature, HaploFind (Vianello et al., 2013) was used. As Haplofind couldn't take in the whole alignment, a custom script (https://git.nci.org.au/tm5170/timmci/blob/master/scripts/Other_Scripts/split_master.py) was written to break the alignment into bins of n=2500, then run these through Haplofind. The resulting split Haplofind documents where merged using the following command in R:
	require(readxl)
	require(plyr)
	
	wd = '/Volumes/MHS/Master_Alignment/SUB_HAPLOGROUPS/HAPLOFIND/'
	setwd(wd)
	files = dir(wd)
	out_file1 = '/Volumes/MHS/Master_Alignment/McInerney_Master_Alignment_Nov30_2017_HAPLOGROUPS_HAPLOFIND2.csv'
	out_file2 = '/Volumes/MHS/Master_Alignment/McInerney_Master_Alignment_Nov30_2017_HAPLOGROUPS_HAPLOFIND2_HGsummary.csv'
	out_file3 = '/Volumes/MHS/Master_Alignment/McInerney_Master_Alignment_Nov30_2017_HAPLOGROUPS_HAPLOFIND2_macroHGsummary.csv'
	
	for (file in 1:length(files)) {
	  if (file == 1) {
	    message(paste0('READING IN FIRST FILE: ', file))
	    df = data.frame(read_xlsx(files[1], sheet = 1, skip = 3))
	  }
	  else {
	    message(paste0('READING IN: ', file))
	    tmp = data.frame(read_xlsx(files[file], sheet = 1, skip = 3))
	    df = rbind.fill(df, tmp)
	    message('FILES BOUND TOGETHER')
	  }
	}
	
	hg.sum = table(df$Haplogroup)
	hg.sum = data.frame(hg.sum)
	colnames(hg.sum) = c('Haplogroup', 'Freq')
	
	df$macrogroup = substring(df$Haplogroup,1,1)
	macro.sum = data.frame(table(df$macrogroup))
	colnames(macro.sum) = c('Macro.Haplogroup', 'Freq')
	
	write.csv(df, out_file1, row.names = F, quote = F)
	write.csv(hg.sum, out_file2, row.names = F, quote = F)
	write.csv(macro.sum, out_file3, row.names = F, quote = F)
	
	hg.sum = read.csv(out_file2)

The Haplogroup assignments were then verified using the HiMC programme (Mitchell et al., 2017 [PRE-PRINT]), using the following script in R:
	install.packages("devtools")
	library("devtools")
	install_github("vserch/himc")
	library("HiMC")
	 
	mp = '/Volumes/MHS/Master_Alignment/sandbox/test.map'
	pd = '/Volumes/MHS/Master_Alignment/sandbox/test.ped'
	 
	snp_data_frame <- HiMC::generate_snp_data(mp,pd)
	data(nodes)
	classifications <- HiMC::getClassifications(snp_data_frame)
	ss2p = data.frame(classifications)
	ssOnly = read.csv('/Volumes/MHS/Master_Alignment/McInerney_Master_Alignment_Nov30_2017_HAPLOGROUPS_Hi_MC.csv', header = T)
	 
	for (i in 1:nrow(ss2p)) {
	  if (ss2p$haplogroup[i] == ssOnly$haplogroup[i]) {
	    ss2p$match.check = TRUE
	  }
	  else {
	    ss2p$match.check = FALSE
	  }
	}
	table(ss2p$match.check)
All 41400 haplogroupings matched, albeit with some terminal node assignments

I combined the exclusion list documents emailed by Chris Patterson, as listed below, into a single document: /Volumes/MHS/Other/Shea_Imputation/metadata/exclusion_lists/ALL.exclude.txt
Documents:
	alzheimers.exclude.txt
	ancient.exclude.txt
	disease.exclude.txt
	partial.exclude.txt
	species.exclude.txt

I removed accession extensions (i.e. ".1") to conform with the VCF file in TextWrangler
I removed duplicate sequence names in R using the following script:
	> require(dplyr)
	> master.seqs = read.table('/Volumes/MHS/Other/Shea_Imputation/metadata/McInerney_Master_Alignment_Nov30_2017_SAMPLE-IDs.txt', header = F)
exclude.seqs = read.table('/Volumes/MHS/Other/Shea_Imputation/metadata/exclusion_lists/ALL.exclude.txt', header = F)
	> exclude.seqs = unique(exclude.seqs)
	> exclude.seqs = data.frame(intersect(exclude.seqs$V1, master.seqs$V1))
	> write.table(exclude.seqs, '/Volumes/MHS/Other/Shea_Imputation/metadata/exclusion_lists/ALL.exclude_REMOVEdup.txt', row.names = F, quote = F, sep = '\t', col.names = F)
I removed these sequences from the master VCF using the following command:
	$ bcftools view -S ^/Volumes/MHS/Other/Shea_Imputation/metadata/exclusion_lists/ALL.exclude_REMOVEdup.txt -Oz -o /Volumes/MHS/Other/Shea_Imputation/data_files/VCF_files/McInerney_Master_Alignment_Nov30_2017_NOalz-anc-dis-par-spe.vcf.gz  /Volumes/MHS/Master_Alignment/McInerney_Master_Alignment_Nov30_2017.vcf.gz
	$ bcftools index /Volumes/MHS/Other/Shea_Imputation/data_files/VCF_files/McInerney_Master_Alignment_Nov30_2017_NOalz-anc-dis-par-spe.vcf.gz
In total, 5,179 sequences were removed, leaving 36,221 remaining

VCF file converted to diploid (to match ADNI format) using plink1.9:
	$ plink1.9 --vcf /Volumes/MHS/Other/Shea_Imputation/data_files/VCF_files/McInerney_Master_Alignment_Nov30_2017_NOalz-anc-dis-par-spe.vcf.gz --recode vcf-iid --out /Volumes/MHS/Other/Shea_Imputation/data_files/VCF_files/McInerney_Master_Alignment_Nov30_2017_NOalz-anc-dis-par-spe_PLINKVCF
	$ bcftools view -Oz -o /Volumes/MHS/Other/Shea_Imputation/data_files/VCF_files/McInerney_Master_Alignment_Nov30_2017_NOalz-anc-dis-par-spe_PLINKVCF.vcf.gz /Volumes/MHS/Other/Shea_Imputation/data_files/VCF_files/McInerney_Master_Alignment_Nov30_2017_NOalz-anc-dis-par-spe_PLINKVCF.vcf
	$ bcftools index /Volumes/MHS/Other/Shea_Imputation/data_files/VCF_files/McInerney_Master_Alignment_Nov30_2017_NOalz-anc-dis-par-spe_PLINKVCF.vcf.gz

In the master VCF, missing data was filled in with the reference allele, as to stop the MaCH algorithm imputing data on the reference:
	$ bcftools +missing2ref -Oz -o /Volumes/MHS/Other/Shea_Imputation/data_files/VCF_files/McInerney_Master_Alignment_Nov30_2017_NOalz-anc-dis-par-spe_PLINKVCF_missingFilled.vcf.gz /Volumes/MHS/Other/Shea_Imputation/data_files/VCF_files/McInerney_Master_Alignment_Nov30_2017_NOalz-anc-dis-par-spe_PLINKVCF.vcf.gz
	$ bcftools index /Volumes/MHS/Other/Shea_Imputation/data_files/VCF_files/McInerney_Master_Alignment_Nov30_2017_NOalz-anc-dis-par-spe_PLINKVCF_missingFilled.vcf.gz
78,796 missing alleles were filled in
The VCF was further subsetted to only include European haplogroups, which was achieved in R using the following script:
	> haps.file = '/Volumes/MHS/Master_Alignment/McInerney_Master_Alignment_Nov30_2017_HAPLOGROUPS_HAPLOFIND2.csv'
	> seq.haps = read.csv(haps.file, header = T)
	> accepted.hgs = c('H', 'I', 'J', 'K', 'T', 'U', 'V', 'W', 'X')
	> Euro.seqs = subset(seq.haps, substr(seq.haps$Haplogroup, 1, 1) == accepted.hgs)
	> Non.Euro.seqs = subset(seq.haps, substr(seq.haps$Haplogroup, 1, 1) != accepted.hgs)
	> write.csv(Euro.seqs, '/Volumes/MHS/Master_Alignment/McInerney_Master_Alignment_Nov30_2017_HAPLOGROUPS_HAPLOFIND2_EuroSubset.csv', row.names = F, quote = F)
	> write.table(Euro.seqs$Sequence.id, '/Volumes/MHS/Other/Shea_Imputation/metadata/exclusion_lists/EuropeanOnlySeqs.txt', row.names = F, quote = F, sep = '\t', col.names = F)
	> write.table(Non.Euro.seqs$Sequence.id, '/Volumes/MHS/Other/Shea_Imputation/metadata/exclusion_lists/NonEuropeanOnlySeqs.txt', row.names = F, quote = F, sep = '\t', col.names = F)
The ALL.exclude_REMOVEdup.txt and NonEuropeanOnlySeqs.txt were joined in TextWrangler and saved as ALL.exclude.NonEuro.txt
Duplicate sequences were removed in R by:
	> exclude.seqs = read.table('/Volumes/MHS/Other/Shea_Imputation/metadata/exclusion_lists/ALL.exclude.NonEuro.txt', header = F)
	> exclude.seqs = unique(exclude.seqs)
	> write.table(exclude.seqs, '/Volumes/MHS/Other/Shea_Imputation/metadata/exclusion_lists/ALL.exclude.NonEuro_REMOVEdup.txt', row.names = F, quote = F, sep = '\t', col.names = F)
The original Master VCF was then re-filtered on using the following command:
	$ bcftools view -S ^/Volumes/MHS/Other/Shea_Imputation/metadata/exclusion_lists/ALL.exclude.NonEuro_REMOVEdup.txt -Oz -o /Volumes/MHS/Other/Shea_Imputation/data_files/VCF_files/McInerney_Master_Alignment_Nov30_2017_EuroONLY_NOalz-anc-dis-par-spe.vcf.gz  /Volumes/MHS/Master_Alignment/McInerney_Master_Alignment_Nov30_2017.vcf.gz
	$ bcftools index /Volumes/MHS/Other/Shea_Imputation/data_files/VCF_files/McInerney_Master_Alignment_Nov30_2017_EuroONLY_NOalz-anc-dis-par-spe.vcf.gz
	$ bcftools +missing2ref -Oz -o /Volumes/MHS/Other/Shea_Imputation/data_files/VCF_files/McInerney_Master_Alignment_Nov30_2017_EuroONLY_NOalz-anc-dis-par-spe_missingfilled.vcf.gz /Volumes/MHS/Other/Shea_Imputation/data_files/VCF_files/McInerney_Master_Alignment_Nov30_2017_EuroONLY_NOalz-anc-dis-par-spe.vcf.gz
	$ bcftools index /Volumes/MHS/Other/Shea_Imputation/data_files/VCF_files/McInerney_Master_Alignment_Nov30_2017_EuroONLY_NOalz-anc-dis-par-spe_missingfilled.vcf.gz
These were converted to diploid VCF format using:
	$ plink1.9 --vcf /Volumes/MHS/Other/Shea_Imputation/data_files/VCF_files/McInerney_Master_Alignment_Nov30_2017_EuroONLY_NOalz-anc-dis-par-spe_missingfilled.vcf.gz --recode vcf-iid --out /Volumes/MHS/Other/Shea_Imputation/data_files/VCF_files/McInerney_Master_Alignment_Nov30_2017_EuroONLY_NOalz-anc-dis-par-spe_missingfilled_DIPLOID
	$ bcftools view -Oz -o /Volumes/MHS/Other/Shea_Imputation/data_files/VCF_files/McInerney_Master_Alignment_Nov30_2017_EuroONLY_NOalz-anc-dis-par-spe_missingfilled_DIPLOID.vcf.gz /Volumes/MHS/Other/Shea_Imputation/data_files/VCF_files/McInerney_Master_Alignment_Nov30_2017_EuroONLY_NOalz-anc-dis-par-spe_missingfilled_DIPLOID.vcf
	$ bcftools index /Volumes/MHS/Other/Shea_Imputation/data_files/VCF_files/McInerney_Master_Alignment_Nov30_2017_EuroONLY_NOalz-anc-dis-par-spe_missingfilled_DIPLOID.vcf.gz
The file was found to have wrong chromosome numbering, so this was changed to MT in TextWrangler
The ADNI files mtSNPs 750, 1438, 2706 and 4769 were manually flipped in TextWrangler so the REF/ALT alelles between the ADNI and McInerney VCF were the same

The filtered reference VCF was merged with the ADNI 1 dataset VCF using bcftools:
	$ bcftools merge -Oz -o /Volumes/MHS/Other/Shea_Imputation/data_files/VCF_files/MERGED_McInerneyNov30_2017_EuroONLY_NOalz-anc-dis-par-spe_missingfilled_ADNI1.vcf.gz /Volumes/MHS/Other/Shea_Imputation/data_files/VCF_files/McInerney_Master_Alignment_Nov30_2017_EuroONLY_NOalz-anc-dis-par-spe_missingfilled_DIPLOID.vcf.gz /Volumes/MHS/Other/Shea_Imputation/data_files/OLD/mito_snps_rcrs_timCORRECTED.vcf.gz
And then into .ped format using plink 1.9:
	$ plink1.9 --vcf /Volumes/MHS/Other/Shea_Imputation/data_files/VCF_files/MERGED_McInerneyNov30_2017_EuroONLY_NOalz-anc-dis-par-spe_missingfilled_ADNI1.vcf.gz --make-bed --const-fid --biallelic-only strict --out /Volumes/MHS/Other/Shea_Imputation/analyses/MACH_IMPUTATION/2018-01-14/MERGED_McInerneyNov30_2017_EuroONLY_NOalz-anc-dis-par-spe_missingfilled_ADNI1_MaCHrecode
The ped file was changed into a MaCH readable format using the following steps:
	Open in TextWrangler
	Change all instances of ' ' to '\t'
	Change all instance of '0\t-9' to 2
	Copy the entire second column and paste it over the first column (first column should be all 0s)
The .dat file was made by the following steps:
	Open the .map file in TextWrangler
	Save it as a new file
	Delete all but the last column
	Before the remaining text, place a 'M\t'
The files were moved onto the GDU cluster and run with the following commands:
	$ qsub -m e -M u5015730@anu.edu.au -b y -N impuMACH_2018-01-14 -o /home/easteallab/tim/logs/ -j y -q hugemem.q -l h_vmem=128g,virtual_free=127.5g /home/easteallab/software/executables/mach1 -d /home/easteallab/tim/Shea_Imputation/analysis/2018-01-14/MERGED_McInerneyNov30_2017_EuroONLY_NOalz-anc-dis-par-spe_missingfilled_ADNI1_MaCHrecode.dat -p /home/easteallab/tim/Shea_Imputation/analysis/2018-01-14/McInerney_ADNIexcludebar1.ped --compact -r 100 --prefix /home/easteallab/tim/Shea_Imputation/analysis/2018-01-14/pleaseWORK_bar1	# ^^^ That first one grabs the .rec and .erate files
	$ FILL IN

