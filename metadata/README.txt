41,846 human mtDNA samples were downloaded from NCBI/GenBank using the search term “Homo [Organism] AND gene_in_mitochondrion[PROP] AND 14000:19000[SLEN] NOT pseudogene[All Fields]”.

The ‘Send to:’ button was used to do this, using the options:
	* ’Complete record’
	* ’Choose Destination: File’
	* Format: Fasta
	* Sort by: Accession number
The file was renamed Homo_sapiens_mtDNA.fasta.txt and transferred to the GDU cluster in the following directory: /home/easteallab/tim/Shea_Imputation/data

The sequences were aligned using MUSCLE v3.8.31 using the following command:
$ qsub -m e -M u5015730@anu.edu.au -b y -N H.sapiens_mtDNA_aln -o /home/easteallab/tim/logs/ -j y -pe threads 24 muscle -in /home/easteallab/tim/Shea_Imputation/data/Homo_sapiens_mtDNA.fasta.txt -out /home/easteallab/tim/Shea_Imputation/data/Homo_sapiens_mtDNA_aligned.fasta

While the sequences were aligning, tests were carried out using Prof. Simon Easteal’s master alignment (2011) (henceforth, the Easteal alignment). This alignment contains 7,747 human mtDNA sequences, aligned to the revised Cambridge Reference Sequence.

The Easteal alignment was split into two fasta files for upload into HaploGrep 2.0. Each alignment contained 3,874 and 3,873 sequences each. These were both uploaded into HaploGrep (https://haplogrep.uibk.ac.at), then exported to .vcf:
	$ /Volumes/mhs/Other/Shea_Imputation/hsapiensCRS7k_HALVED1.vcf
	$ /Volumes/mhs/Other/Shea_Imputation/hsapiensCRS7k_HALVED2.vcf
Using bcftools 1.6, the vcf files were compressed using the following commands:
	$ bcftools view /Volumes/mhs/Other/Shea_Imputation/hsapiensCRS7k_HALVED1.vcf -Oz -o /Volumes/mhs/Other/Shea_Imputation/hsapiensCRS7k_HALVED1.vcf.gz
	$ bcftools view /Volumes/mhs/Other/Shea_Imputation/hsapiensCRS7k_HALVED1.vcf -Oz -o /Volumes/mhs/Other/Shea_Imputation/hsapiensCRS7k_HALVED2.vcf.gz
The vcf.gz files were indexed using the following command:
	$ bcftools index /Volumes/mhs/Other/Shea_Imputation/hsapiensCRS7k_HALVED1.vcf.gz 
	$ bcftools index /Volumes/mhs/Other/Shea_Imputation/hsapiensCRS7k_HALVED2.vcf.gz
The two vcf.gz files were then merged with the following command:
	$ bcftools merge /Volumes/mhs/Other/Shea_Imputation/hsapiensCRS7k_HALVED1.vcf.gz /Volumes/mhs/Other/Shea_Imputation/hsapiensCRS7k_HALVED2.vcf.gz -Oz -o /Volumes/mhs/Other/Shea_Imputation/hsapiensCRS7k_MERGED.vcf.gz
	$ bcftools index /Volumes/mhs/Other/Shea_Imputation/hsapiensCRS7k_MERGED.vcf.gz

A merged dataset of the novel SNPs were merged using the following command:
	$ bcftools merge /Volumes/mhs/Other/Shea_Imputation/hsapiensCRS7k_MERGED.vcf.gz /Volumes/mhs/Other/Shea_Imputation/mitosnps_adnimerge_rcrs.vcf.gz -Oz -o /Volumes/mhs/Other/Shea_Imputation/SHEA_MERGED.vcf.gz
	$ bcftools index /Volumes/mhs/Other/Shea_Imputation/SHEA_MERGED.vcf.gz

The vcf files were converted to PLINK formats using the following command in PLINK 1.9:
	$ plink1.9 --vcf /Volumes/mhs/Other/Shea_Imputation/hsapiensCRS7k_MERGED.vcf.gz --make-bed --double-id --out hsapiensCRS7k_MERGED
The vcf file of novel samples in this study were converted to PLINK formats using the following command:
	$ plink1.9 --vcf /Volumes/mhs/Other/Shea_Imputation/mitosnps_adnimerge_rcrs.vcf.gz --make-bed --double-id --out mitosnps_adnimerge_rcrs

The combined vcf file (Easteal + novel) was converted to PLINK format using the following command:
	$ plink1.9 --vcf /Volumes/mhs/Other/Shea_Imputation/SHEA_MERGED.vcf.gz --make-bed --double-id --out /Volumes/mhs/Other/Shea_Imputation/SHEA_MERGED

As ‘—proxy-assoc all’ has been retired in PLINK 1.9, PLINK 1.07 was used to impute the SNPs using the following command:
	$ plink1.07 --noweb --bfile /Volumes/mhs/Other/Shea_Imputation/SHEA_MERGED --proxy-impute all --make-bed --out /Volumes/mhs/Other/Shea_Imputation/SHEA_IMPUTED
The PLINK formatted imputed files were converted to VCF using the following command:
	$ plink1.9 --bfile /Volumes/mhs/Other/Shea_Imputation/SHEA_IMPUTED --recode vcf --out /Volumes/mhs/Other/Shea_Imputation/SHEA_IMPUTED
	$ bcftools view /Volumes/mhs/Other/Shea_Imputation/SHEA_IMPUTED.vcf -Oz -o /Volumes/mhs/Other/Shea_Imputation/SHEA_IMPUTED.vcf.gz
	$ bcftools index /Volumes/mhs/Other/Shea_Imputation/SHEA_IMPUTED.vcf.gz

As there were file conversion issues with the ADNI and Easteal VCFs, a custom script was written that converts mtDNA VCF files into FASTA files, with the rCRS at the top of the file. In turning the VCF into FASTA, where a position has not been called, it is assigned N. This results in an alignment with n+1 samples, with a length of 16,569 nucelotides. This command was run:
	$ python /Users/u5015730/Documents/GitCode/timmci/scripts/Other_Scripts/mtDNA_VCF2FASTA.py -i /Volumes/mhs/Other/Shea_Imputation/mitosnps_adnimerge_rcrs.vcf.gz -o /Volumes/mhs/Other/Shea_Imputation/mitosnps_adnimerge_rcrs2.fasta

The resulting FASTA file was aligned with the Easteal alignment in Geneious v10.2.3, using the MAFFT algorithm and default settings. The resulting alignment was divided into 2 alignments of roughly equal length, as HaploGrep cannot support files >100MB.
The two alignments were run through HaploGrep and exported as VCF files, then combined in bcftools using the following commands:
	$  bcftools view /Volumes/mhs/Other/Shea_Imputation/Shea_withEASTEAL_H1.vcf -Oz -o /Volumes/mhs/Other/Shea_Imputation/Shea_withEASTEAL_H1.vcf.gz
	$ bcftools index /Volumes/mhs/Other/Shea_Imputation/Shea_withEASTEAL_H2.vcf.gz
	$ bcftools view /Volumes/mhs/Other/Shea_Imputation/Shea_withEASTEAL_H2.vcf -Oz -o /Volumes/mhs/Other/Shea_Imputation/Shea_withEASTEAL_H2.vcf.gz
	$ bcftools index /Volumes/mhs/Other/Shea_Imputation/Shea_withEASTEAL_H2.vcf.gz
	$ bcftools merge /Volumes/mhs/Other/Shea_Imputation/Shea_withEASTEAL_H1.vcf.gz /Volumes/mhs/Other/Shea_Imputation/Shea_withEASTEAL_H2.vcf.gz -Oz -o /Volumes/mhs/Other/Shea_Imputation/SHEA_MERGED_4.vcf.gz
	$ bcftools index /Volumes/mhs/Other/Shea_Imputation/SHEA_MERGED_4.vcf.gz
The resulting combined VCF file was run through plink 1.9 to generate .bed, .bim, and .fam files using the following command:
	$ plink1.9 --vcf /Volumes/mhs/Other/Shea_Imputation/SHEA_MERGED_4.vcf.gz --make-bed --out /Volumes/mhs/Other/Shea_Imputation/SHEA_MERGED_4
Haplotypes were then imputed using the following command in plink 1.07:
	$ plink1.07 --noweb --bfile /Volumes/mhs/Other/Shea_Imputation/SHEA_MERGED_4 --flip-scan --proxy-assoc --proxy-drop --make-bed --out /Volumes/mhs/Other/Shea_Imputation/SHEA_IMPUTED_4
	$ plink1.9 --bfile /Volumes/mhs/Other/Shea_Imputation/SHEA_IMPUTED_4 --recode vcf --out /Volumes/mhs/Other/Shea_Imputation/SHEA_IMPUTED_4p
	$ bcftools view /Volumes/mhs/Other/Shea_Imputation/SHEA_IMPUTED_4p.vcf -Oz -o /Volumes/mhs/Other/Shea_Imputation/SHEA_IMPUTED_4p.vcf.gz
	$ bcftools index /Volumes/mhs/Other/Shea_Imputation/SHEA_IMPUTED_4p.vcf.gz



As the above lines did not work, and it was subsequently discovered that PLINK’s imputation algorithm was inaccurate, the following procedure was carried out:

Merge ADNI and Easteal datasets (VCF files)
	$ bcftools view /Volumes/MHS/Other/Shea_Imputation/data_files/mito_snps_rcrs_CORRECTED.vcf -Oz -o /Volumes/MHS/Other/Shea_Imputation/data_files/mito_snps_rcrs_CORRECTED.vcf.gz
	$ bcftools index /Volumes/MHS/Other/Shea_Imputation/data_files/mito_snps_rcrs_CORRECTED.vcf.gz
	$ bcftools merge /Volumes/MHS/Other/Shea_Imputation/data_files/hsapiensCRS7k_MERGED.vcf.gz /Volumes/MHS/Other/Shea_Imputation/data_files/mito_snps_rcrs_CORRECTED.vcf.gz -Oz -o /Volumes/MHS/Other/Shea_Imputation/data_files/easteal_ADNI_merged.vcf.gz
	$ bcftools index /Volumes/MHS/Other/Shea_Imputation/data_files/easteal_ADNI_merged.vcf.gz

Convert file types to PLINK .ped formats:
	$ plink1.9 --vcf /Volumes/MHS/Other/Shea_Imputation/data_files/easteal_ADNI_merged.vcf.gz --make-bed --const-fid --out /Volumes/MHS/Other/Shea_Imputation/analyses/MACH_IMPUTATION/easteal_ADNI_merged
	$ plink1.9 --vcf /Volumes/MHS/Other/Shea_Imputation/data_files/easteal_ADNI_merged.vcf.gz --recode A --const-fid --out /Volumes/MHS/Other/Shea_Imputation/analyses/MACH_IMPUTATION/easteal_ADNI_merged

Copies of the .ped and .map files were make, with the suffix ‘_ORIGINAL’ coming before the file extension. This ensured that if the .map or .ped files needed to be modified, we could always return to the original

To convert the .map file to a .dat file, the .map file was loaded into R and converted using the following:
	> wd = '/Volumes/MHS/Other/Shea_Imputation/analyses/MACH_IMPUTATION/'
	> setwd(wd)
	> infile = 'easteal_ADNI_merged'

	> map.file = read.table(paste0(infile,'_ORIGINAL.map'), header = F, sep = '\t')
	> mod.map = data.frame(map.file)
	> Index = c(1:nrow(mod.map))
	> new.map = data.frame(Index)

	> for (x in new.map$Index) {
  		new.map$c1[x] = mod.map$V1[x]
  		temp = paste0('MT', x)
  		new.map$c2[x] = temp
  		new.map$c3[x] = mod.map$V3[x]
  		new.map$c4[x] = mod.map$V4[x]
	}
	> new.map$Index = NULL
	> head(new.map)
	> write.table(new.map, paste0(infile,'.map'), sep = '\t', col.names = F, row.names = F, quote = F)
	> dat.map = data.frame(Index)
	> for (x in dat.map$Index) {
  		dat.map$c1[x] = 'M'
  		dat.map$c2[x] = mod.map$V4[x]
	}
	> dat.map$Index = NULL
	> head(dat.map)
	> write.table(dat.map, paste0(infile,'.dat'), sep = '\t', col.names = F, row.names = F, quote = F)

The following step was carried out using MACH for pre-imputation quality control:
	$ mach1 -d /Volumes/MHS/Other/Shea_Imputation/analyses/MACH_IMPUTATION/easteal_ADNI_merged.dat -p /Volumes/MHS/Other/Shea_Imputation/analyses/MACH_IMPUTATION/easteal_ADNI_merged.ped --greedy --prefix /Volumes/MHS/Other/Shea_Imputation/analyses/MACH_IMPUTATION/step1

At this step the following error was thrown, showing that we needed to remove multi allelic sites:
	FATAL ERROR:
	This version of MaCH is designed for bi-allelic markers
	However, the following marker(s) have >2 alleles:
   		7394 12432 10664 9301 6401 2260 16488 10034 327 16265 6158 10843 2359 2779 11044 12702 
   		14905 14580 5606 931 12865 7852 494 15261 1375 7046 16490 8471 9293 4029 12400 8192 4943 
   		10149 13074 11710 14395 417 3438 12019 14793 7241 12621 9653 4736 14497 10739 3105 12375 
   		3924 13422 6806 10742 10361 11047 16131 7257 507 11929 517 4884 14211 1692 8469 10502 3696 
   		9299 5375 15061 10897 9287 15859 3327 6239 3637 4508 14370 12574 16518 16322 8496 16103 
   		9974 188 15752 3411 8781 12684 12026 7058 11329 14550 15927 11229 13285 8701 4580 7927 
   		15861 14272 9944 15331 10303 348 12437 5821 13884 15705 616 7705 13596 12964 3663 3843 3317 
   		1914 6722 6881 11090 16311 8206 4005 9125 345 11437 2766 13351 10181 4203 6022 1664 15316 
   		7595 3873 4733 16242 4161 5267 3392 9950 9275 13368 13387 12028 5436 3450 4158 10892 14894 
   		411 8233 11242 6779 10197 10499 3531 13934 12438 10355 6927 591 14161 4222 7785 6367 15154 
   		1978 9899 14551 2863 7444 3651 5189 6134 3615 8805 12040 11632 15479 13716 8631 9065 6649 
   
	2930 additional markers not listed

	Please remove or recode markers with more than 2 alleles

To remove multi allelic sites, the following commands were run:
	$ bcftools view --max-alleles 2 --exclude-types indels easteal_ADNI_merged.vcf.gz -Oz -o /Volumes/MHS/Other/Shea_Imputation/data_files/easteal_ADNI_merged_biallelic.vcf.gz
	$ bcftools index /Volumes/MHS/Other/Shea_Imputation/data_files/easteal_ADNI_merged_biallelic.vcf.gz

Convert file types to PLINK .ped formats:
	$ plink1.9 --vcf /Volumes/MHS/Other/Shea_Imputation/data_files/easteal_ADNI_merged_biallelic.vcf.gz --recode --const-fid --out /Volumes/MHS/Other/Shea_Imputation/analyses/MACH_IMPUTATION/easteal_ADNI_merged_biallelic
	$ plink1.9 --vcf /Volumes/MHS/Other/Shea_Imputation/data_files/easteal_ADNI_merged_biallelic.vcf.gz --recode A --const-fid --out /Volumes/MHS/Other/Shea_Imputation/analyses/MACH_IMPUTATION/easteal_ADNI_merged_biallelic

Copies of the .ped and .map files were make, with the suffix ‘_ORIGINAL’ coming before the file extension. This ensured that if the .map or .ped files needed to be modified, we could always return to the original

To convert the .map file to a .dat file, the .map file was loaded into R and converted using the following:
	> wd = '/Volumes/MHS/Other/Shea_Imputation/analyses/MACH_IMPUTATION/'
	> setwd(wd)
	> infile = 'easteal_ADNI_merged_biallelic'

	> map.file = read.table(paste0(infile,'_ORIGINAL.map'), header = F, sep = '\t')
	> mod.map = data.frame(map.file)
	> Index = c(1:nrow(mod.map))
	> new.map = data.frame(Index)

	> for (x in new.map$Index) {
  		new.map$c1[x] = mod.map$V1[x]
  		temp = paste0('MT', mod.map$V4[x])
  		new.map$c2[x] = temp
  		new.map$c3[x] = mod.map$V3[x]
  		new.map$c4[x] = mod.map$V4[x]
	}
	> new.map$Index = NULL
	> head(new.map)
	> write.table(new.map, paste0(infile,'.map'), sep = '\t', col.names = F, row.names = F, quote = F)
	> dat.map = data.frame(Index)
	> for (x in dat.map$Index) {
  		dat.map$c1[x] = 'M'
  		dat.map$c2[x] = mod.map$V4[x]
	}
	> dat.map$Index = NULL
	> head(dat.map)
	> write.table(dat.map, paste0(infile,'.dat'), sep = '\t', col.names = F, row.names = F, quote = F)	

The following step was carried out using MACH for pre-imputation quality control:
	$ mach1 -d /Volumes/MHS/Other/Shea_Imputation/analyses/MACH_IMPUTATION/easteal_ADNI_merged_biallelic.dat -p /Volumes/MHS/Other/Shea_Imputation/analyses/MACH_IMPUTATION/easteal_ADNI_merged_biallelic.ped --greedy --prefix /Volumes/MHS/Other/Shea_Imputation/analyses/MACH_IMPUTATION/easteal_ADNI_merged_biallelic_step1

At this step the following error was thrown, showing that we needed to remove multi allelic sites:
	FATAL ERROR:
	This version of MaCH is designed for bi-allelic markers
	However, the following marker(s) have >2 alleles:
   	12432 9301 6401 2260 16488 10034 327 6158 10843 2359 2779 11044 12702 14905 14580 5606 
   	12865 7852 1375 494 15261 12400 9293 16490 8471 491 931 8192 4943 5492 10149 13074 11710 
   	14395 417 3438 12019 7241 12621 9653 4736 14497 10739 3105 12375 6806 10742 10361 11047 
   	7257 14774 11929 14211 8469 3696 9299 5375 15061 10897 9287 3327 6239 4508 12574 16518 8496 
   	16103 188 15752 3411 12684 12026 11329 14550 15927 11229 13285 8701 4580 15861 9944 10303 
   	348 12437 5821 13884 15705 616 7705 13596 12964 3663 3843 3317 1914 6722 6881 11090 8206 
   	4005 9125 11437 13428 13351 4203 10181 6022 15316 7595 4733 4161 5267 9950 13368 13387 
   	12028 5436 3450 4158 10892 14894 411 8233 11242 6779 10197 10499 3531 13934 12438 10355 
   	3429 6927 591 14161 4222 7785 6367 15154 9899 14551 4031 2863 3651 5189 6134 3615 8805 12040 
   	11632 15479 13716 8631 9065 6649 8392 326 12720 16052 15663 5372 7795 16352 5126 5252 11698 
   	6773 287 10235 9055 4561 16070 8170 10895 12610 515 7581 7343 14049 9183 13183 5910 14776 
   
	2452 additional markers not listed

	Please remove or recode markers with more than 2 alleles

At this point, I attempted to use the IMPUTE2 programme:
