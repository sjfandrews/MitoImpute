To change ambiguous character states to missing character state:
$ python ~/GitCode/timmci/scripts/Other_Scripts/ambiguous2missing.py -i /Volumes/MHS/Master_Alignment/FASTA/McInerney_Master_Alignment_Nov30_2017.fasta -v

To convert the resulting fasta file to a VCF file with the reference alleles pegged to the rCRS:
$ python ~/GitCode/timmci/scripts/Other_Scripts/fasta2vcf_mtDNA.py -i /Volumes/MHS/Master_Alignment/FASTA/McInerney_Master_Alignment_Nov30_2017_ambig2missing.fasta -o /Volumes/MHS/Master_Alignment/fasta2vcf_generated/McInerney_Master_Alignment_Nov30_2017_ambig2missing.vcf -v

To compress the VCF file:
$ bcftools view -Oz -o /Volumes/MHS/Master_Alignment/fasta2vcf_generated/McInerney_Master_Alignment_Nov30_2017_ambig2missing.vcf.gz /Volumes/MHS/Master_Alignment/fasta2vcf_generated/McInerney_Master_Alignment_Nov30_2017_ambig2missing.vcf
$ bcftools index /Volumes/MHS/Master_Alignment/fasta2vcf_generated/McInerney_Master_Alignment_Nov30_2017_ambig2missing.vcf.gz

To filter to exclude all ancient samples, all non H.sapiens samples, all non-European samples, and all samples with missing > 5 or gaps > 7
$ bcftools view -S ^/Volumes/MHS/Other/Shea_Imputation/metadata/exclusion_lists/BIG_EXCLUSION_euro_dupsRemoved.txt -Oz -o /Volumes/MHS/Other/Shea_Imputation/data_files/VCF_files/EuroOnly/McInerneyMaster_EuroOnly.vcf.gz /Volumes/MHS/Master_Alignment/fasta2vcf_generated/McInerney_Master_Alignment_Nov30_2017_ambig2missing.vcf.gz
$ bcftools index /Volumes/MHS/Other/Shea_Imputation/data_files/VCF_files/EuroOnly/McInerneyMaster_EuroOnly.vcf.gz

To apply filtration criteria:
$ bcftools +fill-tags /Volumes/TimMcInerney/Other/Shea_Imputation/data_files/VCF_files/EuroOnly/McInerneyMaster_EuroOnly.vcf.gz | bcftools norm -m -any | bcftools view -q 0.01 -Q 0.99 | bcftools view -i 'ALT!="*" && POS!=302 && POS!=303 && POS!=308 && POS!=309 && POS!=310 && POS!=513 && POS!=515 && POS!=522 && POS!=523 && POS!=3106 && POS!=3107' -Oz -o /Volumes/TimMcInerney/Other/Shea_Imputation/data_files/VCF_files/EuroOnly/MOD/McInerneyMaster_EuroOnly_SHEAspecs.vcf.gz
$ bcftools index /Volumes/TimMcInerney/Other/Shea_Imputation/data_files/VCF_files/EuroOnly/MOD/McInerneyMaster_EuroOnly_SHEAspecs.vcf.gz

To extract sample names and assign M sex:
$ bcftools query -l /Volumes/TimMcInerney/Other/Shea_Imputation/data_files/VCF_files/EuroOnly/MOD/McInerneyMaster_EuroOnly_SHEAspecs.vcf.gz > /Volumes/TimMcInerney/Other/Shea_Imputation/data_files/VCF_files/EuroOnly/MOD/McInerneyMaster_EuroOnly_SHEAspecs_SEX.txt
In R:
> x = read.table('/Volumes/TimMcInerney/Other/Shea_Imputation/data_files/VCF_files/EuroNonEuro/MOD/McInerneyMaster_EuroNonEuro_SHEAspecs_SEX.txt', header = F, sep = '\t')
> x$V2 = 'M'
> write.table(x, '/Volumes/TimMcInerney/Other/Shea_Imputation/data_files/VCF_files/EuroNonEuro/MOD/McInerneyMaster_EuroNonEuro_SHEAspecs_SEX.txt', col.names = F, row.names = F, sep = '\t', quote = F)

To convert to OXFORD format:
$ bcftools convert --haplegendsample McInerneyMaster_EuroOnly_SHEAspecs /Volumes/TimMcInerney/Other/Shea_Imputation/data_files/VCF_files/EuroOnly/MOD/McInerneyMaster_EuroOnly_SHEAspecs.vcf.gz --sex /Volumes/TimMcInerney/Other/Shea_Imputation/data_files/VCF_files/EuroOnly/MOD/McInerneyMaster_EuroOnly_SHEAspecs_SEX.txt

IMPUTE2 commands:
$ impute2 -chrX -m /Volumes/TimMcInerney/Other/Shea_Imputation/data_files/plink_files/McI_ref/McInerneyMaster_EuroNonEuro_SHEAspecs.map -h /Volumes/TimMcInerney/Other/Shea_Imputation/data_files/OXFORD/EuroNonEuro/McInerneyMaster_EuroNonEuro_SHEAspecs.hap.gz -l /Volumes/TimMcInerney/Other/Shea_Imputation/data_files/OXFORD/EuroNonEuro/McInerneyMaster_EuroNonEuro_SHEAspecs.legend.gz -g /Volumes/TimMcInerney/Other/Shea_Imputation/data_files/OXFORD/sample_panels/SamplePanel_180221.gen.gz -sample_g /Volumes/TimMcInerney/Other/Shea_Imputation/data_files/OXFORD/sample_panels/SamplePanel_180221.samples -int 1 16569 -Ne 20000 -o /Volumes/TimMcInerney/Other/Shea_Imputation/data_files/IMPUTE2/mito_snps_rcrs_ed_nameMod