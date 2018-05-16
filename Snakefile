'''Snakefile for MitoImpute Version 0.1'''

configfile: "config.yaml"
BPLINK = ["bed", "bim", "fam"]
SAMPLE = config['sample']

rule all:
    input:
        "DerivedData/McInerney_Master_Alignment_Nov30_2017_ambig2missing.vcf.gz"

## 1. Run the ambiguous2missing.py script to change ambiguous character states to missing data:
rule ambiguous2missing:
    input:
        "scripts/PYTHON/ambiguous2missing.py",
        "data/McInerney_Master_Alignment_Nov30_2017.fasta",
    output:
        "DerivedData/McInerney_Master_Alignment_Nov30_2017_ambig2missing.fasta",
    params:
        in_fasta = "data/McInerney_Master_Alignment_Nov30_2017.fasta",
        in_script = "scripts/PYTHON/ambiguous2missing.py",
        out = "DerivedData/McInerney_Master_Alignment_Nov30_2017_ambig2missing.fasta"
    shell:
        'python {params.in_script} -i {params.in_fasta} -o {params.out} -v'

## 2. Run the fasta2vcf_mtDNA.py script
rule fasta2vcf:
    input:
        "scripts/PYTHON/fasta2vcf_mtDNA.py",
        "DerivedData/McInerney_Master_Alignment_Nov30_2017_ambig2missing.fasta"
    output:
        "DerivedData/McInerney_Master_Alignment_Nov30_2017_ambig2missing.vcf.gz"
    params:
        in_fasta = "DerivedData/McInerney_Master_Alignment_Nov30_2017_ambig2missing.fasta",
        in_script = "scripts/PYTHON/fasta2vcf_mtDNA.py",
        out = "DerivedData/McInerney_Master_Alignment_Nov30_2017_ambig2missing.vcf.gz"
    shell:
        'python {params.in_script} -i {params.in_fasta} -o {params.out} -v'

## 3. Pass the resulting VCF through BCFTOOLS to make sure it conforms to all standards and index it
#rule bcf:
#    input:
#        "~/Documents/MitoImpute/DerivedData/McInerney_Master_Alignment_Nov30_2017_ambig2missing.vcf.gz",
#    output:
#        "~/Documents/MitoImpute/DerivedData/McInerney_Master_Alignment_Nov30_2017_ambig2missing.vcf.gz",
#    params:
#        in_vcf = "~/Documents/MitoImpute/DerivedData/McInerney_Master_Alignment_Nov30_2017_ambig2missing.vcf.gz",
#        out = "~/Documents/MitoImpute/DerivedData/McInerney_Master_Alignment_Nov30_2017_ambig2missing.vcf.gz",
#    shell:
#        'bcftools view -Oz -o {params.out} {params.vcf} | bcftools index'

## 4a. Remove low quality samples
#$ bcftools view -S ^low_quality_samples.txt -Oz -o McInerney_Master_Alignment_Nov30_2017_ambig2missing_highQual.vcf.gz McInerney_Master_Alignment_Nov30_2017_ambig2missing.vcf.gz
#$ bcftools index McInerney_Master_Alignment_Nov30_2017_ambig2missing_highQual.vcf.gz

## 4b. Only include high quality samples
#$ bcftools view -S high_quality_samples.txt -Oz -o McInerney_Master_Alignment_Nov30_2017_ambig2missing_highQual.vcf.gz McInerney_Master_Alignment_Nov30_2017_ambig2missing.vcf.gz
#$ bcftools index McInerney_Master_Alignment_Nov30_2017_ambig2missing_highQual.vcf.gz

## 5. Apply filtration criteria
#$ bcftools +fill-tags McInerney_Master_Alignment_Nov30_2017_ambig2missing_highQual.vcf.gz  | bcftools norm -m -any | bcftools view -q 0.01 -Q 0.99 | bcftools view -i 'ALT!="*" && POS!=302 && POS!=303 && POS!=308 && POS!=309 && POS!=310 && POS!=513 && POS!=515 && POS!=522 && POS!=523 && POS!=3106 && POS!=3107' -Oz -o McInerney_Master_Alignment_Nov30_2017_ambig2missing_highQual_filtered.vcf.gz
#$ bcftools index McInerney_Master_Alignment_Nov30_2017_ambig2missing_highQual_filtered.vcf.gz

## 6. Extract sample names and assign M sex label
#$ bcftools query -l McInerney_Master_Alignment_Nov30_2017_ambig2missing_highQual_filtered.vcf.gz > McInerney_Master_Alignment_Nov30_2017_ambig2missing_highQual_filtered_SEX.txt
#$ Rscript assign_sex_labels.R McInerney_Master_Alignment_Nov30_2017_ambig2missing_highQual_filtered_SEX.txt

#$ bcftools query -l ADNI_samples.vcf.gz > ADNI_samples_SEX.txt
#$ Rscript assign_sex_labels.R ADNI_samples_SEX.txt

## 7. Convert to Oxford format
#$ bcftools convert --haplegendsample McInerney_Master_Alignment_Nov30_2017_ambig2missing_highQual_filtered McInerney_Master_Alignment_Nov30_2017_ambig2missing_highQual_filtered.vcf.gz --sex McInerney_Master_Alignment_Nov30_2017_ambig2missing_highQual_filtered_SEX.txt

#$ bcftools convert --haplegendsample ADNI_samples ADNI_samples.vcf.gz --sex ADNI_samples_SEX.txt

## 8. Generate .ped and .map files
#$ plink --vcf McInerney_Master_Alignment_Nov30_2017_ambig2missing_highQual_filtered.vcf.gz --recode --double-id --out McInerney_Master_Alignment_Nov30_2017_ambig2missing_highQual_filtered

#$ plink --vcf ADNI_samples.vcf.gz --recode --double-id --out ADNI_samples

## 9. Generate .gen and .sample files
#$ bcftools convert --gensample McInerney_Master_Alignment_Nov30_2017_ambig2missing_highQual_filtered McInerney_Master_Alignment_Nov30_2017_ambig2missing_highQual_filtered.vcf.gz --sex McInerney_Master_Alignment_Nov30_2017_ambig2missing_highQual_filtered_SEX.txt

#$ bcftools convert --gensample ADNI_samples ADNI_samples.vcf.gz --sex ADNI_samples_SEX.txt

## 10. Run IMPUTE2
#$ impute2 -chrX -m McInerney_Master_Alignment_Nov30_2017_ambig2missing_highQual_filtered map -h McInerney_Master_Alignment_Nov30_2017_ambig2missing_highQual_filtered hap.gz -l McInerney_Master_Alignment_Nov30_2017_ambig2missing_highQual_filtered legend.gz -g ADNI_samples gen.gz -sample_g ADNI_samples samples -int 1 16569 -Ne 20000 -o ADNI_samples_IMPUTED

## 11. Run R scripts to determine haplogroup concordance
#$ <FILL THIS OUT ONCE COMMAND-LINE VERSIONS ARE WRITTEN>
