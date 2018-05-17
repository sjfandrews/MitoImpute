'''Snakefile for MitoImpute Version 0.1'''

configfile: "config.yaml"
BPLINK = ["bed", "bim", "fam"]
SAMPLE = config['sample']

rule all:
    input:
        expand("DerivedData/RefencePanel.{ext}", ext = ['hap.gz', 'legend.gz']),
        expand("DerivedData/RefencePanel.{ext}", ext = ['ped', 'map']),
        expand("DerivedData/RefencePanel.{ext}", ext = ['gen.gz', 'samples']),
        "DerivedData/MtMap.txt"

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

# 3. Pass the resulting VCF through BCFTOOLS to make sure it conforms to all standards and index it
rule VcfCheck:
    input:
        "DerivedData/McInerney_Master_Alignment_Nov30_2017_ambig2missing.vcf.gz",
    output:
        "DerivedData/Refence_panal.vcf.gz",
    params:
        in_vcf = "DerivedData/McInerney_Master_Alignment_Nov30_2017_ambig2missing.vcf.gz",
        out_vcf = "DerivedData/Refence_panal.vcf.gz",
    run:
        shell('bcftools view -Oz -o {params.out_vcf} {params.in_vcf}')
        shell('bcftools index {params.out_vcf}')

# 4a. Identify samples with highQuality sequences
rule LowQualitySequences:
    input:
        "DerivedData/McInerney_Master_Alignment_Nov30_2017_ambig2missing.fasta",
        "scripts/R/removeLowQuality_cmdline.R",
    output:
        "scripts/INFORMATION_LISTS/ReferencePanel_highQualitySequences.txt",
    params:
        in_fasta = "DerivedData/McInerney_Master_Alignment_Nov30_2017_ambig2missing.fasta",
        in_script = "scripts/R/removeLowQuality_cmdline.R",
        out = "scripts/INFORMATION_LISTS/ReferencePanel_highQualitySequences.txt",
    shell:
        'Rscript {params.in_script} {params.in_fasta} {params.out}'

## 4b. Remove low quality samples
rule RemoveLowQuality:
    input:
        "scripts/INFORMATION_LISTS/ReferencePanel_highQualitySequences.txt",
        "DerivedData/Refence_panal.vcf.gz",
    output:
        "DerivedData/RefencePanel_highQual.vcf.gz",
    params:
        in_vcf = "DerivedData/Refence_panal.vcf.gz",
        quality = "scripts/INFORMATION_LISTS/ReferencePanel_highQualitySequences.txt",
        out_vcf = "DerivedData/RefencePanel_highQual.vcf.gz",
    run:
        shell('bcftools view -S ^{params.quality} -Oz -o {params.out_vcf} {params.in_vcf}')
        shell('bcftools index {params.out_vcf}')

## 5. Apply filtration criteria
rule SiteFiltration:
    input:
        "DerivedData/RefencePanel_highQual.vcf.gz",
    output:
        "DerivedData/RefencePanel_highQual_filtered.vcf.gz",
    params:
        in_vcf = "DerivedData/RefencePanel_highQual.vcf.gz",
        out_vcf = "DerivedData/RefencePanel_highQual_filtered.vcf.gz",
    run:
        shell('bcftools +fill-tags {params.in_vcf} | bcftools norm -m -any | bcftools view -q 0.01 -Q 0.99 | bcftools view -i \'ALT!="*" && POS!=302 && POS!=303 && POS!=308 && POS!=309 && POS!=310 && POS!=513 && POS!=515 && POS!=522 && POS!=523 && POS!=3106 && POS!=3107\' -Oz -o {params.out_vcf}')
        shell('bcftools index {params.out_vcf}')

## 6a. Extract sample names from Reference Panel
rule RefSampleNames:
    input:
        "DerivedData/RefencePanel_highQual_filtered.vcf.gz",
    output:
        "scripts/INFORMATION_LISTS/RefSampleList.txt",
    params:
        in_vcf = "DerivedData/RefencePanel_highQual.vcf.gz",
        out_samples = "scripts/INFORMATION_LISTS/RefSampleList.txt",
    shell:
        'bcftools query -l {params.in_vcf} > {params.out_samples}'

## 6b. Assign M sex label to reference Samples
rule RefSampleSex:
    input:
        "scripts/INFORMATION_LISTS/RefSampleList.txt",
        "scripts/R/assign_sex_label.R"
    output:
        "DerivedData/RefSampleList_sex.txt",
    params:
        in_script = "scripts/R/assign_sex_label.R",
        in_samples = "scripts/INFORMATION_LISTS/RefSampleList.txt",
        outfile = "DerivedData/RefSampleList_sex.txt"
    shell:
        'Rscript {params.in_script} {params.in_samples} {params.outfile}'

## 7. Convert to Oxford format
rule Oxford:
    input:
        "DerivedData/RefencePanel_highQual_filtered.vcf.gz",
        "DerivedData/RefSampleList_sex.txt"
    output:
        expand("DerivedData/RefencePanel.{ext}", ext = ['hap.gz', 'legend.gz'])
    params:
        in_vcf = "DerivedData/RefencePanel_highQual_filtered.vcf.gz",
        in_sex = "DerivedData/RefSampleList_sex.txt",
        out = "DerivedData/RefencePanel"
    shell:
        'bcftools convert --haplegendsample {params.out} {params.in_vcf} --sex {params.in_sex}'

## 8. Generate .ped and .map files
rule Plink:
    input:
        "DerivedData/RefencePanel_highQual_filtered.vcf.gz",
    output:
        expand("DerivedData/RefencePanel.{ext}", ext = ['ped', 'map'])
    params:
        in_vcf = "DerivedData/RefencePanel_highQual_filtered.vcf.gz",
        out = "DerivedData/RefencePanel"
    shell:
        'plink --vcf {params.in_vcf} --recode --double-id --out {params.out}'

## 9. Generate .gen and .sample files
rule GenSample:
    input:
        "DerivedData/RefencePanel_highQual_filtered.vcf.gz",
        "DerivedData/RefSampleList_sex.txt"
    output:
        expand("DerivedData/RefencePanel.{ext}", ext = ['gen.gz', 'samples'])
    params:
        in_vcf = "DerivedData/RefencePanel_highQual_filtered.vcf.gz",
        in_sex = "DerivedData/RefSampleList_sex.txt",
        out = "DerivedData/RefencePanel"
    shell:
        'bcftools convert --gensample {params.out} {params.in_vcf} --sex {params.in_sex}'

## 10. Construct .map file for IMPUTE2
rule MakeMapFile:
    input:
        "scripts/R/mt_recombination_map.R",
        "DerivedData/RefencePanel_highQual_filtered.vcf.gz"
    output:
        "DerivedData/MtMap.txt"
    params:
        in_vcf = "DerivedData/RefencePanel_highQual_filtered.vcf.gz",
        in_script = "scripts/R/mt_recombination_map.R",
        out = "DerivedData/MtMap.txt"
    shell:
        'Rscript {params.in_script} {params.in_vcf} {params.out}'


#$ bcftools convert --haplegendsample ADNI_samples ADNI_samples.vcf.gz --sex ADNI_samples_SEX.txt
#$ plink --vcf ADNI_samples.vcf.gz --recode --double-id --out ADNI_samples
#$ bcftools convert --gensample ADNI_samples ADNI_samples.vcf.gz --sex ADNI_samples_SEX.txt

## 10. Run IMPUTE2
#$ impute2 -chrX -m McInerney_Master_Alignment_Nov30_2017_ambig2missing_highQual_filtered map -h McInerney_Master_Alignment_Nov30_2017_ambig2missing_highQual_filtered hap.gz -l McInerney_Master_Alignment_Nov30_2017_ambig2missing_highQual_filtered legend.gz -g ADNI_samples gen.gz -sample_g ADNI_samples samples -int 1 16569 -Ne 20000 -o ADNI_samples_IMPUTED

## 11. Run R scripts to determine haplogroup concordance
#$ <FILL THIS OUT ONCE COMMAND-LINE VERSIONS ARE WRITTEN>
