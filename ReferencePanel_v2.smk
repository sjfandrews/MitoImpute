'''Snakefile for Construction of Reference Panel 0.1'''
# snakemake -s ReferencePanel_v2.smk
# snakemake -s ReferencePanel_v2.smk --dag | dot -Tsvg > dag_ReferencePanel_v2.svg

configfile: "ReferencePanel_config.yaml"
BPLINK = ["bed", "bim", "fam"]
SAMPLE = config['sample']

rule all:
    input:
        expand("/g/data1a/te53/MitoImpute/data/REF_PANEL/Reference_panel_v2.{ext}", ext = ['hap.gz', 'legend.gz']),
        expand("/g/data1a/te53/MitoImpute/data/REF_PANEL/Reference_panel_v2.{ext}", ext = ['ped', 'map']),
        expand("/g/data1a/te53/MitoImpute/data/REF_PANEL/Reference_panel_v2.{ext}", ext = ['gen.gz', 'samples']),
        "DerivedData/ReferencePanel/MtMap.txt",
        "DerivedData/ReferencePanel/MtStrand.txt"

## 1. Run the ambiguous2missing.py script to change ambiguous character states to missing data:
rule ambiguous2missing:
    input:
        "scripts/PYTHON/ambiguous2missing.py",
        "/g/data1a/te53/MitoImpute/data/FASTA/masters/McInerney_Master_Alignment_July18_2018.fasta",
    output:
        temp("/g/data1a/te53/MitoImpute/data/FASTA/ambiguous2missing/McInerney_Master_Alignment_July18_2018_ambig2missing.fasta"),
    params:
        in_fasta = "/g/data1a/te53/MitoImpute/data/FASTA/masters/McInerney_Master_Alignment_July18_2018.fasta",
        in_script = "scripts/PYTHON/ambiguous2missing.py",
        out = "/g/data1a/te53/MitoImpute/data/FASTA/ambiguous2missing/McInerney_Master_Alignment_July18_2018_ambig2missing.fasta"
    shell:
        'python {params.in_script} -i {params.in_fasta} -o {params.out} -v'
        
## 2. Run the fasta2vcf_mtDNA.py script
rule fasta2vcf:
    input:
        "scripts/PYTHON/fasta2vcf_mtDNA.py",
        "/g/data1a/te53/MitoImpute/data/FASTA/ambiguous2missing/McInerney_Master_Alignment_July18_2018_ambig2missing.fasta"
    output:
        "/g/data1a/te53/MitoImpute/data/VCF/McInerney_Master_Alignment_July18_2018_ambig2missing.vcf.gz"
    params:
        in_fasta = "/g/data1a/te53/MitoImpute/data/FASTA/ambiguous2missing/McInerney_Master_Alignment_July18_2018_ambig2missing.fasta",
        in_script = "scripts/PYTHON/fasta2vcf_mtDNA.py",
        out = "/g/data1a/te53/MitoImpute/data/VCF/McInerney_Master_Alignment_July18_2018_ambig2missing.vcf.gz"
    shell:
        'python {params.in_script} -i {params.in_fasta} -o {params.out} -v'

# 3. Pass the resulting VCF through BCFTOOLS to make sure it conforms to all standards and index it
rule VcfCheck:
    input:
        "/g/data1a/te53/MitoImpute/data/VCF/McInerney_Master_Alignment_July18_2018_ambig2missing.vcf.gz",
    output:
        temp("/g/data1a/te53/MitoImpute/data/REF_PANEL/Reference_panel_v2.vcf.gz"),
    params:
        in_vcf = "/g/data1a/te53/MitoImpute/data/VCF/McInerney_Master_Alignment_July18_2018_ambig2missing.vcf.gz",
        out_vcf = "/g/data1a/te53/MitoImpute/data/REF_PANEL/Reference_panel_v2.vcf.gz",
    run:
        shell('bcftools view -Oz -o {params.out_vcf} {params.in_vcf}')
        shell('bcftools index {params.out_vcf}')

# 4a. Identify samples with highQuality sequences
rule LowQualitySequences:
    input:
        "/g/data1a/te53/MitoImpute/data/FASTA/ambiguous2missing/McInerney_Master_Alignment_July18_2018_ambig2missing.fasta",
        "scripts/R/removeLowQuality_cmdline.R",
    output:
        temp("/g/data1a/te53/MitoImpute/metadata/Reference_panel_v2_highQualitySequences.txt"),
    params:
        in_fasta = "/g/data1a/te53/MitoImpute/data/FASTA/ambiguous2missing/McInerney_Master_Alignment_July18_2018_ambig2missing.fasta",
        in_script = "scripts/R/removeLowQuality_cmdline.R",
        out = "/g/data1a/te53/MitoImpute/metadata/Reference_panel_v2_highQualitySequences.txt",
    shell:
        'Rscript {params.in_script} {params.in_fasta} {params.out}'

## 4b. Remove low quality sequences
rule RemoveLowQuality:
    input:
        "/g/data1a/te53/MitoImpute/metadata/Reference_panel_v2_highQualitySequences.txt",
        "/g/data1a/te53/MitoImpute/data/REF_PANEL/Reference_panel_v2.vcf.gz",
    output:
        temp("/g/data1a/te53/MitoImpute/data/REF_PANEL/Reference_panel_v2_highQual.vcf.gz"),
    params:
        in_vcf = "/g/data1a/te53/MitoImpute/data/REF_PANEL/Reference_panel_v2.vcf.gz",
        quality = "/g/data1a/te53/MitoImpute/metadata/Reference_panel_v2_highQualitySequences.txt",
        out_vcf = "/g/data1a/te53/MitoImpute/data/REF_PANEL/Reference_panel_v2_highQual.vcf.gz",
    run:
        shell('bcftools view --force-samples -S {params.quality} -Oz -o {params.out_vcf} {params.in_vcf}')
        shell('bcftools index {params.out_vcf}')


## 5. Apply filtration criteria
rule SiteFiltration:
    input:
        "/g/data1a/te53/MitoImpute/data/REF_PANEL/Reference_panel_v2_highQual.vcf.gz",
    output:
        "/g/data1a/te53/MitoImpute/data/REF_PANEL/Reference_panel_v2_highQual_filtered.vcf.gz",
    params:
        in_vcf = "/g/data1a/te53/MitoImpute/data/REF_PANEL/Reference_panel_v2_highQual.vcf.gz",
        out_vcf = "/g/data1a/te53/MitoImpute/data/REF_PANEL/Reference_panel_v2_highQual_filtered.vcf.gz",
    run:
        shell('vt decompose {params.in_vcf} | bcftools +fill-tags | bcftools view -i \'ALT!="-" \' | bcftools view -q 0.01 -Q 0.99 | bcftools view -Oz -o {params.out_vcf}')
        shell('bcftools index {params.out_vcf}')

## 6a. Extract sample names from Reference Panel
rule RefSampleNames:
    input:
        "/g/data1a/te53/MitoImpute/data/REF_PANEL/Reference_panel_v2_highQual_filtered.vcf.gz",
    output:
        "/g/data1a/te53/MitoImpute/metadata/Reference_panel_v2_SampleList.txt",
    params:
        in_vcf = "/g/data1a/te53/MitoImpute/data/REF_PANEL/Reference_panel_v2_highQual_filtered.vcf.gz",
        out_samples = "/g/data1a/te53/MitoImpute/metadata/Reference_panel_v2_SampleList.txt",
    shell:
        'bcftools query -l {params.in_vcf} > {params.out_samples}'

## 6b. Assign M sex label to reference Samples
rule RefSampleSex:
    input:
        "/g/data1a/te53/MitoImpute/metadata/Reference_panel_v2_SampleList.txt",
        "scripts/R/assign_sex_label.R"
    output:
        "/g/data1a/te53/MitoImpute/metadata/Reference_panel_v2_SampleList_sex.txt",
    params:
        in_script = "scripts/R/assign_sex_label.R",
        in_samples = "/g/data1a/te53/MitoImpute/metadata/Reference_panel_v2_SampleList.txt",
        outfile = "/g/data1a/te53/MitoImpute/metadata/Reference_panel_v2_SampleList_sex.txt"
    shell:
        'Rscript {params.in_script} {params.in_samples} {params.outfile}'

## 7. Convert to Oxford format
rule Oxford:
    input:
        "/g/data1a/te53/MitoImpute/data/REF_PANEL/Reference_panel_v2_highQual_filtered.vcf.gz",
        "/g/data1a/te53/MitoImpute/metadata/Reference_panel_v2_SampleList_sex.txt"
    output:
        expand("/g/data1a/te53/MitoImpute/data/REF_PANEL/Reference_panel_v2.{ext}", ext = ['hap.gz', 'legend.gz'])
    params:
        in_vcf = "/g/data1a/te53/MitoImpute/data/REF_PANEL/Reference_panel_v2_highQual_filtered.vcf.gz",
        in_sex = "/g/data1a/te53/MitoImpute/metadata/Reference_panel_v2_SampleList_sex.txt",
        out = "/g/data1a/te53/MitoImpute/data/REF_PANEL/Reference_panel_v2"
    shell:
        'bcftools convert --haplegendsample {params.out} {params.in_vcf} --sex {params.in_sex}'

## 8. Generate .ped and .map files
rule Plink:
    input:
        "/g/data1a/te53/MitoImpute/data/REF_PANEL/Reference_panel_v2_highQual_filtered.vcf.gz",
    output:
        expand("/g/data1a/te53/MitoImpute/data/REF_PANEL/Reference_panel_v2.{ext}", ext = ['ped', 'map'])
    params:
        in_vcf = "/g/data1a/te53/MitoImpute/data/REF_PANEL/Reference_panel_v2_highQual_filtered.vcf.gz",
        out = "/g/data1a/te53/MitoImpute/data/REF_PANEL/Reference_panel_v2"
    shell:
        'plink --vcf {params.in_vcf} --recode --double-id --out {params.out}'

## 9. Generate .gen and .sample files
rule GenSample:
    input:
        "/g/data1a/te53/MitoImpute/data/REF_PANEL/Reference_panel_v2_highQual_filtered.vcf.gz",
        "/g/data1a/te53/MitoImpute/metadata/Reference_panel_v2_SampleList_sex.txt"
    output:
        expand("/g/data1a/te53/MitoImpute/data/REF_PANEL/Reference_panel_v2.{ext}", ext = ['gen.gz', 'samples'])
    params:
        in_vcf = "/g/data1a/te53/MitoImpute/data/REF_PANEL/Reference_panel_v2_highQual_filtered.vcf.gz",
        in_sex = "g/data1a/te53/MitoImpute/metadata/Reference_panel_v2_SampleList_sex.txt",
        out = "/g/data1a/te53/MitoImpute/data/REF_PANEL/Reference_panel_v2"
    shell:
        'bcftools convert --gensample {params.out} {params.in_vcf} --sex {params.in_sex}'

## 10. Construct .map file for IMPUTE2
rule MakeMapFile:
    input:
        "scripts/R/mt_recombination_map.R",
        "/g/data1a/te53/MitoImpute/data/REF_PANEL/Reference_panel_v2_highQual_filtered.vcf.gz"
    output:
        "DerivedData/ReferencePanel/MtMap.txt",
        "DerivedData/ReferencePanel/MtStrand.txt"
    params:
        in_vcf = "/g/data1a/te53/MitoImpute/data/REF_PANEL/Reference_panel_v2_highQual_filtered.vcf.gz",
        in_script = "scripts/R/mt_recombination_map.R",
        out_map = "DerivedData/ReferencePanel/MtMap.txt",
        out_strand = "DerivedData/ReferencePanel/MtStrand.txt"
    shell:
        'Rscript {params.in_script} {params.in_vcf} {params.out_map} {params.out_strand}'


#$ bcftools convert --haplegendsample ADNI_samples ADNI_samples.vcf.gz --sex ADNI_samples_SEX.txt
#$ plink --vcf ADNI_samples.vcf.gz --recode --double-id --out ADNI_samples
#$ bcftools convert --gensample ADNI_samples ADNI_samples.vcf.gz --sex ADNI_samples_SEX.txt

## 10. Run IMPUTE2
#$ impute2 -chrX -m McInerney_Master_Alignment_Nov30_2017_ambig2missing_highQual_filtered map -h McInerney_Master_Alignment_Nov30_2017_ambig2missing_highQual_filtered hap.gz -l McInerney_Master_Alignment_Nov30_2017_ambig2missing_highQual_filtered legend.gz -g ADNI_samples gen.gz -sample_g ADNI_samples samples -int 1 16569 -Ne 20000 -o ADNI_samples_IMPUTED

## 11. Run R scripts to determine haplogroup concordance
#$ <FILL THIS OUT ONCE COMMAND-LINE VERSIONS ARE WRITTEN>