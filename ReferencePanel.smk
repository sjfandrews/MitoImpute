'''Snakefile for Construction of Reference Panel 0.1'''
# snakemake -s ReferencePanel.smk
# snakemake -s ReferencePanel.smk --dag | dot -Tsvg > dag_ReferencePanel.svg

configfile: "ReferencePanel_config.yaml"
BPLINK = ["bed", "bim", "fam"]
SAMPLE = config['sample']

DATAIN = 'data/ReferencePanel'
DATAOUT = 'DerivedData/test'
FILENAME = 'McInerney_Master_Alignment_July18_2018'

rule all:
    input:
        expand(DATAOUT + "/ReferencePanel.{ext}", ext = ['hap.gz', 'legend.gz']),
        expand(DATAOUT + "/ReferencePanel.{ext}", ext = ['ped', 'map']),
        expand(DATAOUT + "/ReferencePanel.{ext}", ext = ['gen.gz', 'samples']),
        expand(DATAOUT + "/MtMap.txt"),
        expand(DATAOUT + "/MtStrand.txt")

## 1. Run the ambiguous2missing.py script to change ambiguous character states to missing data:
rule ambiguous2missing:
    input:
        in_script = "scripts/PYTHON/ambiguous2missing.py",
        in_fasta = expand(DATAIN + "/{Reference}.fasta", Reference=FILENAME)
    output:
        out = temp(expand(DATAOUT + "/{Reference}_ambig2missing.fasta", Reference=FILENAME)),
    shell:
        'python {input.in_script} -i {input.in_fasta} -o {output.out} -v'

# 2. Identify samples with highQuality sequences
rule LowQualitySequences:
    input:
        in_script = "scripts/R/removeLowQuality_cmdline.R",
        in_fasta = expand(DATAOUT + "/{Reference}_ambig2missing.fasta", Reference=FILENAME),
    output:
        out = expand(DATAOUT + "/ReferencePanel_highQualitySequences.txt", Reference=FILENAME),
    shell:
        'Rscript {input.in_script} {input.in_fasta} {output.out}'

## 3a. Run the fasta2vcf_mtDNA.py script
rule fasta2vcf:
    input:
        in_script = "scripts/PYTHON/fasta2vcf_mtDNA.py",
        in_fasta = expand(DATAOUT + "/{Reference}_ambig2missing.fasta", Reference=FILENAME)
    output:
        out_vcf = temp(expand(DATAOUT + "/{Reference}_ambig2missing.vcf.gz", Reference=FILENAME))
    shell:
        'python {input.in_script} -i {input.in_fasta} -o {output.out_vcf} -v'

# 3b. Pass the resulting VCF through BCFTOOLS to make sure it conforms to all standards and index it
rule VcfCheck:
    input:
        in_vcf = expand(DATAOUT + "/{Reference}_ambig2missing.vcf.gz", Reference=FILENAME),
    output:
        out_vcf = temp(DATAOUT + "/Reference_panal.vcf.gz"),
    run:
        shell('bcftools view -Oz -o {output.out_vcf} {input.in_vcf}')
        shell('bcftools index {output.out_vcf}')

## 4. Remove low quality sequences from VCF
rule RemoveLowQuality:
    input:
        quality = DATAOUT + "/ReferencePanel_highQualitySequences.txt",
        in_vcf = DATAOUT + "/Reference_panal.vcf.gz",
    output:
        out_vcf = temp(DATAOUT + "/ReferencePanel_highQual.vcf.gz"),
    run:
        shell('bcftools view --force-samples -S {input.quality} -Oz -o {output.out_vcf} {input.in_vcf}')
        shell('bcftools index {output.out_vcf}')

## 5. Apply filtration criteria
rule SiteFiltration:
    input:
        in_vcf = DATAOUT + "/ReferencePanel_highQual.vcf.gz",
    output:
        out_vcf = DATAOUT + "/ReferencePanel.vcf.gz",
    run:
        shell('vt decompose {input.in_vcf} | bcftools +fill-tags | bcftools view -i \'ALT!="-" \' | bcftools view -q 0.01 -Q 0.99 | bcftools view -Oz -o {output.out_vcf}')
        shell('bcftools index {output.out_vcf}')

## 6a. Extract sample names from Reference Panel
rule RefSampleNames:
    input:
        in_vcf = DATAOUT + "/ReferencePanel.vcf.gz",
    output:
        out_samples = DATAOUT + "/RefSampleList.txt",
    shell:
        'bcftools query -l {input.in_vcf} > {output.out_samples}'

## 6b. Assign M sex label to reference Samples
rule RefSampleSex:
    input:
        in_samples = DATAOUT + "/RefSampleList.txt",
        in_script = "scripts/R/assign_sex_label.R"
    output:
        outfile = DATAOUT + "/RefSampleList_sex.txt",
    shell:
        'Rscript {input.in_script} {input.in_samples} {output.outfile}'

## 7. Convert to Oxford format
rule Oxford:
    input:
        in_vcf = DATAOUT + "/ReferencePanel.vcf.gz",
        in_sex = DATAOUT + "/RefSampleList_sex.txt"
    output:
        expand(DATAOUT + "/ReferencePanel.{ext}", ext = ['hap.gz', 'legend.gz'])
    params:
        out = DATAOUT + "/ReferencePanel"
    shell:
        'bcftools convert --haplegendsample {params.out} {input.in_vcf} --sex {input.in_sex}'

## 8. Generate .ped and .map files
rule Plink:
    input:
        in_vcf = DATAOUT + "/ReferencePanel.vcf.gz",
    output:
        expand(DATAOUT + "/ReferencePanel.{ext}", ext = ['ped', 'map'])
    params:
        out = DATAOUT + "/ReferencePanel"
    shell:
        'plink --vcf {input.in_vcf} --recode --double-id --keep-allele-order --out {params.out}'

## 9. Generate .gen and .sample files
rule GenSample:
    input:
        in_vcf = DATAOUT + "/ReferencePanel.vcf.gz",
        in_sex = DATAOUT + "/RefSampleList_sex.txt"
    output:
        expand(DATAOUT + "/ReferencePanel.{ext}", ext = ['gen.gz', 'samples'])
    params:
        out = DATAOUT + "/ReferencePanel"
    shell:
        'bcftools convert --gensample {params.out} {input.in_vcf} --sex {input.in_sex}'

## 10. Construct .map file for IMPUTE2
rule MakeMapFile:
    input:
        in_script = "scripts/R/mt_recombination_map.R",
        in_vcf = DATAOUT + "/ReferencePanel.vcf.gz"
    output:
        out_map = DATAOUT + "/MtMap.txt",
        out_strand = DATAOUT + "/MtStrand.txt"
    shell:
        'Rscript {input.in_script} {input.in_vcf} {output.out_map} {output.out_strand}'
