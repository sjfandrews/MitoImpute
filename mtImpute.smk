'''Snakefile for MitoImpute Version 0.1'''
# snakemake -s mtImpute.smk
# snakemake -s mtImpute.smk --dag | dot -Tsvg > dag_mtImpute.svg

import os

configfile: 'mtImpute_config.yaml'
SAMPLE = config['SAMPLE']
DATAIN = config['DATAIN']
DATAOUT = config['DATAOUT']
REFDATA = config['REFDATA']
INFOCUT = config['INFOCUT']

BPLINK = ["bed", "bim", "fam"]
PLINK = ["map", "ped"]
OXFORD = ["gen", "sample"]
RWD = os.getcwd()

## For running on cluster
#mkdir .snakejob; snakejob -s mtImpute.smk -j 100 --until oxford2vcf oxford2ped oxford2bed bplink2plink --max-jobs-per-second 1 --keep-going
#shell.prefix('module load plink/1.90 impute2 R/3.4.3; ')

rule all:
    input:
        expand(DATAOUT + "/{sample}/Imputed_{sample}.{ext}", ext=BPLINK, sample=SAMPLE),
        expand(DATAOUT + "/{sample}/Imputed_{sample}.{ext}", ext=PLINK, sample=SAMPLE),
        expand(DATAOUT + "/{sample}/Imputed_{sample}.vcf", sample=SAMPLE),
        expand(DATAOUT + "/{sample}/stats/{sample}_mtImputed_QC.html", sample=SAMPLE),

rule SexFam:
    input: DATAIN + "/{sample}.fam"
    output: DATAOUT + "/{sample}/{sample}_maleOnly.fam"
    shell: """
        awk '{{$5 = "1"; print}}' {input} > {output}
    """
## Extract SNPs on mitochondrial genome (26)
## Set missing genotype to 'z' default is 0 remove SNPs where allels are non-ACTG, including 0
rule chrMT:
    input:
        bplink = expand(DATAIN + "/{{sample}}.{ext}", ext=BPLINK),
        fam = DATAOUT + "/{sample}/{sample}_maleOnly.fam"
    output:
        expand(DATAOUT + "/{{sample}}/chrMT_{{sample}}.{ext}", ext=BPLINK)
    params:
        inFile = DATAIN + "/{sample}",
        out = DATAOUT + "/{sample}/chrMT_{sample}"
    shell:
        'plink --bfile {params.inFile} --fam {input.fam} \
        --chr 26 --output-chr 26 --missing-genotype z --snps-only just-acgt \
        --keep-allele-order --make-bed --out {params.out}'

rule yri2rcrs_flip:
    input:
        script = 'scripts/yri_to_rcrs_flip.R',
        referenceSnps = 'ReferencePanel/ReferenceSNPs.txt',
        bim = DATAOUT + "/{sample}/chrMT_{sample}.bim"
    output:
        bim = DATAOUT + "/{sample}/chrMT_{sample}_rcrsFlipped.bim"
    shell:
        'Rscript {input.script} {input.referenceSnps} {input.bim} {output.bim}'

rule bplink2oxford:
    input:
        bplink = expand(DATAOUT + "/{{sample}}/chrMT_{{sample}}.{ext}", ext=BPLINK),
        bim = DATAOUT + "/{sample}/chrMT_{sample}_rcrsFlipped.bim"
    output:
        expand(DATAOUT + "/{{sample}}/{{sample}}.{ext}", ext=OXFORD)
    params:
        inFile = DATAOUT + "/{sample}/chrMT_{sample}",
        out = DATAOUT + "/{sample}/{sample}"
    shell:
        'plink --bfile {params.inFile} --bim {input.bim} \
        --recode oxford --keep-allele-order --out {params.out}'

rule bplink2plink:
    input:
        bplink = expand(DATAOUT + "/{{sample}}/chrMT_{{sample}}.{ext}", ext=BPLINK),
        bim = DATAOUT + "/{sample}/chrMT_{sample}_rcrsFlipped.bim"
    output:
        expand(DATAOUT + "/{{sample}}/{{sample}}_typedOnly.{ext}", ext=PLINK)
    params:
        inFile = DATAOUT + "/{sample}/chrMT_{sample}",
        out = DATAOUT + "/{sample}/{sample}_typedOnly"
    shell:
        'plink --bfile {params.inFile} --bim {input.bim} \
        --recode --keep-allele-order --out {params.out}'

rule Impute2:
    input:
        m = expand('{RefData}/MtMap.txt', RefData=REFDATA),
        h = expand('{RefData}/ReferencePanel.hap.gz', RefData=REFDATA),
        l = expand('{RefData}/ReferencePanel.legend.gz', RefData=REFDATA),
        g = DATAOUT + "/{sample}/{sample}.gen",
        sample = DATAOUT + "/{sample}/{sample}.sample",
    output:
        DATAOUT + "/{sample}/{sample}_imputed",
        DATAOUT + "/{sample}/{sample}_imputed_samples",
        DATAOUT + "/{sample}/{sample}_imputed_info"
    params:
        out = DATAOUT + "/{sample}/{sample}_imputed"
    shell:
        'impute2 -chrX -m {input.m} -h {input.h} -l {input.l} -g {input.g} \
        -sample_g {input.sample} -int 1 16569 -Ne 20000 -o {params.out}'

rule FixChromName:
    input:
        InFile = DATAOUT + "/{sample}/{sample}_imputed"
    output:
        OutFile = DATAOUT + "/{sample}/{sample}_imputed_ChromFixed"
    shell:"""
        awk '{{$1 = "26"; print}}' {input.InFile} > {output.OutFile}
    """

rule oxford2bed:
    input:
        gen = DATAOUT + "/{sample}/{sample}_imputed_ChromFixed",
        sample = DATAOUT + "/{sample}/{sample}_imputed_samples"
    output:
        expand(DATAOUT + "/{{sample}}/Imputed_{{sample}}.{ext}", ext=BPLINK)
    params:
        out = DATAOUT + "/{sample}/Imputed_{sample}"
    shell:
        'plink --gen {input.gen} --sample {input.sample} --hard-call-threshold 0.49 \
        --keep-allele-order --make-bed --output-chr 26 --out {params.out}'

rule oxford2ped:
    input:
        gen = DATAOUT + "/{sample}/{sample}_imputed_ChromFixed",
        sample = DATAOUT + "/{sample}/{sample}_imputed_samples"
    output:
        expand(DATAOUT + "/{{sample}}/Imputed_{{sample}}.{ext}", ext=PLINK)
    params:
        out = DATAOUT + "/{sample}/Imputed_{sample}"
    shell:
        'plink --gen {input.gen} --sample {input.sample} --hard-call-threshold 0.49 \
        --keep-allele-order --output-chr 26 --recode --out {params.out}'

rule oxford2vcf:
    input:
        gen = DATAOUT + "/{sample}/{sample}_imputed_ChromFixed",
        sample = DATAOUT + "/{sample}/{sample}_imputed_samples"
    output:
        DATAOUT + "/{sample}/Imputed_{sample}.vcf"
    params:
        out = DATAOUT + "/{sample}/Imputed_{sample}"
    shell:
        'plink --gen {input.gen} --sample {input.sample} --hard-call-threshold 0.49 \
        --keep-allele-order --output-chr 26 --recode vcf --out {params.out}'

rule html_Report:
    input:
        script = 'scripts/MT_imputation_report.Rmd',
        typ_map = DATAOUT + "/{sample}/{sample}_typedOnly.map",
        typ_ped = DATAOUT + "/{sample}/{sample}_typedOnly.ped",
        imp_map = DATAOUT + "/{sample}/Imputed_{sample}.map",
        imp_ped = DATAOUT + "/{sample}/Imputed_{sample}.ped",
        imp_info = DATAOUT + "/{sample}/{sample}_imputed_info",
    output:
        DATAOUT + "/{sample}/stats/{sample}_mtImputed_QC.html"
    params:
        rwd = RWD,
        output_dir = DATAOUT + "/{sample}/stats",
        info_cut = INFOCUT,
        sample = '{sample}'

    shell:
        "R -e 'rmarkdown::render("
        """"{input.script}", output_file = "{output}", output_dir = "{params.output_dir}", \
params = list(rwd = "{params.rwd}", info_cut = "{params.info_cut}", sample = "{params.sample}", \
typ_map = "{input.typ_map}", typ_ped = "{input.typ_ped}", \
imp_map = "{input.imp_map}", imp_ped = "{input.imp_ped}", \
imp_info = "{input.imp_info}", out = "{params.output_dir}"))' --slave
        """
