'''Snakefile for MitoImpute Version 0.1'''
# snakemake -s mtImpute.smk
# snakemake -s mtImpute.smk --dag | dot -Tsvg > dag_mtImpute.svg

import os

configfile: 'mtImpute_config.yaml'
SAMPLE = config['SAMPLE']
DATAIN = config['DATAIN']
DATAOUT = config['DATAOUT']
REFDATA = config['REFDATA']

BPLINK = ["bed", "bim", "fam"]
PLINK = ["map", "ped"]
OXFORD = ["gen", "sample"]
RWD = os.getcwd()

## For running on cluster
#shell.prefix('module load plink/1.90 impute2 R/3.4.3; ')

rule all:
    input:
        expand(DATAOUT + "/Imputed_{sample}.{ext}", ext=BPLINK, sample=SAMPLE),
        expand(DATAOUT + "/Imputed_{sample}.{ext}", ext=PLINK, sample=SAMPLE),
        expand(DATAOUT + "/Imputed_{sample}.vcf", sample=SAMPLE),
        expand(DATAOUT + "/stats/{sample}_mtImputed_QC.html", sample=SAMPLE),

rule SexFam:
    input: DATAIN + "/{sample}.fam"
    output: DATAOUT + "/{sample}_maleOnly.fam"
    shell: """
        awk '{{$5 = "1"; print}}' {input} > {output}
    """

rule plink2oxford:
    input:
        bplink = expand(DATAIN + "/{{sample}}.{ext}", ext=BPLINK),
        fam = DATAOUT + "/{sample}_maleOnly.fam"
    output:
        expand(DATAOUT + "/{{sample}}.{ext}", ext=OXFORD)
    params:
        inFile = expand(DATAIN + "/{sample}", sample = SAMPLE),
        out = DATAOUT + "/{sample}"
    shell:
        'plink --bfile {params.inFile} --fam {input.fam} \
        --recode oxford --chr 26 --output-chr 26 --keep-allele-order --out {params.out}'

rule bplink2plink:
    input:
        bplink = expand(DATAIN + "/{{sample}}.{ext}", ext=BPLINK),
        fam = DATAOUT + "/{sample}_maleOnly.fam"
    output:
        expand(DATAOUT + "/{{sample}}_typedOnly.{ext}", ext=PLINK)
    params:
        inFile = expand(DATAIN + "/{sample}", sample = SAMPLE),
        out = DATAOUT + "/{sample}_typedOnly"
    shell:
        'plink --bfile {params.inFile} --fam {input.fam} \
        --recode --chr 26 --output-chr 26 --keep-allele-order --out {params.out}'

rule Impute2:
    input:
        m = expand('{RefData}/MtMap.txt', RefData=REFDATA),
        h = expand('{RefData}/ReferencePanel.hap.gz', RefData=REFDATA),
        l = expand('{RefData}/ReferencePanel.legend.gz', RefData=REFDATA),
        g = DATAOUT + "/{sample}.gen",
        sample = DATAOUT + "/{sample}.sample",
    output:
        DATAOUT + "/{sample}_imputed",
        DATAOUT + "/{sample}_imputed_samples",
        DATAOUT + "/{sample}_imputed_info"
    params:
        out = DATAOUT + "/{sample}_imputed"
    shell:
        'impute2 -chrX -m {input.m} -h {input.h} -l {input.l} -g {input.g} \
        -sample_g {input.sample} -int 1 16569 -Ne 20000 -o {params.out}'

rule FixChromName:
    input:
        InFile = DATAOUT + "/{sample}_imputed"
    output:
        OutFile = DATAOUT + "/{sample}_imputed_ChromFixed"
    shell:"""
        awk '{{$1 = "26"; print}}' {input.InFile} > {output.OutFile}
    """

rule oxford2bed:
    input:
        gen = DATAOUT + "/{sample}_imputed_ChromFixed",
        sample = DATAOUT + "/{sample}_imputed_samples"
    output:
        expand(DATAOUT + "/Imputed_{{sample}}.{ext}", ext=BPLINK)
    params:
        out = DATAOUT + "/Imputed_{sample}"
    shell:
        'plink --gen {input.gen} --sample {input.sample} --hard-call-threshold 0.49 \
        --keep-allele-order --make-bed --output-chr 26 --out {params.out}'

rule oxford2ped:
    input:
        gen = DATAOUT + "/{sample}_imputed_ChromFixed",
        sample = DATAOUT + "/{sample}_imputed_samples"
    output:
        expand(DATAOUT + "/Imputed_{{sample}}.{ext}", ext=PLINK)
    params:
        out = DATAOUT + "/Imputed_{sample}"
    shell:
        'plink --gen {input.gen} --sample {input.sample} --hard-call-threshold 0.49 \
        --keep-allele-order --output-chr 26 --recode --out {params.out}'

rule oxford2vcf:
    input:
        gen = DATAOUT + "/{sample}_imputed_ChromFixed",
        sample = DATAOUT + "/{sample}_imputed_samples"
    output:
        DATAOUT + "/Imputed_{sample}.vcf"
    params:
        out = DATAOUT + "/Imputed_{sample}"
    shell:
        'plink --gen {input.gen} --sample {input.sample} --hard-call-threshold 0.49 \
        --keep-allele-order --output-chr 26 --recode vcf --out {params.out}'

rule Imputation_QC_Report:
    input:
        script = 'scripts/R/MT_imputation_QC_examples.Rmd',
        typ_map = DATAOUT + "/{sample}_typedOnly.map",
        typ_ped = DATAOUT + "/{sample}_typedOnly.ped",
        imp_map = DATAOUT + "/Imputed_{sample}.map",
        imp_ped = DATAOUT + "/Imputed_{sample}.ped",
        imp_info = DATAOUT + "/{sample}_imputed_info",
    output:
        DATAOUT + "/stats/{sample}_mtImputed_QC.html"
    params:
        rwd = RWD,
        output_dir = DATAOUT + "/stats",
        info_cut = '0'

    shell:
        "R -e 'rmarkdown::render("
        """"{input.script}", output_file = "{output}", output_dir = "{params.output_dir}", \
params = list(rwd = "{params.rwd}", info_cut = "{params.info_cut}", \
typ_map = "{input.typ_map}", typ_ped = "{input.typ_ped}", \
imp_map = "{input.imp_map}", imp_ped = "{input.imp_ped}", \
imp_info = "{input.imp_info}"))' --slave
        """
