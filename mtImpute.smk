'''Snakefile for MitoImpute Version 0.1'''
# snakemake -s mtImpute.smk
# snakemake -s mtImpute.smk --dag | dot -Tsvg > dag_mtImpute.svg

configfile: 'mtImpute_config.yaml'
SAMPLE = config['SAMPLE']
DATAIN = config['DATAIN']
DATAOUT = config['DATAOUT']
REFDATA = config['REFDATA']

BPLINK = ["bed", "bim", "fam"]
PLINK = ["map", "ped"]
OXFORD = ["gen", "sample"]

rule all:
    input:
        expand("{DataOut}/Imputed_{sample}.{ext}", ext=BPLINK, sample=SAMPLE, DataOut=DATAOUT),
        expand("{DataOut}/Imputed_{sample}.{ext}", ext=PLINK, sample=SAMPLE, DataOut=DATAOUT),
        expand("{DataOut}/Imputed_{sample}.vcf", sample=SAMPLE, DataOut=DATAOUT)

rule SexFam:
    input: expand("{DataIn}/{{sample}}.fam", DataIn=DATAIN)
    output: "{DataOut}/{sample}_maleOnly.fam"
    shell: """
        awk '{{$5 = "1"; print}}' {input} > {output}
    """

rule plink2oxford:
    input:
        bplink = expand("{DataIn}/{{sample}}.{ext}", ext=BPLINK, DataIn=DATAIN),
        fam = "{DataOut}/{sample}_maleOnly.fam"
    output:
        expand("{{DataOut}}/{{sample}}.{ext}", ext=OXFORD)
    params:
        inFile = expand("{DataIn}/{sample}", DataIn=DATAIN, sample = SAMPLE),
        out = "{DataOut}/{sample}"
    shell:
        'plink --bfile {params.inFile} --fam {input.fam} \
        --recode oxford --chr 26 --output-chr 26 --keep-allele-order --out {params.out}'

rule Impute2:
    input:
        m = expand('{RefData}/MtMap.txt', RefData=REFDATA),
        h = expand('{RefData}/ReferencePanel.hap.gz', RefData=REFDATA),
        l = expand('{RefData}/ReferencePanel.legend.gz', RefData=REFDATA),
        g = '{DataOut}/{sample}.gen',
        sample = '{DataOut}/{sample}.sample',
    output:
        '{DataOut}/{sample}_imputed',
        '{DataOut}/{sample}_imputed_samples'
    params:
        out = '{DataOut}/{sample}_imputed'
    shell:
        'impute2 -chrX -m {input.m} -h {input.h} -l {input.l} -g {input.g} \
        -sample_g {input.sample} -int 1 16569 -Ne 20000 -o {params.out}'

rule FixChromName:
    input:
        InFile = '{DataOut}/{sample}_imputed'
    output:
        OutFile = '{DataOut}/{sample}_imputed_ChromFixed'
    shell:"""
        awk '{{$1 = "26"; print}}' {input.InFile} > {output.OutFile}
    """

rule oxford2bed:
    input:
        gen = '{DataOut}/{sample}_imputed_ChromFixed',
        sample = '{DataOut}/{sample}_imputed_samples'
    output:
        expand("{{DataOut}}/Imputed_{{sample}}.{ext}", ext=BPLINK)
    params:
        out = '{DataOut}/Imputed_{sample}'
    shell:
        'plink --gen {input.gen} --sample {input.sample} --hard-call-threshold 0.49 \
        --keep-allele-order --make-bed --output-chr 26 --out {params.out}'

rule oxford2ped:
    input:
        gen = '{DataOut}/{sample}_imputed_ChromFixed',
        sample = '{DataOut}/{sample}_imputed_samples'
    output:
        expand("{{DataOut}}/Imputed_{{sample}}.{ext}", ext=PLINK)
    params:
        out = '{DataOut}/Imputed_{sample}'
    shell:
        'plink --gen {input.gen} --sample {input.sample} --hard-call-threshold 0.49 \
        --keep-allele-order --output-chr 26 --recode --out {params.out}'

rule oxford2vcf:
    input:
        gen = '{DataOut}/{sample}_imputed_ChromFixed',
        sample = '{DataOut}/{sample}_imputed_samples'
    output:
        "{DataOut}/Imputed_{sample}.vcf"
    params:
        out = '{DataOut}/Imputed_{sample}'
    shell:
        'plink --gen {input.gen} --sample {input.sample} --hard-call-threshold 0.49 \
        --keep-allele-order --output-chr 26 --recode vcf --out {params.out}'
