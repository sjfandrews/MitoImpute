'''Snakefile for MitoImpute Version 0.1'''

BPLINK = ["bed", "bim", "fam"]
PLINK = ["map", "ped"]
OXFORD = ["gen", "sample"]
SAMPLE = 'ADNI1_MTsnps'
DATAIN = '/Users/sheaandrews/LOAD_minerva/dummy/shea'

rule all:
    input:
        expand("DerivedData/SamplePanel/Imputed_{sample}.{ext}", ext=BPLINK, sample = SAMPLE),
        expand("DerivedData/SamplePanel/Imputed_{sample}.{ext}", ext=PLINK, sample = SAMPLE),
        expand("DerivedData/SamplePanel/Imputed_{sample}.vcf", sample = SAMPLE)

rule SexFam:
    input: expand("{DataIn}/{{sample}}.fam", DataIn=DATAIN)
    output: "DerivedData/SamplePanel/{sample}_maleOnly.fam"
    shell: """
        awk '{{$5 = "1"; print}}' {input} > {output}
    """

rule plink2oxford:
    input:
        bplink = expand("{DataIn}/{{sample}}.{ext}", ext=BPLINK, DataIn=DATAIN),
        fam = "DerivedData/SamplePanel/{sample}_maleOnly.fam"
    output:
        expand("DerivedData/SamplePanel/{{sample}}.{ext}", ext=OXFORD)
    params:
        inFile = expand("{DataIn}/{sample}", DataIn=DATAIN, sample = SAMPLE),
        out = "DerivedData/SamplePanel/{sample}"
    shell:
        'plink --bfile {params.inFile} --fam {input.fam} \
        --recode oxford --chr 26 --output-chr 26 --keep-allele-order --out {params.out}'

rule Impute2:
    input:
        m = 'DerivedData/ReferencePanel/MtMap.txt',
        h = 'DerivedData/ReferencePanel/ReferencePanel.hap.gz',
        l = 'DerivedData/ReferencePanel/ReferencePanel.legend.gz',
        g = 'DerivedData/SamplePanel/{sample}.gen',
        strand_g_ref = 'DerivedData/ReferencePanel/MtStrand.txt',
        sample = 'DerivedData/SamplePanel/{sample}.sample',
    output:
        'DerivedData/SamplePanel/{sample}_imputed',
        'DerivedData/SamplePanel/{sample}_imputed_samples'
    params:
        out = 'DerivedData/SamplePanel/{sample}_imputed'
    shell:
        'impute2 -chrX -m {input.m} -h {input.h} -l {input.l} -g {input.g} -strand_g_ref {input.strand_g_ref} \
        -sample_g {input.sample} -int 1 16569 -Ne 20000 -o {params.out}'

rule FixChromName:
    input:
        InFile = 'DerivedData/SamplePanel/{sample}_imputed'
    output:
        OutFile = 'DerivedData/SamplePanel/{sample}_imputed_ChromFixed'
    shell:"""
        awk '{{$1 = "26"; print}}' {input.InFile} > {output.OutFile}
    """

rule oxford2bed:
    input:
        gen = 'DerivedData/SamplePanel/{sample}_imputed_ChromFixed',
        sample = 'DerivedData/SamplePanel/{sample}_imputed_samples'
    output:
        expand("DerivedData/SamplePanel/Imputed_{{sample}}.{ext}", ext=BPLINK)
    params:
        out = 'DerivedData/SamplePanel/Imputed_{sample}'
    shell:
        'plink --gen {input.gen} --sample {input.sample} --hard-call-threshold 0.49 \
        --keep-allele-order --make-bed --output-chr 26 --out {params.out}'

rule oxford2ped:
    input:
        gen = 'DerivedData/SamplePanel/{sample}_imputed_ChromFixed',
        sample = 'DerivedData/SamplePanel/{sample}_imputed_samples'
    output:
        expand("DerivedData/SamplePanel/Imputed_{{sample}}.{ext}", ext=PLINK)
    params:
        out = 'DerivedData/SamplePanel/Imputed_{sample}'
    shell:
        'plink --gen {input.gen} --sample {input.sample} --hard-call-threshold 0.49 \
        --keep-allele-order --output-chr 26 --recode --out {params.out}'

rule oxford2vcf:
    input:
        gen = 'DerivedData/SamplePanel/{sample}_imputed_ChromFixed',
        sample = 'DerivedData/SamplePanel/{sample}_imputed_samples'
    output:
        "DerivedData/SamplePanel/Imputed_{sample}.vcf"
    params:
        out = 'DerivedData/SamplePanel/Imputed_{sample}'
    shell:
        'plink --gen {input.gen} --sample {input.sample} --hard-call-threshold 0.49 \
        --keep-allele-order --output-chr 26 --recode vcf --out {params.out}'
