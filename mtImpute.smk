'''Snakefile for MitoImpute Version 0.1'''

configfile: 'mtImpute_config.yaml'
SAMPLE = config['SAMPLE']
DATAIN = config['DATAIN']
DATAOUT = config['DATAOUT']
REFDATA = config['REFDATA']
INFOCUT = config['INFOCUT']
ITER = config['ITER']
BURNIN = config['BURNIN']
KHAP = config['KHAP']

BPLINK = ["bed", "bim", "fam"]
PLINK = ["map", "ped"]
OXFORD = ["gen", "sample"]

## For running on cluster
#mkdir .snakejob; snakejob -s mtImpute.smk -j 100 --until oxford2vcf oxford2ped oxford2bed bplink2plink --max-jobs-per-second 1 --keep-going
shell.prefix('module load plink/1.90b6.1 impute2 R/3.4.3; ')

rule all:
    input:
        expand(DATAOUT + "/{sample}/Imputed_{sample}.{ext}", ext=BPLINK, sample=SAMPLE),
        expand(DATAOUT + "/{sample}/Imputed_{sample}.{ext}", ext=PLINK, sample=SAMPLE),
        expand(DATAOUT + "/{sample}/Imputed_{sample}.vcf", sample=SAMPLE),
        expand(DATAOUT + "/{sample}/stats/{sample}_Info.png", sample=SAMPLE),
        expand(DATAOUT + "/{sample}/stats/{sample}_InfoHiMC.png", sample=SAMPLE),
        expand(DATAOUT + "/{sample}/stats/{sample}_Haplogroups.txt", sample=SAMPLE),
        expand(DATAOUT + "/{sample}/stats/{sample}_HaplogroupMatch.png", sample=SAMPLE)

rule SexFam:
    input: DATAIN + "/{sample}.fam"
    output: DATAOUT + "/{sample}/{sample}_maleOnly.fam"
    shell: """
        awk '{{$5 = "1"; print}}' {input} > {output}
    """

## Extract SNPs on mitochondrial genome (26)
## Remove SNPs where allels are non-ACTG
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
        --chr 26 --output-chr 26 --snps-only just-acgt \
        --keep-allele-order --make-bed --out {params.out}'

## Check if mtSNPs are mapped to rCRS, if not liftover
rule yri2rcrs_flip:
    input:
        script = 'scripts/yri_to_rcrs_flip.R',
        referenceSnps = 'ReferencePanel/ReferenceSNPs.txt',
        bim = DATAOUT + "/{sample}/chrMT_{sample}.bim"
    output:
        bim = DATAOUT + "/{sample}/chrMT_{sample}_rcrsFlipped.bim"
    shell:
        'Rscript {input.script} {input.referenceSnps} {input.bim} {output.bim}'

## Convert sample binary plink files to oxford format
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

## Convert sample binary plink files to .map/.ped files
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

## Use Impute to impute mtSNPs
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
        out = DATAOUT + "/{sample}/{sample}_imputed",
        iter = ITER,
        burnin = BURNIN,
        khap = KHAP
    shell:
        'impute2 -chrX -m {input.m} -h {input.h} -l {input.l} -g {input.g} \
        -sample_g {input.sample} -int 1 16569 -Ne 20000 -o {params.out} \
        -iter {params.iter} -burnin {params.burnin} -k_hap {params.khap}'

## Change MT chromsome name to '26'
rule FixChromName:
    input:
        InFile = DATAOUT + "/{sample}/{sample}_imputed"
    output:
        OutFile = DATAOUT + "/{sample}/{sample}_imputed_ChromFixed"
    shell:"""
        awk '{{$1 = "26"; print}}' {input.InFile} > {output.OutFile}
    """

## Convert Oxford files to binary Plink files
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

## Convert Oxford files to .map/.ped files
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

## Convert Oxford files to .vcf files
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

rule plots:
    input:
        script = 'scripts/plots.R',
        typ_map = DATAOUT + "/{sample}/{sample}_typedOnly.map",
        typ_ped = DATAOUT + "/{sample}/{sample}_typedOnly.ped",
        imp_map = DATAOUT + "/{sample}/Imputed_{sample}.map",
        imp_ped = DATAOUT + "/{sample}/Imputed_{sample}.ped",
        imp_info = DATAOUT + "/{sample}/{sample}_imputed_info",
    output:
        InfoPlot = DATAOUT + "/{sample}/stats/{sample}_Info.png",
        InfoHimcPlot = DATAOUT + "/{sample}/stats/{sample}_InfoHiMC.png",
        haplgroups = DATAOUT + "/{sample}/stats/{sample}_Haplogroups.txt",
        HaplogroupPlot = DATAOUT + "/{sample}/stats/{sample}_HaplogroupMatch.png"
    params:
        output_dir = DATAOUT + "/{sample}/stats",
        sample = '{sample}',
        info_cut = INFOCUT,
    shell:
        'Rscript {input.script} {input.typ_map} {input.typ_ped} \
        {input.imp_map} {input.imp_ped} {input.imp_info} \
        {params.output_dir} {params.sample} {params.info_cut}'
