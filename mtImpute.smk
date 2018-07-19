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
RWD = os.getcwd()

rule all:
    input:
        expand("{DataOut}/Imputed_{sample}.{ext}", ext=BPLINK, sample=SAMPLE, DataOut=DATAOUT),
        expand("{DataOut}/Imputed_{sample}.{ext}", ext=PLINK, sample=SAMPLE, DataOut=DATAOUT),
        expand("{DataOut}/Imputed_{sample}.vcf", sample=SAMPLE, DataOut=DATAOUT),
        expand("{DataOut}/stats/{sample}_mtImputed_QC.html", sample=SAMPLE, DataOut=DATAOUT),

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

rule bplink2plink:
    input:
        bplink = expand("{DataIn}/{{sample}}.{ext}", ext=BPLINK, DataIn=DATAIN),
        fam = "{DataOut}/{sample}_maleOnly.fam"
    output:
        expand("{{DataOut}}/{{sample}}_typedOnly.{ext}", ext=PLINK)
    params:
        inFile = expand("{DataIn}/{sample}", DataIn=DATAIN, sample = SAMPLE),
        out = "{DataOut}/{sample}_typedOnly"
    shell:
        'plink --bfile {params.inFile} --fam {input.fam} \
        --recode --chr 26 --output-chr 26 --keep-allele-order --out {params.out}'

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

rule Imputation_QC_Report:
    input:
        script = 'scripts/R/MT_imputation_QC_examples.Rmd',
        typ_map = "{DataOut}/{sample}_typedOnly.map",
        typ_ped = "{DataOut}/{sample}_typedOnly.ped",
        imp_map = '{DataOut}/Imputed_{sample}.map',
        imp_ped = '{DataOut}/Imputed_{sample}.ped',
        imp_info = '{DataOut}/{sample}_imputed_info',
    output:
        "{DataOut}/stats/{sample}_mtImputed_QC.html"
    params:
        rwd = RWD,
        output_dir = "{DataOut}/stats",
        info_cut = '0'

    shell:
        "R -e 'rmarkdown::render("
        """"{input.script}", output_file = "{output}", output_dir = "{params.output_dir}", \
params = list(rwd = "{params.rwd}", info_cut = "{params.info_cut}", \
typ_map = "{input.typ_map}", typ_ped = "{input.typ_ped}", \
imp_map = "{input.imp_map}", imp_ped = "{input.imp_ped}", \
imp_info = "{input.imp_info}"))' --slave
        """
