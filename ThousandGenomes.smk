'''Snakefile for Cleaning Thousand Genomes'''
# snakemake -s ThousandGenomes.smk
import os
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider

with open('data/platforms/Mt_platforms.txt', "r") as f:
    MtPlatforms = [x.rstrip() for x in f]
FTP = FTPRemoteProvider()

rule all:
    input:
        expand("DerivedData/ThousandGenomes/{MtPlatforms}/chrMT_1kg_{MtPlatforms}.{ext}", ext = ['gen.gz', 'samples'], MtPlatforms=MtPlatforms),
        expand("DerivedData/ThousandGenomes/{MtPlatforms}/chrMT_1kg_{MtPlatforms}.{ext}", ext = ['ped', 'map'], MtPlatforms=MtPlatforms)

# 1. Pull down 1000 genomes mitochondrial vcf file from ftp
rule Get1kgMT_vcf:
    input:
        vcf = FTP.remote("http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz", keep_local=True),
        tbi = FTP.remote("http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz.tbi", keep_local=True)
    output:
        vcf = "data/ThousandGenomes/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz",
        tbi = "data/ThousandGenomes/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz.tbi"
    shell:
        "mv {input.vcf} {output.vcf}; "
        "mv {input.tbi} {output.tbi}"

rule NormaliseVcf:
    input:
        vcf = "data/ThousandGenomes/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz",
        fasta = "data/ReferencePanel/rCRS.fasta"
    output:
        vcf = "DerivedData/ThousandGenomes/chrMT_1kg_norm.vcf.gz"
    shell:
        'bcftools norm -f {input.fasta} -m - {input.vcf} | bcftools view -V indels,mnps | bcftools norm -m + | bcftools +fill-tags -Oz -o {output.vcf}'

rule DecomposeVcf:
    input:
        vcf = "DerivedData/ThousandGenomes/chrMT_1kg_norm.vcf.gz",
    output:
        vcf = "DerivedData/ThousandGenomes/chrMT_1kg_norm_decomposed.vcf.gz"
    shell:
        'vt decompose {input.vcf} | bcftools +fill-tags -Oz -o {output.vcf}'

rule pickFirstAlt:
    input:
        vcf = "DerivedData/ThousandGenomes/chrMT_1kg_norm_decomposed.vcf.gz",
    output:
        vcf = "DerivedData/ThousandGenomes/chrMT_1kg_norm_decomposed_firstAlt.vcf.gz"
    shell:
        "./scripts/PYTHON/pickFirstAlt {input.vcf} | bgzip > {output.vcf}; "
        "bcftools index {output.vcf}"

## 6a. Extract sample names from Reference Panel
rule SampleNames1kg:
    input:
        "DerivedData/ThousandGenomes/chrMT_1kg_norm_decomposed_firstAlt.vcf.gz",
    output:
        "DerivedData/ThousandGenomes/SampleList1kg.txt",
    shell:
        'bcftools query -l {input} > {output}'

## 6b. Assign M sex label to reference Samples
rule SampleSex1kg:
    input:
        in_samples = "DerivedData/ThousandGenomes/SampleList1kg.txt",
        in_script = "scripts/R/assign_sex_label.R"
    output:
        "DerivedData/ThousandGenomes/SampleList1kg_sex.txt",
    shell:
        'Rscript {input.in_script} {input.in_samples} {output}'

rule ExtractPlatformMTsnps:
    input:
        MTSnps = 'data/platforms/{MtPlatforms}/{MtPlatforms}_MT_snps.txt',
        vcf_1kg = "DerivedData/ThousandGenomes/chrMT_1kg_norm_decomposed_firstAlt.vcf.gz"
    output:
        vcf = "DerivedData/ThousandGenomes/{MtPlatforms}/chrMT_1kg_{MtPlatforms}.vcf.gz"
    shell:
        "bcftools view -R {input.MTSnps} {input.vcf_1kg} -Oz -o {output.vcf}; "
        "bcftools index {output.vcf}"

rule vcf2gensample:
    input:
        vcf = "DerivedData/ThousandGenomes/{MtPlatforms}/chrMT_1kg_{MtPlatforms}.vcf.gz",
        sex = "DerivedData/ThousandGenomes/SampleList1kg_sex.txt"
    output:
        expand("DerivedData/ThousandGenomes/{{MtPlatforms}}/chrMT_1kg_{{MtPlatforms}}.{ext}", ext = ['gen.gz', 'samples'])
    params:
        out = "DerivedData/ThousandGenomes/{MtPlatforms}/chrMT_1kg_{MtPlatforms}"
    shell:
        "bcftools convert --gensample {params.out} {input.vcf} --sex {input.sex}"

rule Plink:
    input:
        vcf = "DerivedData/ThousandGenomes/{MtPlatforms}/chrMT_1kg_{MtPlatforms}.vcf.gz",
    output:
        expand("DerivedData/ThousandGenomes/{{MtPlatforms}}/chrMT_1kg_{{MtPlatforms}}.{ext}", ext = ['ped', 'map'])
    params:
        out = "DerivedData/ThousandGenomes/{MtPlatforms}/chrMT_1kg_{MtPlatforms}"
    shell:
        'plink --vcf {input.vcf} --recode --double-id --out {params.out}'
