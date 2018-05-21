'''Snakefile for Cleaning Thousand Genomes'''

from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()

rule all:
    input:
        "DerivedData/ThousandGenomes/chrMT_1kg_norm_decomposed_firstAlt.vcf.gz"

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
        './scripts/PYTHON/pickFirstAlt {input.vcf} | bgzip > {output.vcf}'
