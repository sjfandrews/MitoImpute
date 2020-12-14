'''Snakefile for Constructing example data sets'''
# snakemake -s workflow/rules/ExampleData.smk
# snakemake -s workflow/rules/ExampleData.smk --dag | dot -Tsvg > resources/example/dag_Example.svg

##===========================##
##  Example Reference Panel
##===========================##

from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

FTP = FTPRemoteProvider()
HTTP = HTTPRemoteProvider()

rule all:
    input:
        "resources/example/ReferencePanel/MtMap.txt",
        "resources/example/ReferencePanel/MtStrand.txt",
        expand("resources/example/ReferencePanel/ReferencePanel.{ext}", ext = ['hap.gz', 'legend.gz']),
        expand("resources/example/ReferencePanel/ReferencePanel.{ext}", ext = ['ped', 'map']),
        expand("resources/example/ReferencePanel/ReferencePanel.{ext}", ext = ['gen.gz', 'samples']),
        expand("resources/example/SamplePanel/ExampleSamplePanel.{ext}", ext = ['ped', 'map']),
        expand("resources/example/SamplePanel/ExampleSamplePanel.{ext}", ext = ['bed', 'bim', 'fam']),

##============================================================================##
##  Example Reference Panel
##============================================================================##

rule ExampleReferenceSamples:
    input: ref = 'resources/ReferencePanel/ReferencePanel_haplogroups.txt'
    output: out = "resources/example/RefSampleList.txt"
    script: '../scripts/ExampleRefSamples.R'

## 6b. Assign M sex label to reference Samples
rule RefSampleSex:
    input: rules.ExampleReferenceSamples.output.out
    output: "resources/example/ReferencePanel/RefSampleList_sex.txt"
    shell: '''awk '{{$5 = "M"; print}}' {input} > {output}'''

rule ExtractExampleReference:
    input:
        inVCF = "resources/ReferencePanel/ReferencePanel.vcf.gz",
        inSamples = rules.ExampleReferenceSamples.output.out
    output:
        vcf = "resources/example/ReferencePanel/ReferencePanel.vcf.gz"
    conda:
        "../envs/vcf.yaml"
    shell:
        """
        bcftools view -S {input.inSamples} --force-samples -Oz -o {output} {input.inVCF};
        bcftools index {output}
        """

## 7. Convert to Oxford format
rule Oxford:
    input:
        vcf = rules.ExtractExampleReference.output.vcf,
        sex = "resources/example/ReferencePanel/RefSampleList_sex.txt"
    output:
        expand("resources/example/ReferencePanel/ReferencePanel.{ext}", ext = ['hap.gz', 'legend.gz'])
    params:
        out = "resources/example/ReferencePanel/ReferencePanel"
    conda:
        "../envs/vcf.yaml"
    shell:
        'bcftools convert --haplegendsample {params.out} {input.vcf} --sex {input.sex} --haploid2diploid'

## 8. Generate .ped and .map files
rule Plink:
    input: vcf = rules.ExtractExampleReference.output.vcf,
    output:
        expand("resources/example/ReferencePanel/ReferencePanel.{ext}", ext = ['ped', 'map'])
    params:
        out = "resources/example/ReferencePanel/ReferencePanel"
    conda:
        "../envs/plink.yaml"
    shell:
        'plink --vcf {input} --recode --double-id --out {params.out}'

## 9. Generate .gen and .sample files
rule GenSample:
    input:
        vcf = rules.ExtractExampleReference.output.vcf,
        sex = "resources/example/ReferencePanel/RefSampleList_sex.txt"
    output:
        expand("resources/example/ReferencePanel/ReferencePanel.{ext}", ext = ['gen.gz', 'samples'])
    params:
        out = "resources/example/ReferencePanel/ReferencePanel"
    conda:
        "../envs/vcf.yaml"
    shell:
        'bcftools convert --gensample {params.out} {input.vcf} --sex {input.sex}'

## 10. Construct .map file for IMPUTE2
rule MakeMapFile:
    input:
        vcf = rules.ExtractExampleReference.output.vcf,
    output:
        map = "resources/example/ReferencePanel/MtMap.txt",
        strand = "resources/example/ReferencePanel/MtStrand.txt"
    script:
        '../scripts/mt_recombination_map.R'

##===========================##
##  Example Sample Panel
##===========================##

# 1. Pull down 1000 genomes mitochondrial vcf file from ftp
TG_releasedir = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/'
TG_vcf = 'ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz'

rule Get1kgMT_vcf:
    input:
        vcf = FTP.remote(TG_releasedir + TG_vcf),
        tbi = FTP.remote(TG_releasedir + TG_vcf + '.tbi'),
        info = FTP.remote('ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_sample_info.xlsx')
    output:
        vcf = temp('resources/example/data/ThousandGenomes/' + TG_vcf),
        tbi = temp('resources/example/data/ThousandGenomes/' + TG_vcf + '.tbi'),
        info = temp('resources/example/data/ThousandGenomes/20130606_sample_info.xlsx')
    shell:
        "mv {input.vcf} {output.vcf}; "
        "mv {input.tbi} {output.tbi}; "
        "mv {input.info} {output.info}"

rule NormaliseVcf:
    input:
        vcf = rules.Get1kgMT_vcf.output.vcf,
        fasta = "resources/example/ReferencePanel/rCRS.fasta"
    output: temp("resources/example/data/ThousandGenomes/chrMT_1kg_norm_firstAlt.vcf.gz")
    conda:
        "../envs/vcf.yaml"
    shell:
        """
        vt decompose {input.vcf} | bcftools norm -f {input.fasta} | bcftools view -v snps | \
           bcftools norm -d all | bcftools +fill-tags -o {output} -Oz;
           vt index {output}
        """

rule StrandFiles:
    input:
        HTTP.remote("https://www.well.ox.ac.uk/~wrayner/strand/Human610-Quadv1_B-b37-strand.zip", allow_redirects=True)
    output:
        temp('resources/example/data/Human610-Quadv1_B-b37.strand')
    params:
        directory = 'resources/example/data'
    shell:
        "unzip {input} -d {params.directory} *.strand; "

rule StrandFilesMT:
    input:
        strand = 'resources/example/data/Human610-Quadv1_B-b37.strand'
    output:
        out = temp('resources/example/data/platforms/Human610-Quadv1_B-b37_MT_snps.txt')
    script:
        '../scripts/StrandFiles_ExtractMTsnps.R'

rule ExtractPlatformMTsnps:
    input:
        MTSnps = 'resources/example/data/platforms/Human610-Quadv1_B-b37_MT_snps.txt',
        vcf_1kg = rules.NormaliseVcf.output
    output: temp("resources/example/data/ThousandGenomes/chrMT_1kg_Human610-Quadv1_B-b37.vcf.gz")
    conda:
        "../envs/vcf.yaml"
    shell:
        "bcftools view -R {input.MTSnps} {input.vcf_1kg} -Oz -o {output}; "
        "bcftools index {output}"

rule ExampleSampleIDs:
    input: info = rules.Get1kgMT_vcf.output.info,
    output: out = "resources/example/SampleIDList.txt"
    conda:
        "../envs/r.yaml"
    script: '../scripts/ExampleSampleIDs.R'

rule ExampleSamplePanel:
    input:
        inVCF = "resources/example/data/ThousandGenomes/chrMT_1kg_Human610-Quadv1_B-b37.vcf.gz",
        inSamples = rules.ExampleSampleIDs.output.out
    output:
        "resources/example/SamplePanel/ExampleSamplePanel.vcf.gz"
    conda:
        "../envs/vcf.yaml"
    shell:
        "bcftools view -S {input.inSamples} --force-samples -Oz -o {output} {input.inVCF}; "
        "bcftools index {output}"

rule vcf2Plink:
    input:
        vcf = "resources/example/SamplePanel/ExampleSamplePanel.vcf.gz",
    output:
        expand("resources/example/SamplePanel/ExampleSamplePanel.{ext}", ext = ['ped', 'map'])
    params:
        out = "resources/example/SamplePanel/ExampleSamplePanel"
    conda:
        "../envs/plink.yaml"
    shell:
        'plink --vcf {input.vcf} --recode --double-id --out {params.out}'

rule vcf2BPlink:
    input:
        vcf = "resources/example/SamplePanel/ExampleSamplePanel.vcf.gz",
    output:
        expand("resources/example/SamplePanel/ExampleSamplePanel.{ext}", ext = ['bed', 'bim', 'fam'])
    params:
        out = "resources/example/SamplePanel/ExampleSamplePanel"
    conda:
        "../envs/plink.yaml"
    shell:
        'plink --vcf {input.vcf} --make-bed --double-id --out {params.out}'
