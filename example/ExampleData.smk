# snakemake -s example/ExampleData.smk
# snakemake -s example/ExampleData.smk --dag | dot -Tsvg > dag_Example.svg

##===========================##
##  Example Reference Panel
##===========================##

from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

shell.executable("/bin/bash")

FTP = FTPRemoteProvider()
HTTP = HTTPRemoteProvider()

rule all:
    input:
        expand("example/ReferencePanel/ReferencePanel.{ext}", ext = ['hap.gz', 'legend.gz']),
        expand("example/ReferencePanel/ReferencePanel.{ext}", ext = ['ped', 'map']),
        expand("example/ReferencePanel/ReferencePanel.{ext}", ext = ['gen.gz', 'samples']),
        "example/ReferencePanel/MtMap.txt",
        "example/ReferencePanel/MtStrand.txt",
        expand("example/SamplePanel/ExampleSamplePanel.{ext}", ext = ['ped', 'map']),
        expand("example/SamplePanel/ExampleSamplePanel.{ext}", ext = ['bed', 'bim', 'fam']),


## 6b. Assign M sex label to reference Samples
rule RefSampleSex:
    input: "example/ExampleReferenceSamples.txt"
    output: "example/ReferencePanel/RefSampleList_sex.txt"
    shell: '''awk '{{$5 = "1"; print}}' {input} > {output}'''

rule ExampleReferenceSamples:
    input:
        inVCF = "DerivedData/ReferencePanel/ReferencePanel_highQual_filtered.vcf.gz",
        inSamples = "example/ExampleReferenceSamples.txt"
    output:
        "example/ReferencePanel/ReferencePanel.vcf.gz"
    run:
        shell('bcftools view -S {input.inSamples} -Oz -o {output} {input.inVCF}')
        shell('bcftools index {output}')


## 7. Convert to Oxford format
rule Oxford:
    input:
        vcf = "example/ReferencePanel/ReferencePanel.vcf.gz",
        sex = "example/ReferencePanel/RefSampleList_sex.txt"
    output:
        expand("example/ReferencePanel/ReferencePanel.{ext}", ext = ['hap.gz', 'legend.gz'])
    params:
        out = "example/ReferencePanel/ReferencePanel"
    shell:
        'bcftools convert --haplegendsample {params.out} {input.vcf} --sex {input.sex}'

## 8. Generate .ped and .map files
rule Plink:
    input: "example/ReferencePanel/ReferencePanel.vcf.gz",
    output:
        expand("example/ReferencePanel/ReferencePanel.{ext}", ext = ['ped', 'map'])
    params:
        out = "example/ReferencePanel/ReferencePanel"
    shell:
        'plink --vcf {input} --recode --double-id --out {params.out}'

## 9. Generate .gen and .sample files
rule GenSample:
    input:
        vcf = "example/ReferencePanel/ReferencePanel.vcf.gz",
        sex = "example/ReferencePanel/RefSampleList_sex.txt"
    output:
        expand("example/ReferencePanel/ReferencePanel.{ext}", ext = ['gen.gz', 'samples'])
    params:
        out = "example/ReferencePanel/ReferencePanel"
    shell:
        'bcftools convert --gensample {params.out} {input.vcf} --sex {input.sex}'

## 10. Construct .map file for IMPUTE2
rule MakeMapFile:
    input:
        script = "scripts/mt_recombination_map.R",
        panel = "example/ReferencePanel/ReferencePanel.vcf.gz"
    output:
        map = "example/ReferencePanel/MtMap.txt",
        strand = "example/ReferencePanel/MtStrand.txt"
    shell:
        'Rscript {input.script} {input.panel} {output.map} {output.strand}'

##===========================##
##  Example Sample Panel
##===========================##

# 1. Pull down 1000 genomes mitochondrial vcf file from ftp
TG_releasedir = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/'
TG_vcf = 'ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz'

rule Get1kgMT_vcf:
    input:
        vcf = FTP.remote(TG_releasedir + TG_vcf),
        tbi = FTP.remote(TG_releasedir + TG_vcf + '.tbi')
    output:
        vcf = temp('example/data/ThousandGenomes/' + TG_vcf),
        tbi = temp('example/data/ThousandGenomes/' + TG_vcf + '.tbi')
    shell:
        "mv {input.vcf} {output.vcf}; "
        "mv {input.tbi} {output.tbi}"

rule NormaliseVcf:
    input:
        vcf = rules.Get1kgMT_vcf.output.vcf,
        fasta = "example/ReferencePanel/rCRS.fasta"
    output: "example/data/ThousandGenomes/chrMT_1kg_norm_firstAlt.vcf.gz"
    shell:
        '''
vt decompose {input.vcf} | bcftools norm -f {input.fasta} | bcftools view -v snps | \
   bcftools norm -d all | bcftools +fill-tags -o {output} -Oz
'''

rule StrandFiles:
    input:
        HTTP.remote("https://www.well.ox.ac.uk/~wrayner/strand/Human610-Quadv1_B-b37-strand.zip", allow_redirects=True)
    output:
        temp('example/data/Human610-Quadv1_B-b37.strand')
    params:
        directory = 'example/data'
    shell:
        "unzip {input} -d {params.directory} *.strand; "

rule StrandFilesMT:
    input:
        script = 'scripts/StrandFiles_ExtractMTsnps.R',
        strand = 'example/data/Human610-Quadv1_B-b37.strand'
    output:
        out = temp('example/data/platforms/Human610-Quadv1_B-b37_MT_snps.txt')
    shell:
        'Rscript {input.script} {input.strand} {output.out}'

rule ExtractPlatformMTsnps:
    input:
        MTSnps = 'example/data/platforms/Human610-Quadv1_B-b37_MT_snps.txt',
        vcf_1kg = rules.NormaliseVcf.output
    output: "example/data/ThousandGenomes/chrMT_1kg_Human610-Quadv1_B-b37.vcf.gz"
    shell:
        "bcftools view -R {input.MTSnps} {input.vcf_1kg} -Oz -o {output}; "
        "bcftools index {output}"

rule ExampleSamplePanel:
    input:
        inVCF = "example/data/ThousandGenomes/chrMT_1kg_Human610-Quadv1_B-b37.vcf.gz",
        inSamples = "example/ExampleSamplePanel.txt"
    output:
        "example/SamplePanel/ExampleSamplePanel.vcf.gz"
    shell:
        "bcftools view -S {input.inSamples} -Oz -o {output} {input.inVCF}; "
        "bcftools index {output}"

rule vcf2Plink:
    input:
        vcf = "example/SamplePanel/ExampleSamplePanel.vcf.gz",
    output:
        expand("example/SamplePanel/ExampleSamplePanel.{ext}", ext = ['ped', 'map'])
    params:
        out = "example/SamplePanel/ExampleSamplePanel"
    shell:
        'plink --vcf {input.vcf} --recode --double-id --out {params.out}'

rule vcf2BPlink:
    input:
        vcf = "example/SamplePanel/ExampleSamplePanel.vcf.gz",
    output:
        expand("example/SamplePanel/ExampleSamplePanel.{ext}", ext = ['bed', 'bim', 'fam'])
    params:
        out = "example/SamplePanel/ExampleSamplePanel"
    shell:
        'plink --vcf {input.vcf} --make-bed --double-id --out {params.out}'
