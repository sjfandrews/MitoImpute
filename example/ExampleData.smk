# snakemake -s example/ExampleData.smk
# snakemake -s example/ExampleData.smk --dag | dot -Tsvg > dag_Example.svg

##===========================##
##  Example Reference Panel
##===========================##

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
    input:
        in_samples = "example/ExampleReferenceSamples.txt",
        in_script = "scripts/R/assign_sex_label.R"
    output:
        outfile = "example/ReferencePanel/RefSampleList_sex.txt",
    shell:
        'Rscript {input.in_script} {input.in_samples} {output.outfile}'

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
        "example/ReferencePanel/ReferencePanel.vcf.gz",
        "example/ReferencePanel/RefSampleList_sex.txt"
    output:
        expand("example/ReferencePanel/ReferencePanel.{ext}", ext = ['hap.gz', 'legend.gz'])
    params:
        in_vcf = "example/ReferencePanel/ReferencePanel.vcf.gz",
        in_sex = "example/ReferencePanel/RefSampleList_sex.txt",
        out = "example/ReferencePanel/ReferencePanel"
    shell:
        'bcftools convert --haplegendsample {params.out} {params.in_vcf} --sex {params.in_sex}'

## 8. Generate .ped and .map files
rule Plink:
    input:
        "example/ReferencePanel/ReferencePanel.vcf.gz",
    output:
        expand("example/ReferencePanel/ReferencePanel.{ext}", ext = ['ped', 'map'])
    params:
        in_vcf = "example/ReferencePanel/ReferencePanel.vcf.gz",
        out = "example/ReferencePanel/ReferencePanel"
    shell:
        'plink --vcf {params.in_vcf} --recode --double-id --out {params.out}'

## 9. Generate .gen and .sample files
rule GenSample:
    input:
        "example/ReferencePanel/ReferencePanel.vcf.gz",
        "example/ReferencePanel/RefSampleList_sex.txt"
    output:
        expand("example/ReferencePanel/ReferencePanel.{ext}", ext = ['gen.gz', 'samples'])
    params:
        in_vcf = "example/ReferencePanel/ReferencePanel.vcf.gz",
        in_sex = "example/ReferencePanel/RefSampleList_sex.txt",
        out = "example/ReferencePanel/ReferencePanel"
    shell:
        'bcftools convert --gensample {params.out} {params.in_vcf} --sex {params.in_sex}'

## 10. Construct .map file for IMPUTE2
rule MakeMapFile:
    input:
        "scripts/R/mt_recombination_map.R",
        "example/ReferencePanel/ReferencePanel.vcf.gz"
    output:
        "example/ReferencePanel/MtMap.txt",
        "example/ReferencePanel/MtStrand.txt"
    params:
        in_vcf = "example/ReferencePanel/ReferencePanel.vcf.gz",
        in_script = "scripts/R/mt_recombination_map.R",
        out_map = "example/ReferencePanel/MtMap.txt",
        out_strand = "example/ReferencePanel/MtStrand.txt"
    shell:
        'Rscript {params.in_script} {params.in_vcf} {params.out_map} {params.out_strand}'

##===========================##
##  Example Sample Panel
##===========================##

rule ExampleSamplePanel:
    input:
        inVCF = "DerivedData/ThousandGenomes/Human610-Quadv1_B-b37/chrMT_1kg_Human610-Quadv1_B-b37.vcf.gz",
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
