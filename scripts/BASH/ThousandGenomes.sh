#!/bin/bash

# EXTRACT PLATFORM SNPs
echo "EXTRACTING PLATFORM SNPs"
MtPlatforms=BDCHP-1X10-HUMANHAP240S_11216501_A-b37
MTSnps=~/GitCode/MitoImpute/data/platforms/${MtPlatforms}/${MtPlatforms}_MT_snps.txt
vcf_1kg=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/chrMT_1kg_norm_decomposed_firstAlt.vcf.gz
vcf=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/${MtPlatforms}/chrMT_1kg_${MtPlatforms}.vcf.gz

bcftools view -R ${MTSnps} ${vcf_1kg} -Oz -o ${vcf}
bcftools index ${vcf}

# GENERATE GEN SAMPLE
echo "GENERATING GEN SAMPLE"
vcf=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/${MtPlatforms}/chrMT_1kg_${MtPlatforms}.vcf.gz
sex=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/SampleList1kg_sex.txt
out=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/${MtPlatforms}/chrMT_1kg_${MtPlatforms}

bcftools convert --gensample ${out} ${vcf} --sex ${sex}
Rscript ~/GitCode/MitoImpute/scripts/R/FixSamplesFile.R ${out}.samples

# GENERATE PLINK FILES
echo "GENERATING PLINK FILES"
out=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/${MtPlatforms}/chrMT_1kg_${MtPlatforms}

plink1.9 --vcf ${vcf} --recode --double-id --keep-allele-order --out ${out}

# RUN IMPUTE2
echo "RUNNING IMPUTE2 ON ${MtPlatforms}"
m=/Volumes/MHS/MitoImpute/data/REF_PANEL/Reference_panel_v2_MtMap.txt 
h=/Volumes/MHS/MitoImpute/data/OXFORD/Reference_Panel_v2.hap.gz
l=/Volumes/MHS/MitoImpute/data/OXFORD/Reference_Panel_v2.legend.gz
g=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/${MtPlatforms}/chrMT_1kg_${MtPlatforms}.gen.gz
s=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/${MtPlatforms}/chrMT_1kg_${MtPlatforms}.samples
out=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/${MtPlatforms}/chrMT_1kg_${MtPlatforms}_imputed

impute2 -chrX -m ${m} -h ${h} -l ${l} -g ${g} -sample_g ${s} -int 1 16569 -Ne 20000 -o ${out}

#rule Impute2:
#    input:
#        m = expand('{RefData}/MtMap.txt', RefData=REFDATA),
#        h = expand('{RefData}/ReferencePanel.hap.gz', RefData=REFDATA),
#        l = expand('{RefData}/ReferencePanel.legend.gz', RefData=REFDATA),
#        g = "DerivedData/ThousandGenomes/{MtPlatforms}/chrMT_1kg_{MtPlatforms}.gen.gz",
#        sample = "DerivedData/ThousandGenomes/{MtPlatforms}/chrMT_1kg_{MtPlatforms}.samples",
#    output:
#        'DerivedData/ThousandGenomes/{MtPlatforms}/chrMT_1kg_{MtPlatforms}_imputed',
#        'DerivedData/ThousandGenomes/{MtPlatforms}/chrMT_1kg_{MtPlatforms}_imputed_samples'
#    params:
#        out = 'DerivedData/ThousandGenomes/{MtPlatforms}/chrMT_1kg_{MtPlatforms}_imputed'
#    shell:
#        'impute2 -chrX -m {input.m} -h {input.h} -l {input.l} -g {input.g} \
#        -sample_g {input.sample} -int 1 16569 -Ne 20000 -o {params.out}'
#REFDATA = "example/ReferencePanel"