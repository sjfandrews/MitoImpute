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