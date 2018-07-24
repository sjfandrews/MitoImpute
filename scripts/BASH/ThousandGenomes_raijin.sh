#!/bin/bash

# LOAD MODULE FILES
module load bcftools/1.8
module load plink/1.9

# EXTRACT PLATFORM SNPs
echo "EXTRACTING PLATFORM SNPs"
MtPlatforms=BDCHP-1X10-HUMANHAP240S_11216501_A-b37
MTSnps=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${MtPlatforms}_MT_snps.txt
vcf_1kg=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/chrMT_1kg_norm_decomposed_firstAlt.vcf.gz
vcf=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/chrMT_1kg_${MtPlatforms}.vcf.gz

bcftools view -R ${MTSnps} ${vcf_1kg} -Oz -o ${vcf}
bcftools index ${vcf}

# GENERATE GEN SAMPLE
echo "GENERATING GEN SAMPLE"
#vcf=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/${MtPlatforms}/chrMT_1kg_${MtPlatforms}.vcf.gz
sex=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/SampleList1kg_sex.txt 
out=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/chrMT_1kg_${MtPlatforms}

bcftools convert --gensample ${out} ${vcf} --sex ${sex}
Rscript ~/GitCode/MitoImpute/scripts/R/FixSamplesFile.R ${out}.samples

# GENERATE PLINK FILES
echo "GENERATING PLINK FILES"
out=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/chrMT_1kg_${MtPlatforms}

plink1.9 --vcf ${vcf} --recode --double-id --keep-allele-order --out ${out}

# RUN IMPUTE2
echo "RUNNING IMPUTE2 ON ${MtPlatforms}"
m=~/GitCode/MitoImpute/DerivedData/ReferencePanel_v2/Reference_panel_v2_MtMap.txt 
h=~/GitCode/MitoImpute/DerivedData/ReferencePanel_v2/Reference_Panel_v2.hap.gz
l=~/GitCode/MitoImpute/DerivedData/ReferencePanel_v2/Reference_Panel_v2.legend.gz
g=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/chrMT_1kg_${MtPlatforms}.gen.gz
s=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/chrMT_1kg_${MtPlatforms}.samples
out=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/chrMT_1kg_${MtPlatforms}_imputed

echo impute2 -chrX -m ${m} -h ${h} -l ${l} -g ${g} -sample_g ${s} -int 1 16569 -Ne 20000 -o ${out}

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