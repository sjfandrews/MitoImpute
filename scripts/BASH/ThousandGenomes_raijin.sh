#!/bin/bash
#PBS -P te53
#PBS -q normalbw
#PBS -l walltime=12:00:00
#PBS -l mem=32GB
#PBS -l ncpus=1
#PBS -N impute_SNPchip_1kGP
#PBS -m be
#PBS -M u5015730@anu.edu.au
#PBS -j oe
#PBS -o /g/data1a/te53/MitoImpute/logs/

# LOAD THE MODULE
module unload intel-fc intel-cc
module load intel-fc/16.0.3.210
module load intel-cc/16.0.3.210
module load Rpackages/3.4.3
module load bcftools/1.8
module load plink/1.9
module load impute2/2.3.2
echo
echo "LOADED R v3.4.3"
echo "LOADED bcftools v1.8"
echo "LOADED plink v1.9"
echo "LOADED IMPUTE2 v2.3.2"

# EXTRACT PLATFORM SNPs
echo
echo "EXTRACTING PLATFORM SNPs"
#MtPlatforms=BDCHP-1X10-HUMANHAP240S_11216501_A-b37
MTSnps=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${MtPlatforms}_MT_snps.txt
vcf_1kg=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/chrMT_1kg_norm_decomposed_firstAlt.vcf.gz
vcf=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/chrMT_1kg_${MtPlatforms}.vcf.gz

bcftools view -R ${MTSnps} ${vcf_1kg} -Oz -o ${vcf}
bcftools index ${vcf}

# GENERATE GEN SAMPLE
echo
echo "GENERATING GEN SAMPLE"
#vcf=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/${MtPlatforms}/chrMT_1kg_${MtPlatforms}.vcf.gz
sex=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/SampleList1kg_sex.txt 
out=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/chrMT_1kg_${MtPlatforms}

bcftools convert --gensample ${out} ${vcf} --sex ${sex}
Rscript ~/GitCode/MitoImpute/scripts/R/FixSamplesFile_raijin.R ${out}.samples

# GENERATE PLINK FILES
echo
echo "GENERATING PLINK FILES"
out=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/chrMT_1kg_${MtPlatforms}

plink1.9 --vcf ${vcf} --recode --double-id --keep-allele-order --out ${out}

# RUN IMPUTE2
echo
echo "RUNNING IMPUTE2 ON ${MtPlatforms}"
m=~/GitCode/MitoImpute/DerivedData/ReferencePanel_v2/Reference_panel_v2_MtMap.txt 
h=~/GitCode/MitoImpute/DerivedData/ReferencePanel_v2/Reference_Panel_v2.hap.gz
l=~/GitCode/MitoImpute/DerivedData/ReferencePanel_v2/Reference_Panel_v2.legend.gz
g=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/chrMT_1kg_${MtPlatforms}.gen.gz
s=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/chrMT_1kg_${MtPlatforms}.samples
out=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/chrMT_1kg_${MtPlatforms}_imputed

if [ -f ${out} ]
then
	echo "$ {out} FOUND! ... PASSING"
else
	echo "$ {out} NOT FOUND! ... RUNNING IMPUTE2"
	echo impute2 -chrX -m ${m} -h ${h} -l ${l} -g ${g} -sample_g ${s} -int 1 16569 -Ne 20000 -o ${out}
fi

# FIX CHROMOSOME NAMES
echo
echo "FIXING CHROMOSOME NAMES"
InFile=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/${MtPlatforms}/chrMT_1kg_${MtPlatforms}_imputed
OutFile=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/${MtPlatforms}/chrMT_1kg_${MtPlatforms}_imputed_ChromFixed