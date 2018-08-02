#!/bin/bash
#PBS -P te53
#PBS -q normal
#PBS -l walltime=48:00:00
#PBS -l mem=32GB
#PBS -l ncpus=1
#PBS -N impute_SNPchip_1kGP
#PBS -m e
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

# SPECIFY REFERENCE PANEL
REFpanel="ReferencePanel_v3"
echo
echo "REFERENCE PANEL: ${REFpanel}"

# CREATE DIRECTORY
if [ -d  /g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/ ]
then
	echo
	echo "/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/ EXISTS ... PASSING"
else
	echo
	echo "/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/ NOT FOUND ... CREATING DIRECTORY"
	mkdir /g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/
fi

# EXTRACT PLATFORM SNPs
echo
echo "EXTRACTING PLATFORM SNPs"
#MtPlatforms=BDCHP-1X10-HUMANHAP240S_11216501_A-b37
MTSnps=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${MtPlatforms}_MT_snps.txt
vcf_1kg=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/chrMT_1kg_norm_decomposed_firstAlt.vcf.gz
vcf=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}.vcf.gz

bcftools view -R ${MTSnps} ${vcf_1kg} -Oz -o ${vcf}
bcftools index ${vcf}

# GENERATE GEN SAMPLE
echo
echo "GENERATING GEN SAMPLE"
#vcf=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/${MtPlatforms}/chrMT_1kg_${MtPlatforms}.vcf.gz
sex=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/SampleList1kg_sex.txt 
out=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}

bcftools convert --gensample ${out} ${vcf} --sex ${sex}
Rscript ~/GitCode/MitoImpute/scripts/R/FixSamplesFile_raijin.R ${out}.samples

# GENERATE PLINK FILES
echo
echo "GENERATING PLINK FILES"
out=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}

plink1.9 --vcf ${vcf} --recode --double-id --keep-allele-order --out ${out}

# RUN IMPUTE2
echo
echo "RUNNING IMPUTE2 ON ${MtPlatforms}"
m=~/GitCode/MitoImpute/DerivedData/${REFpanel}/${REFpanel}_MtMap.txt 
h=~/GitCode/MitoImpute/DerivedData/${REFpanel}/${REFpanel}.hap.gz
l=~/GitCode/MitoImpute/DerivedData/${REFpanel}/${REFpanel}.legend.gz
g=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}.gen.gz
s=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}.samples
out=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed

if [ -f ${out} ]
then
	echo "${out} FOUND! ... PASSING"
else
	echo "${out} NOT FOUND! ... RUNNING IMPUTE2"
	impute2 -chrX -m ${m} -h ${h} -l ${l} -g ${g} -sample_g ${s} -int 1 16569 -Ne 20000 -o ${out}
fi

# FIX CHROMOSOME NAMES
echo
echo "FIXING CHROMOSOME NAMES"
InFile=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed
OutFile=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed_ChromFixed
awk '{{$1 = "26"; print}}' ${InFile} > ${OutFile}

# CONVERT OXFORD TO PEDIGREE
echo
echo "CONVERTING OXFORD TO PEDIGREE"
gen=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed_ChromFixed
sam=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed_samples
out=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed

plink1.9 --gen ${gen} --sample ${sam} --hard-call-threshold 0.49 --keep-allele-order --output-chr 26 --recode --out ${out}

# CONVERT OXFORD TO VCF
echo
echo "CONVERTING OXFORD TO VCF"
gen=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed_ChromFixed
sam=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed_samples
out=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed

plink1.9 --gen ${gen} --sample ${sam} --hard-call-threshold 0.49 --keep-allele-order --output-chr 26 --recode vcf --out ${out}


# GENERATE QC REPORT
echo
echo "GENERATING QC REPORT"
s=~/GitCode/MitoImpute/scripts/R/MT_imputation_QC.Rmd
wgs_map=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/chrMT_1kg_norm_decomposed_firstAlt.map
wgs_ped=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/chrMT_1kg_norm_decomposed_firstAlt.ped
wgs_vcf=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/chrMT_1kg_norm_decomposed_firstAlt.vcf.gz
typ_map=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}.map
typ_ped=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}.ped
typ_vcf=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}.vcf.gz
imp_map=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed.map
imp_ped=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed.ped
imp_vcf=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed.vcf
imp_info=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed_info

output=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_mtImputed_QC.html

rwd = `pwd`/
output_dir=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/${MtPlatforms}/
info_cut='0'