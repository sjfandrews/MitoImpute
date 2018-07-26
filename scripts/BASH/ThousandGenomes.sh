#!/bin/bash

MtPlatforms=BDCHP-1X10-HUMANHAP240S_11216501_A-b37
# EXTRACT PLATFORM SNPs
echo
echo "EXTRACTING PLATFORM SNPs"
MTSnps=~/GitCode/MitoImpute/data/platforms/${MtPlatforms}/${MtPlatforms}_MT_snps.txt
vcf_1kg=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/chrMT_1kg_norm_decomposed_firstAlt.vcf.gz
vcf=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/${MtPlatforms}/chrMT_1kg_${MtPlatforms}.vcf.gz

bcftools view -R ${MTSnps} ${vcf_1kg} -Oz -o ${vcf}
bcftools index ${vcf}

# GENERATE GEN SAMPLE
echo
echo "GENERATING GEN SAMPLE"
vcf=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/${MtPlatforms}/chrMT_1kg_${MtPlatforms}.vcf.gz
sex=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/SampleList1kg_sex.txt
out=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/${MtPlatforms}/chrMT_1kg_${MtPlatforms}

bcftools convert --gensample ${out} ${vcf} --sex ${sex}
Rscript ~/GitCode/MitoImpute/scripts/R/FixSamplesFile.R ${out}.samples

# GENERATE PLINK FILES
echo
echo "GENERATING PLINK FILES"
out=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/${MtPlatforms}/chrMT_1kg_${MtPlatforms}

plink1.9 --vcf ${vcf} --recode --double-id --keep-allele-order --out ${out}

# RUN IMPUTE2
echo
echo "RUNNING IMPUTE2 ON ${MtPlatforms}"
m=/Volumes/MHS/MitoImpute/data/REF_PANEL/Reference_panel_v2_MtMap.txt 
h=/Volumes/MHS/MitoImpute/data/OXFORD/Reference_Panel_v2.hap.gz
l=/Volumes/MHS/MitoImpute/data/OXFORD/Reference_Panel_v2.legend.gz
g=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/${MtPlatforms}/chrMT_1kg_${MtPlatforms}.gen.gz
s=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/${MtPlatforms}/chrMT_1kg_${MtPlatforms}.samples
out=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/${MtPlatforms}/chrMT_1kg_${MtPlatforms}_imputed

if [ -f ${out} ]
then
	echo "${out} FOUND! ... PASSING"
else
	echo "${out} NOT FOUND! ... RUNNING IMPUTE2"
	echo impute2 -chrX -m ${m} -h ${h} -l ${l} -g ${g} -sample_g ${s} -int 1 16569 -Ne 20000 -o ${out}
fi

# FIX CHROMOSOME NAMES
echo
echo "FIXING CHROMOSOME NAMES"
InFile=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/${MtPlatforms}/chrMT_1kg_${MtPlatforms}_imputed
OutFile=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/${MtPlatforms}/chrMT_1kg_${MtPlatforms}_imputed_ChromFixed
awk '{{$1 = "26"; print}}' ${InFile} > ${OutFile}

# CONVERT OXFORD TO PEDIGREE
echo
echo "CONVERTING OXFORD TO PEDIGREE"
gen=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/${MtPlatforms}/chrMT_1kg_${MtPlatforms}_imputed_ChromFixed
sam=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/${MtPlatforms}/chrMT_1kg_${MtPlatforms}_imputed_samples
out=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/${MtPlatforms}/chrMT_1kg_${MtPlatforms}_imputed

plink1.9 --gen ${gen} --sample ${sam} --hard-call-threshold 0.49 --keep-allele-order --output-chr 26 --recode --out ${out}

# CONVERT OXFORD TO VCF
echo
echo "CONVERTING OXFORD TO VCF"
gen=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/${MtPlatforms}/chrMT_1kg_${MtPlatforms}_imputed_ChromFixed
sam=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/${MtPlatforms}/chrMT_1kg_${MtPlatforms}_imputed_samples
out=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/${MtPlatforms}/chrMT_1kg_${MtPlatforms}_imputed

plink1.9 --gen ${gen} --sample ${sam} --hard-call-threshold 0.49 --keep-allele-order --output-chr 26 --recode vcf --out ${out}


# GENERATE QC REPORT
echo
echo "GENERATING QC REPORT"
s=~/GitCode/MitoImpute/scripts/R/MT_imputation_QC.Rmd
wgs_map=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/chrMT_1kg_norm_decomposed_firstAlt.map
wgs_ped=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/chrMT_1kg_norm_decomposed_firstAlt.ped
wgs_vcf=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/chrMT_1kg_norm_decomposed_firstAlt.vcf.gz
typ_map=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/${MtPlatforms}/chrMT_1kg_${MtPlatforms}.map
typ_ped=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/${MtPlatforms}/chrMT_1kg_${MtPlatforms}.ped
typ_vcf=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/${MtPlatforms}/chrMT_1kg_${MtPlatforms}.vcf.gz
imp_map=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/${MtPlatforms}/chrMT_1kg_${MtPlatforms}_imputed.map
imp_ped=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/${MtPlatforms}/chrMT_1kg_${MtPlatforms}_imputed.ped
imp_vcf=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/${MtPlatforms}/chrMT_1kg_${MtPlatforms}_imputed.vcf
imp_info=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/${MtPlatforms}/chrMT_1kg_${MtPlatforms}_imputed_info

output=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/${MtPlatforms}/chrMT_1kg_{MtPlatforms}_mtImputed_QC.html

rwd=`pwd`/
output_dir=~/GitCode/MitoImpute/DerivedData/ThousandGenomes/${MtPlatforms}/
info_cut='0'

#R -e 'rmarkdown::render(rmarkdown::render(${s}, output_file = ${output}, output_dir = ${output_dir}, params = list(rwd = ${rwd}, info.cut = ${info_cut}, wgs.map = ${wgs_map}, wgs.ped = ${wgs_ped}, wgs.vcf = ${wgs_vcf}, typ.map = "{input.typ_map}", typ.ped = ${typ_ped}, typ.vcf = ${typ_vcf}, imp.map = ${imp_map}, imp.ped = ${imp_ped}, imp.vcf= ${imp_vcf}, imp.info = ${input.imp_info}))' --slave
echo
echo "rmarkdown::render(${s}, output_file = ${output}, output_dir = ${output_dir}, params = list(rwd = ${rwd}, info.cut = ${info_cut}, wgs.map = ${wgs_map}, wgs.ped = ${wgs_ped}, wgs.vcf = ${wgs_vcf}"#, typ.map = {input.typ_map}, typ.ped = ${typ_ped}, typ.vcf = ${typ_vcf}, imp.map = ${imp_map}, imp.ped = ${imp_ped}, imp.vcf= ${imp_vcf}, imp.info = ${input.imp_info}))"
