#!/bin/bash
#PBS -P te53
#PBS -q normal
#PBS -l walltime=48:00:00
#PBS -l mem=32GB
#PBS -l ncpus=1
#PBS -N IMPUTE2_ADNI
#PBS -m be
#PBS -M u5015730@anu.edu.au
#PBS -j oe
#PBS -o /g/data1a/te53/MitoImpute/logs/

# LOAD THE MODULE
module load impute2/2.3.2

# RUN THE COMMANDS
REF_DIR="/g/data1a/te53/MitoImpute/data/REF_PANEL/"
REF_FILE="${REF_DIR}Reference_panel_v2"
ADNI_FILE="${REF_DIR}ADNI_samples"

echo "COMMAND RUN:"
echo impute2 -chrX -m ${REF_FILE}.map -h ${REF_FILE}.hap.gz -l ${REF_FILE}.legend.gz -g ${ADNI_FILE}.gen.gz -sample_g ${ADNI_FILE}.samples -int 1 16569 -Ne 20000 -o ${ADNI_FILE}_IMPUTED
impute2 -chrX -m ${REF_FILE}.map -h ${REF_FILE}.hap.gz -l ${REF_FILE}.legend.gz -g ${ADNI_FILE}.gen.gz -sample_g ${ADNI_FILE}.samples -int 1 16569 -Ne 20000 -o ${ADNI_FILE}_IMPUTED
