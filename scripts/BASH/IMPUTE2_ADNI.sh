#!/bin/bash
#PBS -P te53
#PBS -q express
#PBS -l walltime=00:01:00
#PBS -l mem=1GB
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
echo ${REF_FILE}

#impute2 \
#-chrX \
#-m McInerney_Master_Alignment_Nov30_2017_ambig2missing_highQual_filtered map -h McInerney_Master_Alignment_Nov30_2017_ambig2missing_highQual_filtered hap.gz -l McInerney_Master_Alignment_Nov30_2017_ambig2missing_highQual_filtered legend.gz -g ADNI_samples gen.gz -sample_g ADNI_samples samples -int 1 16569 -Ne 20000 -o ADNI_samples_IMPUTED
