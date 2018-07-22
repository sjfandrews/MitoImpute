#!/bin/bash
#PBS -P te53
#PBS -q express
#PBS -l walltime=01:00:00
#PBS -l mem=32GB
#PBS -l ncpus=1
#PBS -N mtRefPanel
#PBS -m be
#PBS -M u5015730@anu.edu.au
#PBS -j oe
#PBS -o /g/data1a/te53/MitoImpute/logs/

# LOAD THE MODULE
module unload intel-fc intel-cc
module load intel-fc/16.0.3.210
module load intel-cc/16.0.3.210
module load Rpackages/3.4.3
module load R/3.4.3
module load python3/3.6.2
pyvenv /g/data1a/te53/MitoImpute/Py3virt
source /g/data1a/te53/MitoImpute/Py3virt/bin/activate
module load vt
module load bcftools/1.8
module avail impute2

# RUN THE COMMANDS
cd ~/GitCode/MitoImpute
snakemake -s ReferencePanel_v2.smk