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
#PBS -o /g/data1a/te53/haploco/logs/

# LOAD THE MODULE
initPy3virt_MitoImpute
initR
module load vt
module load bcftools/1.8

# RUN THE COMMANDS
cd ~/GitCode/MitoImpute
snakemake -s ReferencePanel_v2.smk