# MitoImpute
A snakemake pipeline for the imputation of mitochondrial genomes using Impute2.

Currenlty there are four snakemake files:

1.  **ReferencePanel.smk:** Pipeline for processing a FASTA alignment of mitochondrial sequnces to be used as the reference panel in imputation. A pre-build reference panel will be provided so that is not needed to be ran.
2.  **PlatformStrandFiles.smk:** Script to downloand [files](http://www.well.ox.ac.uk/~wrayner/strand/) giving the strand orientation and position of variants of the most common genotyping chips on build 37. These files are used to create 'in silico' microarray chips of mitochondrial variants to evaluate imputation quality in Tousand Genomes.
3.  **ThousandGenomes.smk:** Pipline to download and process mitochondrial sequences from thousand genomes and creat 'in silico' microarrys based on the downloaded strand files.
4.  **mtImpute.smk:** Pipeline for imputing untyped mitochondrial SNPs in a study dataset from a custom mitochondrial referen panel (generated in ReferencePanel.smk)

## Getting Started
