# MitoImpute
A snakemake pipeline for the imputation of mitochondrial genomes using Impute2.

Currenlty there are four snakemake files:

1.  **ReferencePanel.smk:** Pipeline for processing a FASTA alignment of mitochondrial sequnces to be used as the reference panel in imputation. A pre-build reference panel will be provided so that is not needed to be ran.
2.  **PlatformStrandFiles.smk:** Script to downloand [files](http://www.well.ox.ac.uk/~wrayner/strand/) giving the strand orientation and position of variants of the most common genotyping chips on build 37. These files are used to create 'in silico' microarray chips of mitochondrial variants to evaluate imputation quality in Tousand Genomes.
3.  **ThousandGenomes.smk:** Pipline to download and process mitochondrial sequences from thousand genomes and creat 'in silico' microarrys based on the downloaded strand files.
4.  **mtImpute.smk:** Pipeline for imputing untyped mitochondrial SNPs in a study dataset from a custom mitochondrial referen panel (generated in ReferencePanel.smk)

The mtImpute.smk pipeline is the sole pipeline needed for imputation of mitochondrial SNPs.

## Getting Started
### Installation
Be sure to download and install the latest versions of the following software packages:
1. [Python 3](https://www.python.org/downloads/)
2. [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
3. [R](https://cran.r-project.org/)
4. [PLINK](https://www.cog-genomics.org/plink2)
5. [Impute2](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#download)

The following R packages are also required:
1. tidyverse
2. ggplot2
3. rmarkdown
4. [Hi-MC](https://github.com/vserch/himc)
5. [ggforce](https://github.com/thomasp85/ggforce)

Note that the development versions of ggforce (required for plotting alluvial diagrams) and Hi-MC (required for mitochondrial haplogroup assignment) are required. These packages can be isntalled directly from github using devtools (see their respective pages).

Once all the prerequiste software is isntalled, MitoImpute can be installed on a git-enabled machine by typeing:

```bash
git clone https://github.com/sjfandrews/MitoImpute
```

### Usage Overview
To impute mitochondrial SNPs in a study dataset, run the following code:

```bash
snakemake -s mtImpute.smk
```

Options for the snakemake file are set in the corresponding config file ```mtImpute_config.yaml``` file. The avaliable options are:

```bash
SAMPLE: 'name of binary plink file'
DATAIN: 'path/to/plink/file'
DATAOUT: 'path/to/output/file'
REFDATA: 'path/to/reference/panel'
```

The defualt options are for the example dataset.

### Reference panel
A custom reference panel for imputation can be found in the ```MitoImpute/DerivedData/ReferencePanel/``` directory. The key files consist of:
1. -h: A file of known haplotypes ```(ReferencePanel.hap.gz)```.
2. -l: Legend file(s) with information about the SNPs in the -h file ```(ReferencePanel.legend.gz)```
3. -m: A fine-scale recombination map for the region to be analyzed ```(MtMap.txt)```

setting REFDATA in the ```mtImpute_config.yaml``` file to ```path/to/MitoImpute/DerivedData/ReferencePanel``` will automaticlay call these files.
