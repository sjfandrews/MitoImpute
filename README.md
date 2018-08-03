# MitoImpute
**mtImpute.smk** is a snakemake pipeline for the imputation of mitochondrial genomes using Impute2 Chromosome X protocol. The steps in the pipline include:
1. Change sex of all samples to male (as males are haploid for the X chromsome)
2. Converts Bplink (.bed/.bim/.fam) files to:
   - oxford format (.gen/.sample)
   - plink format (.map/.ped)
3. Runs the chromsome X Impute2 imputation protocol. This step uses a custom mitochondrial reference panel constructed using the MitoImputePrep pipeline - see below.
4. Fixes chromosome label on the Impute2 output
5. Converts the Imputed files to:
   - Bplink format
   - Plink format
   - vcf format
6. Generates a html rmarkdown report

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

The following Python modules are required:
1. pysam

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
SAMPLE: 'name of input binary plink file'
DATAIN: 'path/to/input/directory'
DATAOUT: 'path/to/output/directory'
REFDATA: 'path/to/reference/panel'
```

The default options are for the example dataset.

### Reference panel
A custom reference panel for imputation can be found in the ```MitoImpute/DerivedData/ReferencePanel/``` directory. The key files consist of:
1. -h: A file of known haplotypes ```(ReferencePanel.hap.gz)```.
2. -l: Legend file(s) with information about the SNPs in the -h file ```(ReferencePanel.legend.gz)```
3. -m: A fine-scale recombination map for the region to be analyzed ```(MtMap.txt)```

setting REFDATA in the ```mtImpute_config.yaml``` file to ```path/to/MitoImpute/DerivedData/ReferencePanel``` will automaticlay call these files.

To construct the reference panel, 44,299 sequences were downloaded from NCBI using the search term:

```(016500[SLEN]:016600[SLEN]) AND Homo[Organism] AND mitochondrion[FILT] AND complete genome NOT (Homo sp. Altai OR Denisova hominin OR neanderthalensis OR heidelbergensis OR consensus OR ancient human remains OR shotgun)```

The 44,299 sequences were aligned to Easteal, Jermiin, and Ott mitochondrial master alignment (~8,000 manually aligned sequneces...) using Geneious 10.2.6.
