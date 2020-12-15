# MitoImpute
**MitoImpute** is a snakemake pipeline for the imputation of mitochondrial genomes using Impute2 Chromosome X protocol. The steps in the pipline include:
1. Change sex of all samples to male (as males are haploid for the X chromsome)
2. Extract mtSNPs from Bplink (.bed/.bim/.fam)
3. Check reference alignment (hg19, Yoruba, or b37, rCRS) of mtSNPs - converts YRI to rCRS.
4. Converts Bplink files to:
   - oxford format (.gen/.sample)
   - plink format (.map/.ped)
5. Runs the chromsome X Impute2 imputation protocol. This step uses a custom mitochondrial reference panel constructed using the MitoImputePrep pipeline - see below.
6. Fixes chromosome label on the Impute2 output
7. Converts the Imputed files to:
   - Bplink format
   - Plink format
   - vcf format
8. Generates plots for Info Score and assigned Haplogroups

<img align="center" src=images/rulegraph.svg alt="DAG">

## Getting Started
### Installation

Be sure to download and install the latest versions of the following software packages:
1. [Miniconda Python 3](https://conda.io/en/latest/miniconda.html)
2. [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

MitoImpute can be installed on a git-enabled machine by typeing:

```bash
git clone https://github.com/sjfandrews/MitoImpute
```

## Usage Overview
### MitoImpute
To impute mitochondrial SNPs in a study dataset, run the following code:

```bash
snakemake -j --use-conda
```

Options for the snakemake file are set in the corresponding config file ```config/config.yaml``` file. The avaliable options are:

```bash
SAMPLE: 'name of input binary plink file'
DATAIN: 'path/to/input/directory'
DATAOUT: 'path/to/output/directory'
REFDATA: 'path/to/reference/panel'
INFOCUT: Info score threshold
ITER: Total number of MCMC iterations to perform, including burn-in.
BURNIN: Number of MCMC iteractions to discard as burn-in
```

The default options are for the example dataset.

#### Reference panel
A custom reference panel for imputation can be found in the ```resources/ReferencePanel/``` directory. The key files consist of:
1. -h: A file of known haplotypes ```(ReferencePanel.hap.gz)```.
2. -l: Legend file(s) with information about the SNPs in the -h file ```(ReferencePanel.legend.gz)```
3. -m: A fine-scale recombination map for the region to be analyzed ```(MtMap.txt)```

setting REFDATA in the ```config/config.yaml``` file to ```resources/ReferencePanel``` will automaticlay call these files.

#### Parallelization
Snakemake handles parallelization of jobs using [wildcards](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#wildcards). Defining a list of sample names in the config.yaml file and specifing the [number of avaliable cores](https://snakemake.readthedocs.io/en/stable/executable.html#useful-command-line-arguments) in the command line will result in snakemake submitting jobs in parallel.

config file with multiple samples
```
SAMPLE: ['ExampleSamplePanel', 'ExampleSamplePanel_2', 'ExampleSamplePanel_3']
DATAIN: 'example/SamplePanel'
DATAOUT: "example/imputed"
REFDATA: "example/ReferencePanel"
INFOCUT: 0
ITER: 2
BURNIN: 1
KHAP: 1000
```

the corresponding comand line argument
```
snakemake -j 3 --use-conda
```

#### Cluster Execution
For excuting this pipline on a cluster computing environment, refer to [SnakeMake's readme](https://snakemake.readthedocs.io/en/stable/executable.html#cluster-execution).

## Dependencies

MitoImpute uses snakemake conda environments to handle package dependencies - if you dont want to use conda envs (which is not recomended) ensure the following software and packages are installed.

1. [R](https://cran.r-project.org/)
2. [PLINK](https://www.cog-genomics.org/plink2)
3. [Impute2](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#download)

Plink and Impute2 executibles should be located within the the /usr/local/bin/ directory. The following code can be used to move the executibles: ```cp </path/to/executible> /usr/local/bin/```

The following R packages are also required:
1. [tidyverse](https://www.tidyverse.org/packages/)
2. [Hi-MC](https://github.com/vserch/himc)
3. [ggfittext](https://cran.r-project.org/web/packages/ggfittext/index.html)

Note that the development versions of Hi-MC (required for mitochondrial haplogroup assignment) are required.

```r
## Install tidyverse, rmarkdown, and devtools
install.packages(c("tidyverse", "ggfittext", "devtools"))

## Install HiMC
devtools::install_github(c("vserch/himc/HiMC"))
```
