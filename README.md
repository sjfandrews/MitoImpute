[![DOI](https://zenodo.org/badge/143345792.svg)](https://zenodo.org/badge/latestdoi/143345792)

# Reference Alignments

The reference alignments included in this repository are:
*	`resources/alignments/McInerney_Master_Alignment_July18_2018.fasta.gz`
*	`resources/alignments/hsapiensCRS7k.fasta.gz`

`McInerney_Master_Alignment_July18_2018.fasta.gz` is the novel reference alignment constructed in 2018 from the sequences downloaded on the 18th of July, 2018. It contains 44,299 aligned complete mitochondrial DNA sequences. These sequences are all 16,569 DNA nucleotide states long (to match the numbering conventions of the revised Cambridge Reference Sequence - Andrews et al., 1999). From this alignment the Reference Panels were filtered down to 36,960 sequences and filtered to thresholds detailed in [McInerney et al. (2020)](https://www.biorxiv.org/content/10.1101/649293v3).

`hsapiensCRS7k.fasta.gz` is the previous reference alignment constructed in 2011 by Dr's Simon Easteal and Lars Jermiin. It contains 7,747 aligned complete mitochondrial DNA sequences. These sequences are all 16,569 DNA nucleotide states long (to match the numbering conventions of the revised Cambridge Reference Sequence - Andrews et al., 1999). This curated alignment was used to align the sequences downloaded on the 18th of July 2018. Novel sequences were aligned in batches of 2,500 sequences. Any gaps forced into `hsapiensCRS7k.fasta.gz` were removed.

# Reference Panels

The reference panels included in this repository are:
*	`ReferencePanel_v1_0.01`
*	`ReferencePanel_v1_0.005`
*	`ReferencePanel_v1_0.001`

Each of these corresponds to a filtering of sites to a minor allele frequency of 1%, 0.5%, and 0.01%, respectively. All reference panels contain variant information for the 36,960 sequences from the `McInerney_Master_Alignment_July18_2018.fasta.gz` reference alignment.

All references panels can be found in the directory: `MitoImpute/resources/` in their own specific subdirectory. Each reference panel is stored in VCF (\*.vcf.gz) and oxford (\*.gen.gz, .*hap.gz, *.legend.gz, *.samples) formats. These files will be necessary for using a reference panel for genotype imputation.

# MitoImpute
**MitoImpute** is a snakemake pipeline for the imputation of mitochondrial genomes using Impute2 Chromosome X framework. The steps in the pipline include:
1. Change sex of all samples to male (as males are haploid for the X chromsome)
2. Extract mtSNPs from Bplink (.bed/.bim/.fam)
3. Check reference alignment (hg19, Yoruba, or b37, rCRS) of mtSNPs - converts YRI to rCRS.
4. Converts Bplink files to:
   - oxford format (.gen/.sample)
   - plink format (.map/.ped)
5. Runs the chromsome X Impute2 imputation protocol. This step uses a custom mitochondrial genome reference panel.
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
REFAF: [0.01, 0.005, 0.001, 'example']
INFOCUT: Info score threshold
ITER: Total number of MCMC iterations to perform, including burn-in.
BURNIN: Number of MCMC iteractions to discard as burn-in
```

The default options are for the example dataset.

#### Reference panel
A custom reference panels for imputation can be found in the ```resources/``` directory. The key files consist of:
1. -h: A file of known haplotypes ```(ReferencePanel.hap.gz)```.
2. -l: Legend file(s) with information about the SNPs in the -h file ```(ReferencePanel.legend.gz)```
3. -m: A fine-scale recombination map for the region to be analyzed ```(MtMap.txt)```

setting REFAF in the `config/config.yaml` file to the desired MAF will automaticlay call these files.

#### Parallelization
Snakemake handles parallelization of jobs using [wildcards](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#wildcards). Defining a list of sample names in the config.yaml file and specifing the [number of avaliable cores](https://snakemake.readthedocs.io/en/stable/executable.html#useful-command-line-arguments) in the command line will result in snakemake submitting jobs in parallel.

config file with multiple samples
```
SAMPLE: ['example', 'example_2', 'example_3']
DATAIN: 'resources/example'
DATAOUT: "results/example"
REFAF: "example"
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


## Citation
McInerney, T. W. et al. (2020) A globally diverse reference alignment and panel for imputation of mitochondrial DNA variants. [BioRxiv](https://www.biorxiv.org/content/10.1101/649293v3).
