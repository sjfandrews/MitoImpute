# MitoImpute
**mtImpute.smk** is a snakemake pipeline for the imputation of mitochondrial genomes using Impute2 Chromosome X protocol. The steps in the pipline include:
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
8. Generates a html rmarkdown report

<img align="center" src=dag_mtImpute.svg alt="DAG">

## Getting Started
### Installation
Be sure to download and install the latest versions of the following software packages:
1. [Python 3](https://www.python.org/downloads/)
2. [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
3. [R](https://cran.r-project.org/)
4. [Rstudio](https://www.rstudio.com/products/rstudio/download/)
5. [PLINK](https://www.cog-genomics.org/plink2)
6. [Impute2](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#download)

Rstudio is recomended for running the MitoImputeShiny app. Plink and Impute2 executibles should be located within the the /usr/local/bin/ directory. The following code can be used to move the executibles: ```cp </path/to/executible> /usr/local/bin/```

The following R packages are also required:
1. [tidyverse](https://www.tidyverse.org/packages/)
2. [rmarkdown](https://cran.r-project.org/web/packages/rmarkdown/index.html)
4. [shiny](https://cran.r-project.org/web/packages/shiny/index.html)
5. [Hi-MC](https://github.com/vserch/himc)
6. [ggforce](https://github.com/thomasp85/ggforce)

Note that the development versions of ggforce (required for plotting alluvial diagrams) and Hi-MC (required for mitochondrial haplogroup assignment) are required. These packages can be isntalled directly from github using devtools.

```r
## Install tidyverse, rmarkdown, and devtools
install.packages(c("tidyverse", "shiny", "rmarkdown", "devtools"))

## Install HiMC and ggforce
devtools::install_github(c("vserch/himc/HiMC", 'thomasp85/ggforce'))
```

Once all the prerequiste software is isntalled, MitoImpute can be installed on a git-enabled machine by typeing:

```bash
git clone https://github.com/sjfandrews/MitoImpute
```

## Usage Overview
### MitoImpute
To impute mitochondrial SNPs in a study dataset, run the following code:

```bash
cd MitoImpute
snakemake -s mtImpute.smk
```

If the rmarkdown report is not required, run the following code:

```bash
cd MitoImpute
snakemake -s mtImpute.smk --until oxford2vcf oxford2ped oxford2bed bplink2plink
```

Options for the snakemake file are set in the corresponding config file ```mtImpute_config.yaml``` file. The avaliable options are:

```bash
SAMPLE: 'name of input binary plink file'
DATAIN: 'path/to/input/directory'
DATAOUT: 'path/to/output/directory'
REFDATA: 'path/to/reference/panel'
INFOCUT: Info score threshold
```

The default options are for the example dataset.

#### Reference panel
A custom reference panel for imputation can be found in the ```MitoImpute/ReferencePanel/``` directory. The key files consist of:
1. -h: A file of known haplotypes ```(ReferencePanel.hap.gz)```.
2. -l: Legend file(s) with information about the SNPs in the -h file ```(ReferencePanel.legend.gz)```
3. -m: A fine-scale recombination map for the region to be analyzed ```(MtMap.txt)```

setting REFDATA in the ```mtImpute_config.yaml``` file to ```ReferencePanel``` will automaticlay call these files.

Further details on the pipeline for constructing the reference panel can be found in the [MitoImputePrep](https://github.com/sjfandrews/MitoImputePrep) repo.

### MitoImpute Shiny
MitoImputeShiny is a R shiny application for displaying the results of the 1000 Genomes MitoImpute validation pipeline. As detailed in the [MitoImputePrep](https://github.com/sjfandrews/MitoImputePrep) pipline, [Strand Files](http://www.well.ox.ac.uk/~wrayner/strand/), which provide the strand orientation and position of variants of the most common genotyping platform on build 37, were used to obtain a list of mtSNPs for each platform. The platform specific mtSNPs are then extracted from the Thousand Genome mitochondrial Sequences and Using these 'Typed SNPs', additional mitochondrial SNPs were then imputed using the reference panel. The Shiny application displays the concordance of haplogroup assignments made using the full mitochondrial sequence data vs using only the Typed mtSNPs or Typed + Imputed mtSNPs.

To run MitoImputeShiny, in Rstudio run the following code:

```r
library(shiny)
setwd("~/path/to/MitoImpute/Shiny")
runApp()
```
