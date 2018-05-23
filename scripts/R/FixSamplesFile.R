suppressPackageStartupMessages(library(tidyverse))

## Read in Arguments 
# 1: Current working directory 
args = commandArgs(trailingOnly=TRUE)
file.path <- args[1]

dat <- read_delim(file.path, delim = " ", col_names = T)
dat[1,4] <- "D"
write_delim(dat, file.path, delim = " ", col_names = T)