library(tidyverse)

args = commandArgs(trailingOnly = TRUE) # Set arguments from the command line

infile = args[1]

outfile = paste0(infile, "_fixed")

imp.bim <- read_delim(infile, delim = " ", col_names = F)
imp.bim$X1 <- '26'
write_delim(imp.bim, outfile, col_names = F)
