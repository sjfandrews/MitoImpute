#!/usr/bin/Rscript
##============================================================================##
##      R script to construct a file containing sex information of
##      Reference samples. .bgen and .GenSamples require sex information
##      and as we are using the X chromosmoe protocol for impute2,
##      we need to speicify all samples are Male due to haploid nature of
##      the mitochondria.
##============================================================================##


args = commandArgs(trailingOnly = TRUE) # Set arguments from the command line

samples.file = args[1] # Samples file
outfile = args[2] # Samples file

x = read.table(samples.file, header = F, sep = '\t')
x$V2 = 'M'
write.table(x, outfile, col.names = F, row.names = F, sep = '\t', quote = F)
