args = commandArgs(trailingOnly = TRUE) # Set arguments from the command line

samples.file = args[1] # Samples file

x = read.table(samples.file, header = F, sep = '\t')
x$V2 = 'M'
write.table(x, samples.file, col.names = F, row.names = F, sep = '\t', quote = F)
