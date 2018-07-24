
## Read in Arguments 
# 1: Current working directory 
args = commandArgs(trailingOnly=TRUE)
file.path <- args[1]

dat <- read.table(file.path, sep = " ", header = T)
dat[1,4] <- "D"
write.table(dat, file.path, sep = " ", col.names = T, row.names = F, quote = F)