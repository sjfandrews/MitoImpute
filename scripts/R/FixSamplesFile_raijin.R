
## Read in Arguments 
# 1: Current working directory 
args = commandArgs(trailingOnly=TRUE)
file.in <- args[1]

dat <- read.table(file.in, sep = " ", header = T)
dat[1,4] <- "D"
write.table(dat, file.in, sep = " ", col.names = T, row.names = F, quote = F)
message(paste0(file.in, " SUCCESSFULLY CORRECTED"))