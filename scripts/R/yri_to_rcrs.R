#!/usr/bin/Rscript
#-----------------------------------------------------##
## 	Converting Yoruba (YRI - AF347015) SNP positions
##	to the revised cambrige reference sequence (rCRS)
##	https://www.mitomap.org/MITOMAP/YorubanConversion
#	Yoruba		Cambridge
#	1-309		Same
#	310			Deleted
#	311-316		Subtract One
#	317			Deleted
#	318-3108	Subtract Two
#	Deleted		Placeholder at 3107
#	3109-16194	Subtract One
#	16195		Deleted
#	16196-16571	Subtract Two
#-----------------------------------------------------##

suppressPackageStartupMessages(library(tidyverse))

args = commandArgs(trailingOnly=TRUE)
filename <- args[1]
outfile <- args[2]

##  Read in vcf file
vcf <- read_tsv(filename, comment = '##')

## function for converting YRI -> rCRS
yri_rcrs <- function(x){
	if(x < 311){
	x
}else{
	if(x >= 311 & x < 318){
	x - 1
}else {
	if(x >= 318 & x < 3109){
	x - 2
}else {
	if(x >= 3109 & x < 16195){
	x - 1
}else {
	if(x >= 16196 & x < 16571){
	x - 2
}else {
	x
}}}}}
}

vcf$POS <- sapply(vcf$POS, yri_rcrs)
write_tsv(vcf, append = T, outfile, col_names = T)
