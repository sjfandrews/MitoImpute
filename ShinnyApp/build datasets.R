library(tidyverse)
library(pbapply)
library(HiMC); data(nodes)
library(taRifx)

##  Function
source('~/Dropbox/Research/PostDoc/MitoWax/3_Scripts/import_snps.R', chdir = TRUE)

setwd("~/Dropbox/STRANDS")

##===============================##
##  WGS plink Files 
##===============================##
wgs.map <- '~/Dropbox/src/MitoImpute/DerivedData/ThousandGenomes/chrMT_1kg_norm_decomposed_firstAlt.map'
wgs.ped <- '~/Dropbox/src/MitoImpute/DerivedData/ThousandGenomes/chrMT_1kg_norm_decomposed_firstAlt.ped'
wgs.dat <- generate_snp_data(wgs.map, wgs.ped)  

# Assign haplogorups
MTwgs.classifications <- HiMC::getClassifications(wgs.dat) 
MTwgs.classifications <- as.tibble(MTwgs.classifications)
  
##===============================##
#  Typed only
##===============================##

# file names
typ.map <- list.files(path = "~/Dropbox/STRANDS", recursive = TRUE, pattern = "*b37.map")
typ.ped <- list.files(path = "~/Dropbox/STRANDS", recursive = TRUE, pattern = "*b37.ped")
typ.names <- typ.map %>% as.tibble() %>% separate(value, c('platform', 'file'), sep = '/')

# read in files
typ.dat <- mapply(generate_snp_data, typ.map, typ.ped, SIMPLIFY = F)

# assign haplogroups
MTtyp.classifications <- pblapply(typ.dat, HiMC::getClassifications)
MTtyp.classifications <- lapply(MTtyp.classifications, as.tibble)

# Join Typed and WGS classifications 
MT_haps.out <- lapply(MTtyp.classifications, function(x){
  out <- x %>% 
    left_join(MTwgs.classifications, by = 'Individual', suffix = c("_typ", "_wgs")) %>% 
    as.tibble()
  out
})
names(MT_haps.out) <- typ.names$platform

MT_haps <- MT_haps.out[imp.names$platform] 

saveRDS(MT_haps, "~/Dropbox/Research/PostDoc-MSSM/3_mitoWAX/3_Scripts/ShinnyApp/MT_haps.rds")

##===============================##
##  Imputed plink files 
##===============================##

# file names
imp.map <- list.files(path = "~/Dropbox/STRANDS", recursive = TRUE, pattern = "*imputed.map")
imp.ped <- list.files(path = "~/Dropbox/STRANDS", recursive = TRUE, pattern = "*imputed.ped")
imp.names <- imp.map %>% as.tibble() %>% separate(value, c('platform', 'file'), sep = '/')

# read in files
imp.dat <- mapply(generate_snp_data, imp.map, imp.ped, SIMPLIFY = F)

imp.dat <- lapply(imp.dat, function(x){
  out <- x[,-c(grep("\\<189\\>", colnames(x)), grep("\\<16183\\>", colnames(x)))]
  out
})
imp.dat <- lapply(imp.dat, as.tibble)

names(imp.dat) <- imp.names$platform

saveRDS(imp.dat, "~/Dropbox/Research/PostDoc-MSSM/3_mitoWAX/3_Scripts/ShinnyApp/imp.dat.rds")

##===============================##
##  Info Score Files 
##===============================##
info <- list.files(path = "~/Dropbox/STRANDS", recursive = TRUE, pattern = "*_info")
info.names <- info %>% as.tibble() %>% separate(value, c('platform', 'file'), sep = '/')
info.dat <- lapply(info, read_delim, delim = " ")
names(info.dat) <- info.names$platform

imp.info <- lapply(info.dat, function(x){
  out <- mutate(x, info_comb = ifelse(info_type0 == -1, info,info_type0 ))
  out <- mutate(out, himc = ifelse(position %in% c(825, 1018, 1438, 1719, 1736, 2092, 3505, 3552, 3594, 4580, 4769, 4917, 4977, 5178, 5442, 6371, 7028, 8251, 8414, 8468, 8703, 9042, 9055, 9347, 9950, 10115, 10398, 10398, 10400, 10550, 11177, 11251, 11947, 12007, 12308, 12705, 13263, 13368, 13506, 13708, 13789, 14178, 14318, 14470, 14560, 14668, 14766, 15043, 15326, 15452, 15535, 16111, 16189, 16391), 'yes', 'no'))
  out
})

saveRDS(imp.info, "~/Dropbox/Research/PostDoc-MSSM/3_mitoWAX/3_Scripts/ShinnyApp/imp.info.rds")


##===============================##
##   
##===============================##


x <- list( A=list(p=runif(5)), B=list(q=runif(5)) )
y <- list( A=list(r=runif(5)), C=list(s=runif(5)) )
merge.list(x,y)


MT_haps <- MT_haps.out[imp.names$platform] 


test <- merge.list(MT_haps[c(1,2)], imp.dat[c(1,2)])

anti_join(typ, imp, by = 'platform')


# GSA-24v1-0_A2-b37
imp.info <- read_delim('~/Dropbox/src/MitoImpute/DerivedData/ThousandGenomes/GSA-24v1-0_A2-b37/chrMT_1kg_GSA-24v1-0_A2-b37_imputed_info', delim = " ")
imp.map <- '~/Dropbox/src/MitoImpute/DerivedData/ThousandGenomes/GSA-24v1-0_A2-b37/chrMT_1kg_GSA-24v1-0_A2-b37_imputed.map'
imp.ped <- "~/Dropbox/src/MitoImpute/DerivedData/ThousandGenomes/GSA-24v1-0_A2-b37/chrMT_1kg_GSA-24v1-0_A2-b37_imputed.ped"
typ.map <- "~/Dropbox/src/MitoImpute/DerivedData/ThousandGenomes/GSA-24v1-0_A2-b37/chrMT_1kg_GSA-24v1-0_A2-b37.map"
typ.ped <- "~/Dropbox/src/MitoImpute/DerivedData/ThousandGenomes/GSA-24v1-0_A2-b37/chrMT_1kg_GSA-24v1-0_A2-b37.ped"

# Human610-Quadv1_B-b37
imp.info <- read_delim('~/Dropbox/src/MitoImpute/DerivedData/ThousandGenomes/Human610-Quadv1_B-b37/chrMT_1kg_Human610-Quadv1_B-b37_imputed_info', delim = " ")
imp.map <- '~/Dropbox/src/MitoImpute/DerivedData/ThousandGenomes/Human610-Quadv1_B-b37/chrMT_1kg_Human610-Quadv1_B-b37_imputed.map'
imp.ped <- "~/Dropbox/src/MitoImpute/DerivedData/ThousandGenomes/Human610-Quadv1_B-b37/chrMT_1kg_Human610-Quadv1_B-b37_imputed.ped"
typ.map <- "~/Dropbox/src/MitoImpute/DerivedData/ThousandGenomes/Human610-Quadv1_B-b37/chrMT_1kg_Human610-Quadv1_B-b37.map"
typ.ped <- "~/Dropbox/src/MitoImpute/DerivedData/ThousandGenomes/Human610-Quadv1_B-b37/chrMT_1kg_Human610-Quadv1_B-b37.ped"

# NeuroX_15036164_A-b37
imp.info <- read_delim('~/Dropbox/src/MitoImpute/DerivedData/ThousandGenomes/NeuroX_15036164_A-b37/chrMT_1kg_NeuroX_15036164_A-b37_imputed_info', delim = " ")
imp.map <- '~/Dropbox/src/MitoImpute/DerivedData/ThousandGenomes/NeuroX_15036164_A-b37/chrMT_1kg_NeuroX_15036164_A-b37_imputed.map'
imp.ped <- "~/Dropbox/src/MitoImpute/DerivedData/ThousandGenomes/NeuroX_15036164_A-b37/chrMT_1kg_NeuroX_15036164_A-b37_imputed.ped"
typ.map <- "~/Dropbox/src/MitoImpute/DerivedData/ThousandGenomes/NeuroX_15036164_A-b37/chrMT_1kg_NeuroX_15036164_A-b37.map"
typ.ped <- "~/Dropbox/src/MitoImpute/DerivedData/ThousandGenomes/NeuroX_15036164_A-b37/chrMT_1kg_NeuroX_15036164_A-b37.ped"


##=========================##
##  Readin in plink files 
##=========================##
imp.dat <- generate_snp_data(imp.map, imp.ped)

typ.dat <- generate_snp_data(typ.map, typ.ped) 


##========================##
##  Haplogroup assignment
##    typ and wgs only
##========================##

MTtyp.classifications <- HiMC::getClassifications(typ.dat)
MTwgs.classifications <- HiMC::getClassifications(wgs.dat)

MT_haps <- MTtyp.classifications %>% 
  left_join(MTwgs.classifications, by = 'Individual', suffix = c("_typ", "_wgs")) %>% 
  as.tibble()

##========================##
##  Info Scores
##========================##
imp.info <- mutate(imp.info, info_comb = ifelse(info_type0 == -1, info,info_type0 ))
imp.info <- mutate(imp.info, himc = ifelse(position %in% c(825, 1018, 1438, 1719, 1736, 2092, 3505, 3552, 3594, 4580, 4769, 4917, 4977, 5178, 5442, 6371, 7028, 8251, 8414, 8468, 8703, 9042, 9055, 9347, 9950, 10115, 10398, 10398, 10400, 10550, 11177, 11251, 11947, 12007, 12308, 12705, 13263, 13368, 13506, 13708, 13789, 14178, 14318, 14470, 14560, 14668, 14766, 15043, 15326, 15452, 15535, 16111, 16189, 16391), 'yes', 'no'))

dat.ls <- list(MT_haps, imp.dat, imp.info)
names(dat.ls) <- c('MT_haps', 'imp.dat', 'imp.info')

#ls.out <- list()
ls.out[[3]] <- dat.ls

names(ls.out) <- c('GSA-24v1-0_A2-b37', 'Human610-Quadv1_B-b37', 'NeuroX_15036164_A-b37')

saveRDS(ls.out, "~/Dropbox/Research/PostDoc-MSSM/3_mitoWAX/3_Scripts/ShinnyApp/mttest.rds")


ls.out[[input$select]]$imp.info
rm.info <- filter(ls.out[[input$select]]$imp.info, info > info.cut)
ls.out[[input$select]]$imp.dat <- 
  ls.out[[input$select]]$imp.dat[ ,colnames(ls.out[[input$select]]$imp.dat) %in% 
                                    c('Individual', rm.info$position)]

MTimp.classifications <- HiMC::getClassifications(ls.out[[input$select]]$imp.dat)

mt.haps <- ls.out[[input$select]]$MT_haps %>% 
  left_join(MTimp.classifications, by = 'Individual') %>% 
  rename(haplogroup_wgs = haplogroup, full_path_wgs = full_path)
  














