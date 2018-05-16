generate_snp_data <- function (map_file, ped_file) 
{
  map <- read.csv(map_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  header_row <- c("Family", "Individual", "Father", "Mother", 
                  "Sex", "Phenotype")
  snps = map[, 4]
  new_header = c(header_row, snps)
  ped <- read.csv(ped_file, sep = " ", header = FALSE, stringsAsFactors = FALSE, colClasses = 'character')
  range1 = seq(1, 6, by = 1)
  snp_data = data.frame(seq(1, nrow(ped), by = 1))
  for (i in range1) {
    snp_data[, i] = ped[, i]
  }
  range2 = seq(7, ncol(ped), by = 2)
  for (i in range2) {
    index = ((i - 7)/2) + 1
    snp_data[, index + 6] = paste(ped[, i], ped[, i + 1])
  }
  names(snp_data) <- new_header
  return(snp_data)
}


snp_data <- function (map_file, ped_file) 
{ 
  #read in map/ped files define coloumn headers
  map <- read.csv(map_file, sep = "\t", header = F, stringsAsFactors = T)
  header_row <- c("Family", "Individual", "Father", "Mother", 
                  "Sex", "Phenotype")
  snps = as.character(map[, 2])
  new_header = c(header_row, snps)
  ped <- read.csv(ped_file, sep = " ", header = F, stringsAsFactors = T, colClasses = 'character')
  #construct dataframe for snp data
  range1 = seq(1, 6, by = 1)
  snp_data = data.frame(seq(1, nrow(ped), by = 1))
  for (i in range1) {
    snp_data[, i] = ped[, i]
  }
  #colapse alleles to define participant genotypes
  range2 = seq(7, ncol(ped), by = 2)
  for (i in range2) {
    index = ((i - 7)/2) + 1
    snp_data[, index + 6] = paste0(ped[, i], '/', ped[, i + 1])
  }
  #apply header to SNP data
  names(snp_data) <- new_header
  #replacing 0/0 with NA to define missing data
  snp_data[snp_data == '0/0'] <- NA
  return(snp_data)
}

mito_data <- function (map_file, ped_file) 
{
  map <- read.csv(map_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  header_row <- c("Family", "Individual", "Father", "Mother", 
                  "Sex", "Phenotype")
  snps = map[, 4]
  new_header = c(header_row, snps)
  ped <- read.csv(ped_file, sep = " ", header = FALSE, stringsAsFactors = FALSE, colClasses = 'character')
  range1 = seq(1, 6, by = 1)
  snp_data = data.frame(seq(1, nrow(ped), by = 1))
  for (i in range1) {
    snp_data[, i] = ped[, i]
  }
  range2 = seq(7, ncol(ped), by = 2)
  for (i in range2) {
    index = ((i - 7)/2) + 1
    snp_data[, index + 6] = paste0(ped[, i], '/', ped[, i + 1])
  }
  names(snp_data) <- new_header
  snp_data[snp_data == '0/0'] <- NA
  return(snp_data)
}
