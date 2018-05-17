suppressPackageStartupMessages(library(tidyverse))

args = commandArgs(trailingOnly = TRUE) # Set arguments from the command line

full.file = args[1] # List of sequences
exclude.file = args[2] # List of sequences to exclude
out.file = args[3] # SPECIFY THE OUTPUT FILE
message("")
if (is.na(out.file) == TRUE | is.null(out.file) == TRUE) {
  message('OUTPUT FILE NOT SPECIFIED')
  out.file = paste0(sub("^([^.]*).*", "\\1", full.file), "_FINAL.txt")
  message(paste('OUTPUT FILE ASSIGNED TO ', out.file))
} else {
  message(paste('OUTPUT FILE ASSIGNED TO ', out.file))
}
message("")

full.obj = read_table(full.file,  col_names = F)
exclude.obj = read_table(exclude.file,  col_names = F)

exclude.obj = unique(exclude.obj)
exclude.obj = exclude.obj$X1[exclude.obj$X1 %in% full.obj$X1]
message("EXCLUDING ", length(exclude.obj), " SEQUENCES")

final = data.frame(full.obj$X1[!full.obj$X1 %in% exclude.obj])
message("SEQUENCES REMAINING AFTER FILTERING: ", nrow(final))

write.table(final, out.file, col.names = F, row.names = F, quote = F, sep = "\t")
