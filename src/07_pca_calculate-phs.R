library(tidyverse)
library(genio)

source("src/config.R")

# Takes a list of sumstats (or just 1 I guess)
# Takes bfile

args = commandArgs(trailingOnly=TRUE)

model_file = args[1] # Sumstats file
model = basename(args[1]) # filename for output
genotype_bfiles = args[2] # PHS data

sumstats = readr::read_tsv(model_file)

# Import genotypes
genotype = genio::read_plink(genotype_bfiles)$X

# Pre-process data so RSIDs match in sumstats and genotype
filter_sort_sumstats = function(sumstats, genotype) {
  sumstats_filtered = dplyr::filter(sumstats, RSID %in% rownames(genotype))
  sumstats_filtered_sorted = sumstats_filtered[order(sumstats_filtered$RSID), ]
  return(sumstats_filtered_sorted)
}
# sumstats_filtered_sorted <- lapply(sumstats, filter_sort_sumstats, antegenes)
sumstats_filtered_sorted = filter_sort_sumstats(sumstats, genotype)

filter_sort_genotype = function(sumstats, genotype) {
  genotype_filtered = genotype[sumstats$RSID,]
  genotype_filtered_sorted = t(genotype_filtered[sort(rownames(genotype_filtered)), ])
  return(genotype_filtered_sorted)
}
 
# genotype_filtered_sorted = lapply(sumstats_filtered_sorted, filter_sort_genotype, antegenes)
genotype_filtered_sorted = filter_sort_genotype(sumstats_filtered_sorted, genotype)

# Some SNPs from PHS models missing in genotype data imputed from moba
impute_col_means = function(x) {
  for (i in seq_len(ncol(x))) {
    x[is.na(x[, i]), i] = mean(x[, i], na.rm = TRUE)
  }
  return(x)
}
# 
# genotype_filtered_sorted_imputed = lapply(genotype_filtered_sorted, impute_col_means)
genotype_filtered_sorted_imputed = impute_col_means(genotype_filtered_sorted)

# Check that sumstats rows and matrix columns match order
# mapply(function(x, y) x$RSID != colnames(y), sumstats_filtered_sorted, genotype_filtered_sorted_imputed) %>%
#   lapply(sum)

mismatches = sum(sumstats_filtered_sorted$RSID != colnames(genotype_filtered_sorted_imputed))
if (mismatches > 0) {
  stop("There are ", mismatches, " mismatches between RSIDs in sumstats and genotype.")
}

# Pull out betas into list of vectors
# betas = lapply(sumstats_filtered_sorted, function(x) x$B)

betas = sumstats_filtered_sorted$B

# Calcuate PHSes with matrix multiplication
# phses = mapply('%*%', genotype_filtered_sorted_imputed, betas, SIMPLIFY = FALSE)

phses = genotype_filtered_sorted_imputed %*% betas

# # Save PHSes
# lapply(seq_along(phses), 
#        function(i) { 
#          write_tsv(as_tibble(phses[[i]], 
#                              rownames = "Subject"),
#                    paste0(gsub("\\..*","", names(phses)[[i]]),
#                           "_phs_imputed_antegenes"),
#                    col_names = FALSE) })
readr::write_tsv(tibble::as_tibble(phses, rownames = "Subject"), paste0(PROCESSED_DATA, gsub("\\..*","", model),
       "_phs"), col_names = FALSE)
