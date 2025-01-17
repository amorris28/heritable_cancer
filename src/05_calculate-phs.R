library(tidyverse)
library(genio)

source("src/00_config.R")

# Takes a list of sumstats (or just 1 I guess)
# Takes bfile

args = commandArgs(trailingOnly=TRUE)

model_file = args[1] # Sumstats file
model = gsub("\\..*","", basename(args[1])) # filename for output

genotype_bfiles = args[2] # PHS data

sumstats = readr::read_tsv(model_file)

# Import genotypes
genotype = genio::read_plink(genotype_bfiles)$X

# Check for missing SNPs in genotype data
genotype_filtered = genotype[rownames(genotype) %in% sumstats$RSID, ]
percent_missingness = rowSums(is.na(genotype_filtered))/ncol(genotype_filtered)
pm_tbl = as_tibble(percent_missingness, rownames = "RSID")
pm_tbl = rename(pm_tbl, percent_missingness = value)
pm_tbl = cbind(model, pm_tbl)
write_tsv(pm_tbl, paste0("reports/", model, "_missingness.tsv"))

# Pre-process data so RSIDs match in sumstats and genotype
sumstats_filtered = dplyr::filter(sumstats, RSID %in% rownames(genotype))
sumstats_filtered_sorted = sumstats_filtered[order(sumstats_filtered$RSID), ]

genotype_filtered = genotype[sumstats_filtered_sorted$RSID,]
genotype_filtered_sorted = t(genotype_filtered[sort(rownames(genotype_filtered)), ])

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
readr::write_tsv(tibble::as_tibble(phses, rownames = "Subject"), paste0(PROCESSED_DATA, model,
       "_phs"), col_names = FALSE)


# Write diagnostics

diagnostics = tibble(model,
		     missing_SNPs = sumstats$RSID[!(sumstats$RSID %in% rownames(genotype))], 
		     beta = sumstats$B[!(sumstats$RSID %in% rownames(genotype))])
readr::write_tsv(diagnostics, paste0("reports/", model, "_diagnostics.tsv"))

