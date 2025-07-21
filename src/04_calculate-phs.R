library(tidyverse)
library(genio)

source("src/00_config.R")

# Takes a list of sumstats (or just 1 I guess)
# Takes bfile

args = commandArgs(trailingOnly=TRUE)

model_file = args[1] # Sumstats file
model = gsub("\\..*","", basename(model_file)) # filename for output

genotype_bfiles = args[2] # PHS data

tsd_proj = args[3] # PHS data

sumstats = readr::read_tsv(model_file)

# Import genotypes
genotype = genio::read_plink(genotype_bfiles)$X

# if (tsd_proj == "p212") {
# 	phenotypes = read_tsv(paste0(PROCESSED_DATA, "heritable_cancer_phenotypes.tsv"))
# 	phenotypes %>%
# 		select(bfile_id, hoxb13_snp) %>%
# 		mutate(hoxb13_snp = if_else(hoxb13_snp == 0, 1, 0)) ->
# 		hoxb13
# 	hoxb13_df = as.data.frame(hoxb13)
#     m <- as.matrix(hoxb13_df[,-1])
#     dimnames(m) <- list(hoxb13_df[,1], "rs138213197")
#     genotype = t(genotype)
# 	merged = merge(genotype, m, by = "row.names", all = TRUE)
# 	genotype = as.matrix(merged[-1])
# 	rownames(genotype) = merged[, 1]
# 	genotype = t(genotype)
# }

# Check for missing SNPs in genotype data
genotype_filtered = genotype[rownames(genotype) %in% sumstats$RSID, ]
percent_missingness = rowSums(is.na(genotype_filtered))/ncol(genotype_filtered)
pm_tbl = as_tibble(percent_missingness, rownames = "RSID")
pm_tbl = rename(pm_tbl, percent_missingness = value)
pm_tbl = cbind(model, pm_tbl)
write_tsv(pm_tbl, paste0("reports/", paste(tsd_proj, model, "missingness.tsv", sep = "_")))

# Pre-process data so RSIDs match in sumstats and genotype
sumstats_filtered = dplyr::filter(sumstats, RSID %in% rownames(genotype))
sumstats_filtered_sorted = sumstats_filtered[order(sumstats_filtered$RSID), ]

genotype_filtered = genotype[sumstats_filtered_sorted$RSID,]
genotype_filtered_sorted = t(genotype_filtered[sort(rownames(genotype_filtered)), ])

# # Some SNPs from PHS models missing in genotype data imputed from moba
# impute_col_means = function(x) {
#   for (i in seq_len(ncol(x))) {
#     x[is.na(x[, i]), i] = mean(x[, i], na.rm = TRUE)
#   }
#   return(x)
# }

# # genotype_filtered_sorted_imputed = lapply(genotype_filtered_sorted, impute_col_means)
# genotype_filtered_sorted_imputed = impute_col_means(genotype_filtered_sorted)

# Check that sumstats rows and matrix columns match order
# mapply(function(x, y) x$RSID != colnames(y), sumstats_filtered_sorted, genotype_filtered_sorted_imputed) %>%
#   lapply(sum)

mismatches = sum(sumstats_filtered_sorted$RSID != colnames(genotype_filtered_sorted))
if (mismatches > 0) {
  stop("There are ", mismatches, " mismatches between RSIDs in sumstats and genotype.")
}

# Pull out betas into list of vectors
# betas = lapply(sumstats_filtered_sorted, function(x) x$B)

betas = sumstats_filtered_sorted$B

# Calcuate PHSes with matrix multiplication
# phses = mapply('%*%', genotype_filtered_sorted_imputed, betas, SIMPLIFY = FALSE)
genotype_filtered_sorted[is.na(genotype_filtered_sorted)] = 0
phses = genotype_filtered_sorted %*% betas
# # Save PHSes
# lapply(seq_along(phses), 
#        function(i) { 
#          write_tsv(as_tibble(phses[[i]], 
#                              rownames = "Subject"),
#                    paste0(gsub("\\..*","", names(phses)[[i]]),
#                           "_phs_imputed_antegenes"),
#                    col_names = FALSE) })
readr::write_tsv(tibble::as_tibble(phses, rownames = "Subject"), 
		 paste0(PROCESSED_DATA, paste(model, "phs", sep = "_")), col_names = FALSE)


# Write diagnostics

diagnostics = tibble(model,
		     missing_SNPs = sumstats$RSID[!(sumstats$RSID %in% rownames(genotype))], 
		     beta = sumstats$B[!(sumstats$RSID %in% rownames(genotype))])
readr::write_tsv(diagnostics, paste0("reports/", paste(tsd_proj, model, "diagnostics.tsv", sep = "_")))

