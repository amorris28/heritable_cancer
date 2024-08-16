library(tidyverse)
library(genio)

#######################
# First argument is a txt file containing a list of sumstat filenames
# Second argument is the basename for your Bfiles
# Be sure to filter the RSIds in your Bfiles so as not to load
# a gigantic genotype data file into memory

#######################
# First argument is a sumstat file
# Second argument is the basename for your Bfiles
# Be sure to pre-filter the RSIds in your Bfiles so as not to load
# a gigantic genotype data file into memory

args = commandArgs(trailingOnly=TRUE)

sumstat_files = read_lines(args[1])
sumstats <- lapply(sumstat_files, read_tsv)
names(sumstats) <- sumstat_files

antegenes = read_plink(args[2])$X
seibert_notin_ante = sumstats[[1]]$RSID[!(sumstats[[1]]$RSID %in% rownames(antegenes))]
write_lines(seibert_notin_ante, "test.txt")


# Check SNPs not in genotype data
#seibert_notin_moba = sumstats[[1]]$RSID[!(sumstats[[1]]$RSID %in% rownames(moba))]
seibert_notin_ante = sumstats[[1]]$RSID[!(sumstats[[1]]$RSID %in% rownames(antegenes))]
# write_lines(seibert_notin_moba, "seibert_notin_moba2")

# Replace the missing SNP with an RSID in LD
sumstats[[1]][sumstats[[1]]$RSID == seibert_notin_ante[1], "RSID"] = "rs10905266"

# Pre-process data so RSIDs match in sumstats and genotype
filter_sort_sumstats = function(sumstats, genotype) {
  sumstats_filtered = filter(sumstats, RSID %in% rownames(genotype))
  sumstats_filtered_sorted = sumstats_filtered[order(sumstats_filtered$RSID), ]
  return(sumstats_filtered_sorted)
}

sumstats_filtered_sorted <- lapply(sumstats, filter_sort_sumstats, antegenes)

filter_sort_genotype = function(sumstats, genotype) {
  genotype_filtered = genotype[sumstats$RSID,]
  genotype_filtered_sorted = t(genotype_filtered[sort(rownames(genotype_filtered)), ])
  return(genotype_filtered_sorted)
}

genotype_filtered_sorted = lapply(sumstats_filtered_sorted, filter_sort_genotype, antegenes)

# Some SNPs from PHS models missing in genotype data imputed from moba
impute_col_means = function(x) {
  for (i in seq_len(ncol(x))) {
    x[is.na(x[, i]), i] = mean(x[, i], na.rm = TRUE)
  }
  return(x)
}

genotype_filtered_sorted_imputed = lapply(genotype_filtered_sorted, impute_col_means)

# Check that sumstats rows and matrix columns match order
mapply(function(x, y) x$RSID != colnames(y), sumstats_filtered_sorted, genotype_filtered_sorted_imputed) %>%
  lapply(sum)

# Pull out betas into list of vectors
betas = lapply(sumstats_filtered_sorted, function(x) x$B)

# Calcuate PHSes with matrix multiplication
phses = mapply('%*%', genotype_filtered_sorted_imputed, betas, SIMPLIFY = FALSE)

# Plot histograms of PHSes
ps = lapply(phses, qplot)
n = names(ps)
sapply(n, FUN = function(nameindex, ps) { ggsave(paste0("figs/", gsub("\\..*","",nameindex), "_phs_hist.png"), ps[[nameindex]]) }, ps, simplify = FALSE, USE.NAMES = TRUE)

# Save PHSes
lapply(seq_along(phses), 
       function(i) { 
         write_tsv(as_tibble(phses[[i]], 
                             rownames = "Subject"),
                   paste0(gsub("\\..*","", names(phses)[[i]]),
                          "_phs_imputed_antegenes"),
                   col_names = FALSE) })


