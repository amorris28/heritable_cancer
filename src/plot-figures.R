library(tidyverse)
library(stringr)
library(hazrd)
library(survival)
library(patchwork)
source("src/config.R")
source("src/functions.R")

args = commandArgs(trailingOnly=TRUE)

metadata_file = args[1] # Phenotype data
output_prefix = args[2]

#########################################
# Import PHS data and metadata
combined_data = read_tsv(metadata_file)

#########################################
# Filter over 40 
combined_data %>%
  filter(status == 1) %>%
  summarize(min_age = min(age)) # 40ish

combined_data %>%
  filter(age >= 40) %>%
  mutate(age_minus_40 = age - 40) ->
  combined_data

# Plot KM Curves
any_pca_km_curve = plot_km_curves(combined_data)
any_pca_inc_curve = plot_km_curves(combined_data, style = "risk")

# KM Curves for metastatic disease
metastatic_km_curve = plot_km_curves(combined_data, age = "age_of_metastasis", status = "metastasis")
metastatic_inc_curve = plot_km_curves(combined_data, age = "age_of_metastasis", status = "metastasis", style = "risk")

# Save Output
saveRDS(any_pca_km_curve, paste(output_prefix, "any_pca_km_curve.rds", sep = "_"))
saveRDS(metastatic_km_curve, paste(output_prefix, "metastatic_km_curve.rds", sep = "_"))
saveRDS(any_pca_inc_curve, paste(output_prefix, "any_pca_km_curve.rds", sep = "_"))
saveRDS(metastatic_inc_curve, paste(output_prefix, "metastatic_km_curve.rds", sep = "_"))
