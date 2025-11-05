args = commandArgs(trailingOnly=TRUE)

model_file = args[1] # PHS file
inverse = args[2]
output_file = args[3] # output file name

library(tidyverse)
library(readxl)

source("src/config.R")
source("src/functions.R")

all_data = read_excel("/ess/p212/home/p212-andrewmo/FW_Poly_tables/poly_final_220425_estimated_lastobs.xls", sheet = 1)
case_data = read_excel("/ess/p212/home/p212-andrewmo/FW_Poly_tables/poly_final_220425_estimated_lastobs.xls", sheet = 2)
#fam_file = read_table("data/raw/moba_imputed_all.fam", col_names = FALSE) %>% as_tibble()
cancer_reg = read_delim("/ess/p212/data/durable/brca2_cancer_registry/cancer_registry_20180628.tsv", delim = ";")
id_key = read_excel("/ess/p212/home/p212-andrewmo/FW_Poly_tables/EDTA_REGNR_06.05.2022.xlsx")

#decode_variants = read_csv("/ess/p212/home/p212-andrewmo/heritable_cancer/data/raw/deCODE_variant_data.coding_only.20250408.csv")
#sample_phenotype_map = read_csv("/ess/p212/home/p212-andrewmo/heritable_cancer/data/raw/deCODE_variant_data.coding_only.20250408.sample_phenotype_map.csv")

decode_variants = read_excel("/ess/p212/home/p212-andrewmo/heritable_cancer/data/raw/deCODE_variant_data.coding_only.20250408.xlsx", sheet = 3)
sample_phenotype_map = read_excel("/ess/p212/home/p212-andrewmo/heritable_cancer/data/raw/deCODE_variant_data.coding_only.20250408.xlsx", sheet = 1)

rare_variants = read_excel("/ess/p212/home/p212-andrewmo/FW_Poly_tables/poly_final_220425_estimated_lastobs.xls", sheet = 3)

## Output
###############################################################################
# pid = Subject ID - pid, regnr
# bfile_id = bfile id
# age_lastobs = age last observed
# age_diagnosis = Age of diagnosis
# age = age of diagnosis for cases and age last obs for controls
# age_of_metastasis
# metastasis_code = raw Metastasis code (per date)
# metastasis = metastasis case/control 1/0
# conversion_to_metastasis_years = metastasis cases only, metastasis age - diagnosis age
# brca2_carrier = BRCA2 carrier - 0/1
# hoxb13_carrier = HOXB13 carrier - 0/1
# xxxx_carrier = carrier of other genes
# hoxb13_snp = HOXB13 carrier from imputed genotype data
# wgs = WGS - 0/1
# PHS - calculated later
# age_at_death = Death date
# diagnosis_to_death_years = death age minus diagnosis age
# doed_av_pca = Cause of death 0/1
# death = Dead - 0/1
# status = diagnosis 0/1
#################################################################################
# Filter by ICD10/ICD9 code, male

################################################
# Outline

# 1. filter dataset for metadata
# -1. Save list of bfile IDs in the right format
##################################################

# 1. filter prostate cancer cases/controls

all_data %>%
  filter(SEX == "M") %>%
  mutate(WGS = if_else(is.na(WGS), 0, WGS)) %>%
  select(pid = PID, death = DEATH, wgs = WGS, age = AGELASTOBS) ->
  male_data

case_data %>%
  filter(icd9 == 185) %>%
  mutate(status = 1) %>%
  select(pid, age_diagnosis = age_ca, status) %>%
  group_by(pid) %>%
  arrange(desc(age_diagnosis)) %>%
  filter(age_diagnosis == min(age_diagnosis)) %>%
  ungroup() -> cases

all_data %>%
  filter(SEX == "M") %>%
  mutate(WGS = if_else(is.na(WGS), 0, WGS)) %>%
  select(pid = PID, death = DEATH, wgs = WGS, age_lastobs = AGELASTOBS) %>%
  full_join(cases, by = "pid") %>%
  mutate(status = if_else(is.na(status), 0, status)) %>%
  mutate(age = if_else(status == 1, age_diagnosis, age_lastobs)) %>%
  mutate(pid = as.character(pid)) ->
  combined_data

# SubjectIDs in a text file

id_key %>%
  select(bfile_id = EDTANR, pid = PTN_REGNR) %>%
  mutate(bfile_id = paste("Eivind", bfile_id, sep = "_")) ->
  id_key

# First, get the full dataset. Then filter it to the right subjects?
combined_data %>%
  left_join(id_key, by = "pid") ->
  combined_data

combined_data %>%
  mutate(bfile_id2 = bfile_id) %>%
  select(bfile_id, bfile_id2) %>%
  distinct() ->
  bfile_id_samples

write_tsv(bfile_id_samples, paste0(PROCESSED_DATA, "bfile_id_samples.txt"), col_names = FALSE)

###########################################################################
# Metastasis code
#####################################################################

# There is only one PCA diagnosis per subject in cancer registry data
cancer_reg %>% 
    filter(ICD10_gr == "C61") %>%
    group_by(regnr) %>%
    summarize(n = n()) %>%
    arrange(desc(n))

# identify age at metastasis
cancer_reg %>%
  filter(ICD10_gr == "C61") %>%
    select(pid = regnr, ALDER, METASTASE) %>%
    mutate(pid = as.character(pid)) %>%
    mutate(metastasis = if_else(METASTASE %in% metastase_strict, 1, 0)) %>%
    mutate(age_of_metastasis = if_else(metastasis == 1, ALDER, NA)) %>%
    rename(metastasis_code = METASTASE) %>%
    select(-ALDER) ->
    metastasis

combined_data %>%
  left_join(metastasis, by = join_by(pid)) %>%
  mutate(metastasis = if_else(is.na(metastasis), 0, metastasis)) %>%
  mutate(age_of_metastasis = if_else(is.na(age_of_metastasis), age_lastobs, age_of_metastasis)) %>%
  mutate(conversion_to_metastasis_years = age_of_metastasis - age_diagnosis) ->
  combined_data

##############################################################################
# cause of death
#############################################################################

# Died of PCA
cancer_reg %>%
  filter(DAAR_CANCER == "C61") %>%
  select(pid = regnr, FOEDT_AAR, STATUSAAR) %>%
  mutate(doed_av_pca = 1) %>%
  mutate(age_at_death = STATUSAAR - FOEDT_AAR) %>%
  mutate(pid = as.character(pid)) %>%
  distinct() %>%
  select(pid, doed_av_pca, age_at_death) ->
  doed_av_pca

combined_data %>%
  left_join(doed_av_pca, by = "pid") %>%
  mutate(diagnosis_to_death_years = age_at_death - age_diagnosis) ->
  combined_data

#####################################################################
# WGS
####################################################################

sample_phenotype_map %>%
    select(deCODE_ID, EID) %>%
    full_join(decode_variants) %>%
    mutate(bfile_id = paste0("Eivind_", EID)) ->
    decode_variants

decode_variants %>%
    filter(CPSR_CLASSIFICATION_SOURCE == "ClinVar") %>%
    filter(FINAL_CLASSIFICATION == "Pathogenic") %>%
    select(bfile_id, SYMBOL) ->
    all_decode_variants
variants = c('BRCA1', 'BRCA2',  'ATM', 'PALB2', 'CHEK2', 'HOXB13', 'MLH1', 'MSH2', 'MSH6', 'PMS2', 'TP53', 'EPCAM')

for(i in variants) {
  new <- if_else(all_decode_variants$SYMBOL == i, 1, NA)
  all_decode_variants[ , ncol(all_decode_variants) + 1] <- new
  colnames(all_decode_variants)[ncol(all_decode_variants)] <- paste(i, "carrier", sep = "_")
}

all_decode_variants %>%
#  mutate(brca2_carrier = if_else(SYMBOL == "BRCA2", 1, NA)) %>%
#  mutate(hoxb13_carrier = if_else(SYMBOL == "HOXB13", 1, NA)) %>%
#  mutate(brca1_carrier = if_else(SYMBOL == "BRCA1", 1, NA)) %>%
#  mutate(mlh1_carrier = if_else(SYMBOL == "MLH1", 1, NA)) %>%
#  mutate(msh2_carrier = if_else(SYMBOL == "MSH2", 1, NA)) %>%
#  mutate(msh6_carrier = if_else(SYMBOL == "MSH6", 1, NA)) %>%
#  mutate(atm_carrier = if_else(SYMBOL == "ATM", 1, NA)) %>%
  select(-SYMBOL) %>% 
  filter(if_any(-bfile_id, ~ !is.na(.))) ->
  carriers

combined_data %>%
  left_join(carriers, by = join_by(bfile_id)) ->
  combined_data

# Hoxb13 SNP

#rs138213197 = read_tsv("data/processed/rs138213197.vcf", skip = 6)
#rs138213197[, 10:ncol(rs138213197)] %>%
#  pivot_longer(cols = everything(), names_to = "bfile_id") %>%
#  mutate(pos1 = as.numeric(substr(value, 1, 1))) %>%
#  mutate(pos2 = as.numeric(substr(value, 3, 3))) %>%
#  rowwise() %>%
#  mutate(hoxb13_snp = sum(pos1, pos2)) %>%
#  select(-value, -pos1, -pos2) ->
#  hoxb13
#
#combined_data %>%
#  left_join(hoxb13, by = join_by(bfile_id)) ->
#  combined_data

# Add PHS
phs_data = import_phs_data(model_file, inverse)

combined_data = process_combined_data(combined_data, phs_data)

# Save output
write_tsv(combined_data, output_file)
