library(tidyverse)
library(readxl)

source("src/00_config.R")

all_data = read_excel("/ess/p212/home/p212-andrewmo/FW_Poly_tables/poly_final_220425_estimated_lastobs.xls", sheet = 1)
case_data = read_excel("/ess/p212/home/p212-andrewmo/FW_Poly_tables/poly_final_220425_estimated_lastobs.xls", sheet = 2)
fam_file = read_table("data/raw/moba_imputed_all.fam", col_names = FALSE) %>% as_tibble()
cancer_reg = read_delim("/ess/p212/data/durable/brca2_cancer_registry/cancer_registry_20180628.tsv", delim = ";")
id_key = read_excel("/ess/p212/home/p212-andrewmo/FW_Poly_tables/EDTA_REGNR_06.05.2022.xlsx")

decode_variants = read_csv("data/raw/deCODE_variant_data.coding_only.20250408.csv")
sample_phenotype_map = read_csv("data/raw/deCODE_variant_data.coding_only.20250408.sample_phenotype_map.csv")

rare_variants = read_excel("/ess/p212/home/p212-andrewmo/FW_Poly_tables/poly_final_220425_estimated_lastobs.xls", sheet = 3)

## Output
###############################################################################
# Subject ID - pid, regnr
# bfile id
# Age of diagnosis
# Metastasis code (per date)
# BRCA2 carrier - 0/1
# HOXB13 carrier - 0/1
# WGS - 0/1
# PHS - calculated later
# Death date
# Cause of death
# Dead - 0/1
#################################################################################
# Filter by ICD10/ICD9 code, male

# Create a second dataset with females BRCA2
# BRCA2 0/1
# age of diagnosis
# BRCA2 0/1
# WGS 0/1

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
  distinct()->
  bfile_id_samples

write_tsv(bfile_id_samples, paste0(PROCESSED_DATA, "bfile_id_samples.txt"), col_names = FALSE)

###########################################################################
# Metastasis code
#####################################################################

# metastasis, 50 cases with strict and 100 with loose so I am gonna go with loose
#metastase_strict = c(1,2,4,6,7)
metastase_loose = c(1,2,4,5,6,7)

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
    # mutate(metastase_strict = if_else(METASTASE %in% metastase_strict, 1, 0)) %>%
    # mutate(metastase_loose = if_else(METASTASE %in% metastase_loose, 1, 0)) %>%
    mutate(metastasis = if_else(METASTASE %in% metastase_loose, 1, 0)) %>%
    rename(age_of_metastasis = ALDER,
           metastasis_code = METASTASE) ->
    metastasis

# sum(metastasis$metastase_loose)
# sum(metastasis$metastase_strict)

# nrow(metastasis)

combined_data %>%
  left_join(metastasis, by = join_by(pid)) %>%
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

all_decode_variants %>%
  mutate(brca2_carrier = if_else(SYMBOL == "BRCA2", 1, NA)) %>%
  mutate(hoxb13_carrier = if_else(SYMBOL == "HOXB13", 1, NA)) %>%
  mutate(brca1_carrier = if_else(SYMBOL == "BRCA1", 1, NA)) %>%
  mutate(mlh1_carrier = if_else(SYMBOL == "MLH1", 1, NA)) %>%
  mutate(msh2_carrier = if_else(SYMBOL == "MSH2", 1, NA)) %>%
  mutate(msh6_carrier = if_else(SYMBOL == "MSH6", 1, NA)) %>%
  mutate(atm_carrier = if_else(SYMBOL == "ATM", 1, NA)) %>%
  select(-SYMBOL) %>% 
  filter(if_any(-bfile_id, ~ !is.na(.))) ->
  brca2_hoxb13_carriers

combined_data %>%
  left_join(brca2_hoxb13_carriers, by = join_by(bfile_id)) ->
  combined_data

###################################################################
# Female subjects with BRCA2 and/or breast cancer


# BRCA2 breast cancer data

# Get case IDs and diagnosis date
case_data %>%
    filter(icd9 == 174) %>%
    select(pid, age_diag = age_ca) %>%
    mutate(status = 1) %>%
    mutate(pid = as.character(pid)) ->
    bca_cases

# Get control IDs and date last obs
all_data %>%
    filter(SEX == "F") %>%
    select(pid = PID, age_lob = AGELASTOBS, wgs = WGS) %>%
    mutate(wgs = if_else(is.na(wgs), 0, wgs)) %>%
    mutate(pid = as.character(pid)) ->
    female_subjects

# link cases and controls
female_subjects %>%
    full_join(bca_cases, by = "pid") %>%
    mutate(age = if_else(is.na(age_diag), age_lob, age_diag)) %>%
    mutate(status = if_else(is.na(status), 0, status)) %>%
    left_join(id_key, by = join_by(pid)) ->
    bca_combined_data

# Pull out BRCA2 carriers

all_decode_variants %>%
  mutate(brca2_carrier = if_else(SYMBOL == "BRCA2", 1, NA)) %>%
  select(-SYMBOL) %>%
  filter(!is.na(brca2_carrier)) ->
  brca2_carriers

bca_combined_data %>%
  left_join(brca2_carriers, by = join_by(bfile_id)) %>%
    mutate(brca2_carrier = if_else(is.na(brca2_carrier), 0 , brca2_carrier)) ->
    bca_combined_data

# There's only 4 brca2 carriers??
bca_combined_data %>%
  # filter(wgs == 1) %>%
  select(status, brca2_carrier) %>%
  table()

####################################################################
# Hoxb13 SNP

rs138213197 = read_tsv("data/processed/rs138213197.vcf", skip = 6)
rs138213197[, 10:ncol(rs138213197)] %>%
  pivot_longer(cols = everything(), names_to = "bfile_id") %>%
  mutate(pos1 = as.numeric(substr(value, 1, 1))) %>%
  mutate(pos2 = as.numeric(substr(value, 3, 3))) %>%
  rowwise() %>%
  mutate(hoxb13_snp = sum(pos1, pos2)) %>%
  select(-value, -pos1, -pos2) ->
  hoxb13

combined_data %>%
  left_join(hoxb13, by = join_by(bfile_id)) ->
  combined_data

# Save output
write_tsv(combined_data, paste0(PROCESSED_DATA, "heritable_cancer_phenotypes.tsv"))
write_tsv(bca_combined_data, paste0(PROCESSED_DATA, "brca2_breast_cancer_metadata.tsv"))

#####################################################################
#####################################################################
#####################################################################
#####################################################################

# # Grab IDs of prostate cancer cases
# ## poly_final
# caseData %>% 
#   filter(icd9 == 185) %>% 
#   pull(pid) %>% 
#   unique ->
#   caseIDs

# ## all cancer cases
# caseData %>%
#   pull(pid) %>%
#   unique() ->
#   allCases

# ## cancer registry
# cancer_reg %>%
#   filter(ICD10_gr == "C61") %>%
#   pull(regnr) %>%
#   unique ->
#   crIDs

# paste0(round(sum(caseIDs %in% crIDs) / length(caseIDs) * 100), "% of pca cases are in the cancer registry data")
# paste0(round(sum(crIDs %in% caseIDs) / length(crIDs) * 100), "% of cr pca cases are in the heritable cancer data")

# ## allIDs
# allData %>%
#   pull(PID) %>%
#   unique ->
#   allIDs

# ## all IDs cancer registry
# cancer_reg %>%
#   pull(regnr) %>%
#   unique ->
#   allcrIDs

# paste0(round(sum(allCases %in% allcrIDs) / length(allCases) * 100), "% of all cases are in the cancer registry data")



# # Grab IDs for all cancer-free males
# allData %>% 
#   filter(SEX == "M") %>% 
#   filter(!(PID %in% caseIDs)) %>% 
#   pull(PID) %>% 
#   unique() ->
#   controlIDs

# case = tibble(PTN_REGNR = as.character(caseIDs))
# control = tibble(PTN_REGNR = as.character(controlIDs))

# allData %>%
# 	select(PTN_REGNR = PID) %>%
# 	mutate(PTN_REGNR = as.character(PTN_REGNR)) %>%
# 	left_join(id_key, by = "PTN_REGNR") %>%
#   mutate(bfile_id = paste("Eivind", EDTANR, sep = "_")) ->
#   all_samples

# n_missing_in_geno = sum(!(all_samples$bfile_id %in% fam_file$X1)) # 564
# n_missing_in_pheno = sum(!(fam_file$X1 %in% all_samples$bfile_id)) # 0

# missing_subjects = all_samples$bfile_id[!(all_samples$bfile_id %in% fam_file$X1)]

# case %>%
#   left_join(id_key, by = "PTN_REGNR") %>% 
#   mutate(bfile_id = paste("Eivind", EDTANR, sep = "_")) ->
#   case_samples

# control %>%
#   left_join(id_key, by = "PTN_REGNR") %>% 
#   mutate(bfile_id = paste("Eivind", EDTANR, sep = "_")) ->
#   control_samples

# sum(case_samples$bfile_id %in% missing_subjects) # 21
# sum(control_samples$bfile_id %in% missing_subjects) # 108

# bind_rows(case_samples, control_samples) %>%
# 	mutate(bfile_id2 = bfile_id) %>% 
# 	select(bfile_id, bfile_id2) ->
# 	bfile_id_samples

# bind_rows(case, control) %>% 
#   left_join(id_key, by = "PTN_REGNR") %>% 
#   mutate(bfile_id = paste("Eivind", EDTANR, sep = "_")) ->
#   samples

# samples %>% 
#   mutate(bfile_id2 = bfile_id) %>% 
#   select(bfile_id, bfile_id2) ->
#   bfile_id_samples

# # This is necessary if using the antegenes imputation because the family IDs are zero
# # whereas with the MOBA/HRC imputation, family IDs = individual IDs
# #bfile_id_samples$bfile_id = 0

# # Save bfile ids for later filtering of genotype data
# write_tsv(bfile_id_samples, paste0(PROCESSED_DATA, "bfile_id_samples.txt"), col_names = FALSE)

# ########################################################
# # bfile_id, age, status
# # Eivind_00001, 56, 0/1
# caseData %>%
#   filter(icd9 == 185) %>%
#   rename(PTN_REGNR = pid) %>%
#   rename(age = age_ca) %>%
#   mutate(status = 1) %>%
#   select(PTN_REGNR, age, status) %>%
#   group_by(PTN_REGNR) %>%
#   filter(age == min(age)) %>%
#   ungroup() -> cases

# allData %>%
#   filter(PID %in% controlIDs) %>%
#   rename(PTN_REGNR = PID, age = AGELASTOBS) %>%
#   select(PTN_REGNR, age) %>%
#   mutate(status = 0)
#   group_by(PTN_REGNR) -> controls

# bind_rows(cases, controls) %>%
#   mutate(PTN_REGNR = as.character(PTN_REGNR)) ->
#   combined_data


# ################################################################
# # breast cancer cases and controls for BRCA2 analysis
# caseData %>%
#   filter(icd9 == 174) %>%
#   rename(PTN_REGNR = pid) %>%
#   rename(age = age_ca) %>%
#   mutate(status = 1) %>%
#   select(PTN_REGNR, age, status) %>%
#   group_by(PTN_REGNR) %>%
#   filter(age == min(age)) %>%
#   ungroup() -> brcaCases

# allData %>% 
#   filter(SEX == "F") %>% 
#   filter(!(PID %in% brcaCases$PTN_REGNR)) %>% 
#   pull(PID) %>% 
#   unique() ->
#   brcaControlIDs

# allData %>%
#   filter(PID %in% brcaControlIDs) %>%
#   rename(PTN_REGNR = PID, age = AGELASTOBS) %>%
#   select(PTN_REGNR, age) %>%
#   mutate(status = 0) %>%
#   group_by(PTN_REGNR) -> brcaControls

# bind_rows(brcaCases, brcaControls) %>%
#   mutate(PTN_REGNR = as.character(PTN_REGNR)) ->
#   brca_combined_data

# all_samples %>%
#   select(bfile_id, PTN_REGNR) %>%
#   right_join(brca_combined_data, by = "PTN_REGNR") %>%
#   select(-PTN_REGNR) -> brca_phenotypes

# allData %>%
#   filter(PID %in% controlIDs) %>%
#   rename(PTN_REGNR = PID, age = AGELASTOBS) %>%
#   select(PTN_REGNR, age) %>%
#   mutate(status = 0) %>%
#   group_by(PTN_REGNR) -> controls

# bind_rows(cases, controls) %>%
#   mutate(PTN_REGNR = as.character(PTN_REGNR)) ->
#   combined_data
