args <- commandArgs(trailingOnly = TRUE)

model_file <- args[1] # PHS file
inverse <- args[2]

library(tidyverse)
library(readxl)

source("src/config.R")
source("src/functions.R")

all_data <-
  read_excel("/ess/p212/home/p212-andrewmo/FW_Poly_tables/poly_final_220425_estimated_lastobs.xls",
    sheet = 1
  )
case_data <-
  read_excel("/ess/p212/home/p212-andrewmo/FW_Poly_tables/poly_final_220425_estimated_lastobs.xls",
    sheet = 2
  )
cancer_reg <-
  read_delim("/ess/p212/data/durable/brca2_cancer_registry/cancer_registry_20180628.tsv",
    delim = ";"
  )
id_key <-
  read_excel("/ess/p212/home/p212-andrewmo/FW_Poly_tables/EDTA_REGNR_06.05.2022.xlsx")
decode_variants <-
  read_excel("/ess/p212/home/p212-andrewmo/heritable_cancer/data/raw/deCODE_variant_data.coding_only.20250408.xlsx",
    sheet = 3
  )
sample_phenotype_map <-
  read_excel("/ess/p212/home/p212-andrewmo/heritable_cancer/data/raw/deCODE_variant_data.coding_only.20250408.xlsx",
    sheet = 1
  )
rare_variants <-
  read_excel("/ess/p212/home/p212-andrewmo/FW_Poly_tables/poly_final_220425_estimated_lastobs.xls",
    sheet = 3
  )

## Output
###############################################################################
# pid = Subject ID - pid, regnr
# bfile_id = bfile id
# age_lastobs = age last observed
# age_diagnosis = Age of diagnosis
# age = age of diagnosis for cases and age last obs for controls
# age_of_metastasisrite
# metastasis_code = raw Metastasis code (per date)
# metastasis = metastasis case/control 1/0
# conversion_to_metastasis_years = metastasis cases only, metastasis age - diagnosis age
# brca2_carrier = BRCA2 carrier - 0/1
# hoxb13_carrier = HOXB13 carrier - 0/1
# xxxx_carrier = carrier of other genes
# hoxb13_snp = HOXB13 carrier from imputed genotype data
# wgs = 0/1
# PHS - calculated later
# age_at_death = Death date
# diagnosis_to_death_years = death age minus diagnosis age
# doed_av_pca = Cause of death 0/1
# death = 0/1
# status = diagnosis 0/1
#################################################################################
# Filter by ICD10/ICD9 code, male

################################################
# Outline

# 1. filter dataset for metadata
# -1. Save list of bfile IDs in the right format
##################################################

# 1. filter prostate cancer cases/controls
male_data <-
  all_data %>%
  filter(SEX == "M") %>%
  mutate(WGS = if_else(is.na(WGS), 0, WGS)) %>%
  select(pid = PID, death = DEATH, wgs = WGS, age = AGELASTOBS)

cases <-
  case_data %>%
  filter(icd9 == 185) %>%
  mutate(status = 1) %>%
  select(pid, age_diagnosis = age_ca, status) %>%
  group_by(pid) %>%
  arrange(desc(age_diagnosis)) %>%
  filter(age_diagnosis == min(age_diagnosis)) %>%
  ungroup()

combined_data <-
  all_data %>%
  filter(SEX == "M") %>%
  mutate(WGS = if_else(is.na(WGS), 0, WGS)) %>%
  select(pid = PID, death = DEATH, wgs = WGS, age_lastobs = AGELASTOBS) %>%
  full_join(cases, by = "pid") %>%
  mutate(status = if_else(is.na(status), 0, status)) %>%
  mutate(age = if_else(status == 1, age_diagnosis, age_lastobs)) %>%
  mutate(pid = as.character(pid))

# SubjectIDs in a text file
id_key <-
  id_key %>%
  select(bfile_id = EDTANR, pid = PTN_REGNR) %>%
  mutate(bfile_id = paste("Eivind", bfile_id, sep = "_"))

# First, get the full dataset. Then filter it to the right subjects?
combined_data <-
  combined_data %>%
  left_join(id_key, by = "pid")

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
metastasis <-
  cancer_reg %>%
  filter(ICD10_gr == "C61") %>%
  select(pid = regnr, ALDER, METASTASE) %>%
  mutate(pid = as.character(pid)) %>%
  mutate(metastasis = if_else(METASTASE %in% metastase_strict, 1, 0)) %>%
  mutate(age_of_metastasis = if_else(metastasis == 1, ALDER, NA)) %>%
  rename(metastasis_code = METASTASE) %>%
  select(-ALDER)

combined_data <-
  combined_data %>%
  left_join(metastasis, by = join_by(pid)) %>%
  mutate(metastasis = if_else(is.na(metastasis), 0, metastasis)) %>%
  mutate(age_of_metastasis = if_else(is.na(age_of_metastasis), age_lastobs, age_of_metastasis)) %>%
  mutate(conversion_to_metastasis_years = age_of_metastasis - age_diagnosis)

##############################################################################
# cause of death
#############################################################################

# Died of PCA
doed_av_pca <-
  cancer_reg %>%
  filter(DAAR_CANCER == "C61") %>%
  select(pid = regnr, FOEDT_AAR, STATUSAAR) %>%
  mutate(doed_av_pca = 1) %>%
  mutate(age_at_death = STATUSAAR - FOEDT_AAR) %>%
  mutate(pid = as.character(pid)) %>%
  distinct() %>%
  select(pid, doed_av_pca, age_at_death)

combined_data <-
  combined_data %>%
  left_join(doed_av_pca, by = "pid") %>%
  mutate(diagnosis_to_death_years = age_at_death - age_diagnosis) %>%
  mutate(pca_death = if_else(is.na(doed_av_pca), 0L, 1L)) %>%
  mutate(age_at_pca_death = if_else(pca_death == 1L, age_at_death, age_lastobs))

#####################################################################
# WGS
####################################################################
decode_variants <-
  sample_phenotype_map %>%
  select(deCODE_ID, EID) %>%
  full_join(decode_variants) %>%
  mutate(bfile_id = paste0("Eivind_", EID))

variants <-
  c(
    "BRCA1", "BRCA2", "ATM", "PALB2", "CHEK2", "HOXB13", "MLH1", "MSH2",
    "MSH6", "PMS2", "TP53", "EPCAM"
  )

all_decode_variants <-
  decode_variants %>%
  filter(CPSR_CLASSIFICATION_SOURCE == "ClinVar") %>%
  filter(FINAL_CLASSIFICATION == "Pathogenic") %>%
  filter(SYMBOL %in% variants) %>%
  select(bfile_id, SYMBOL, GENOMIC_CHANGE, ALTERATION, DBSNP_RSID)

for (i in variants) {
  new <- if_else(all_decode_variants$SYMBOL == i, 1, NA)
  all_decode_variants[, ncol(all_decode_variants) + 1] <- new
  colnames(all_decode_variants)[ncol(all_decode_variants)] <-
    paste(i, "carrier", sep = "_")
}

carriers <-
  all_decode_variants %>%
  select(-SYMBOL) %>%
  filter(if_any(-bfile_id, ~ !is.na(.)))

combined_data <-
  combined_data %>%
  left_join(carriers, by = join_by(bfile_id))

# Hoxb13 SNP
raw_file <- "tmp_rs138213197.raw"

dat <-
  read.table(
    raw_file,
    header = TRUE, stringsAsFactors = FALSE, check.names = FALSE
  )

# Find the SNP column name (PLINK typically uses the SNP ID as the column name)
snp_colname <- grep("^rs138213197$|rs138213197", colnames(dat), value = TRUE)

if (length(snp_colname) != 1) {
  stop(
    "SNP column not found or not unique: ",
    paste(snp_colname, collapse = ", ")
  )
}

result <- data.frame(
  bfile_id = dat$IID, # IID from the .fam is the sample identifier
  hoxb13_snp = dat[[snp_colname]],
  stringsAsFactors = FALSE
)

# Convert PLINK missing code (-9) to NA if present
result$hoxb13_snp[result$hoxb13_snp == -9] <- NA

# Ensure integer type 0/1/2
result$hoxb13_snp <- as.integer(result$hoxb13_snp)

combined_data <-
  combined_data %>%
  left_join(result, by = join_by(bfile_id))

# Add PHS
phs_data <- import_phs_data(model_file, inverse)

combined_data <- process_combined_data(combined_data, phs_data)

# Save output
write_tsv(combined_data, "phenotypes.tsv")

# List all HOXB13 variants not classified as pathogenic in clinvar
all_variants <-
  decode_variants %>%
  filter(bfile_id %in% combined_data$bfile_id) %>%
  filter(SYMBOL %in% variants) %>%
  select(
    SYMBOL, GENOMIC_CHANGE, ALTERATION, DBSNP_RSID, FINAL_CLASSIFICATION
  ) %>%
  group_by(
    SYMBOL, GENOMIC_CHANGE, ALTERATION, DBSNP_RSID, FINAL_CLASSIFICATION
  ) %>%
  summarize(n = n())

write_tsv(all_variants, "all_variants.tsv")

#########################################################################
## Summary statistics: median and IQR
birth_year_p212 <-
  cancer_reg %>%
  filter(ICD10_gr == "C61") %>%
  select(pid = regnr, birth_year = FOEDT_AAR) %>%
  mutate(pid = as.character(pid)) %>%
  distinct()

age_summary <- combined_data %>%
  left_join(birth_year_p212, by = "pid") %>%
  mutate(
    age_first_follow_up = NA_real_,
    age_of_metastasis = if_else(metastasis == 1, age_of_metastasis, NA_real_),
    age_lastobs_censored = if_else(status == 0, age_lastobs, NA_real_)
  ) %>%
  summarize(
    across(
      c(birth_year, age_first_follow_up, age_diagnosis, age_of_metastasis,
        age_lastobs_censored),
      list(
        median = ~ median(.x, na.rm = TRUE),
        q25    = ~ quantile(.x, 0.25, na.rm = TRUE),
        q75    = ~ quantile(.x, 0.75, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    )
  ) %>%
  pivot_longer(
    everything(),
    names_to  = c("variable", ".value"),
    names_sep = "_(?=median|q25|q75)"
  ) %>%
  mutate(iqr = q75 - q25) %>%
  select(variable, median, q25, q75, iqr)

write_tsv(age_summary, "age_summary.tsv")