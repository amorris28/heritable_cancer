args <- commandArgs(trailingOnly = TRUE)

model_file <- args[1] # PHS file
inverse <- args[2]
genotype_bfile <- args[3]
phs_bfile <- args[4]
sample_counts_file <- args[5]

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

show_duplicate_keys(
  filter(all_data, SEX == "M"), "PID", "all_data (males)"
)
show_duplicate_keys(
  filter(case_data, icd9 == 185), "pid", "case_data PCa rows"
)

cases <-
  case_data %>%
  filter(icd9 == 185) %>%
  mutate(status = 1) %>%
  select(pid, age_diagnosis = age_ca, status) %>%
  group_by(pid) %>%
  arrange(desc(age_diagnosis)) %>%
  filter(age_diagnosis == min(age_diagnosis)) %>%
  ungroup()

# Registry (ICD10 C61) is a second source of PCa diagnoses and is used below
# for metastasis and PCa death, so it also has to define case status.
reg_pca <-
  cancer_reg %>%
  filter(ICD10_gr == "C61") %>%
  mutate(pid = as.character(regnr))

show_duplicate_keys(
  select(reg_pca, pid, ALDER, METASTASE, FOEDT_AAR, STATUSAAR, DAAR_CANCER),
  "pid", "registry C61 rows"
)

registry_dx <-
  reg_pca %>%
  filter(!is.na(ALDER)) %>%
  group_by(pid) %>%
  summarize(age_diagnosis_reg = min(ALDER)) %>%
  ungroup()

combined_data <-
  male_data %>%
  full_join(cases, by = "pid") %>%
  mutate(pid = as.character(pid)) %>%
  left_join(registry_dx, by = "pid") %>%
  mutate(
    in_case_data = !is.na(status),
    in_registry = !is.na(age_diagnosis_reg),
    diagnosis_source = case_when(
      in_case_data & in_registry ~ "both",
      in_case_data ~ "case_data only",
      in_registry ~ "registry only",
      TRUE ~ "none"
    ),
    # A registry C61 diagnosis counts as a case even if sheet 2 lacks it.
    status = if_else(in_case_data | in_registry, 1, 0),
    age_diagnosis = coalesce(age_diagnosis, age_diagnosis_reg),
    age = if_else(status == 1, age_diagnosis, age_lastobs)
  )

message("PCa diagnosis source:")
print(count(combined_data, diagnosis_source))
registry_only <- filter(combined_data, diagnosis_source == "registry only")
if (nrow(registry_only) > 0) {
  message("Registry C61 diagnoses missing from case_data (treated as cases):")
  print(registry_only, n = Inf, width = Inf)
}
age_mismatch <-
  filter(combined_data, diagnosis_source == "both", abs(age_diagnosis - age_diagnosis_reg) > 1)
if (nrow(age_mismatch) > 0) {
  message("Diagnosis age differs >1 year between case_data and registry:")
  print(age_mismatch, n = Inf, width = Inf)
}
combined_data <- select(combined_data, -in_case_data, -in_registry, -age_diagnosis_reg, -diagnosis_source)

# SubjectIDs in a text file
id_key <-
  id_key %>%
  select(bfile_id = EDTANR, pid = PTN_REGNR) %>%
  mutate(bfile_id = paste("Eivind", bfile_id, sep = "_")) %>%
  mutate(pid = as.character(pid)) %>%
  resolve_duplicate_ids(c("pid", "bfile_id"), "p212 id_key")

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
  reg_pca %>%
  select(pid, ALDER, METASTASE) %>%
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

# Registry ages and case_data ages can disagree by a year. An event dated up to
# 1 year before the diagnosis is set to the diagnosis age; larger gaps are
# real conflicts and are caught by the checks below.
clamp_to_diagnosis <- function(event_age, event, diagnosis_age, tolerance = 1) {
  early <- event == 1 & !is.na(diagnosis_age) & event_age < diagnosis_age &
    diagnosis_age - event_age <= tolerance
  if_else(early, diagnosis_age, event_age)
}

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

n_clamped <- c(
  metastasis = sum(combined_data$metastasis == 1 &
    combined_data$age_of_metastasis < combined_data$age_diagnosis, na.rm = TRUE),
  pca_death = sum(combined_data$pca_death == 1 &
    combined_data$age_at_pca_death < combined_data$age_diagnosis, na.rm = TRUE)
)
message(
  "Events dated before diagnosis (<=1 year are set to the diagnosis age): ",
  paste(names(n_clamped), n_clamped, sep = "=", collapse = ", ")
)
combined_data <-
  combined_data %>%
  mutate(
    age_of_metastasis = clamp_to_diagnosis(age_of_metastasis, metastasis, age_diagnosis),
    age_at_pca_death = clamp_to_diagnosis(age_at_pca_death, pca_death, age_diagnosis),
    age_at_death = clamp_to_diagnosis(age_at_death, pca_death, age_diagnosis),
    conversion_to_metastasis_years = age_of_metastasis - age_diagnosis,
    diagnosis_to_death_years = age_at_death - age_diagnosis
  )

# Consistency checks between sources.
impossible <- filter(combined_data, status == 0 & (metastasis == 1 | pca_death == 1))
if (nrow(impossible) > 0) {
  print(impossible, n = Inf, width = Inf)
  stop("Controls with metastasis or PCa death: ", nrow(impossible), " subjects")
}
inconsistent <- bind_rows(
  filter(combined_data, metastasis == 1, age_of_metastasis < age_diagnosis) %>%
    mutate(problem = "metastasis before diagnosis"),
  filter(combined_data, pca_death == 1, age_at_death < age_diagnosis) %>%
    mutate(problem = "PCa death before diagnosis"),
  filter(combined_data, pca_death == 1, death != 1) %>%
    mutate(problem = "registry PCa death but all_data death != 1")
)
if (nrow(inconsistent) > 0) {
  message("Inconsistent records between sources:")
  print(inconsistent, n = Inf, width = Inf)
}
stopifnot(
  "metastasis or PCa death more than 1 year before diagnosis" =
    !any(inconsistent$problem %in% c("metastasis before diagnosis", "PCa death before diagnosis"))
)

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
raw_file <- "data/processed/tmp_rs138213197.raw"

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

combined_data <- resolve_duplicate_ids(combined_data, c("pid", "bfile_id"), "p212 phenotypes")

genotype_ids <- read_fam_iids(genotype_bfile)
phs_bfile_ids <- read_fam_iids(phs_bfile)
male_ids <- unique(as.character(male_data$pid))
male_bfile_ids <- unique(as.character(id_key$bfile_id[match(male_ids, id_key$pid)]))
male_bfile_ids <- male_bfile_ids[!is.na(male_bfile_ids)]
wgs_full_cohort_n <- all_data %>%
  filter(WGS == 1) %>%
  summarize(n = n_distinct(PID, na.rm = TRUE)) %>%
  pull(n)
wgs_male_n <- male_data %>%
  filter(wgs == 1) %>%
  summarize(n = n_distinct(pid, na.rm = TRUE)) %>%
  pull(n)
sample_counts <- sample_count_rows(
  "NFCS",
  stage = c(
    "source_cohort",
    "source_males",
    "post_qc_genotype_file_samples",
    "genotype_qc_passed_cohort",
    "male_with_genotype_qc",
    "phs_variant_bfile_samples",
    "phs_variant_bfile_cohort",
    "phs_score_file_samples",
    "male_with_phs_score",
    "male_with_phs_and_complete_phenotype",
    "male_with_phs_and_complete_phenotype_age40plus",
    "wgs_flagged_source_cohort",
    "wgs_flagged_males",
    "wgs_males_with_phs_and_complete_phenotype",
    "wgs_males_with_phs_and_complete_phenotype_age40plus"
  ),
  n = c(
    dplyr::n_distinct(all_data$PID, na.rm = TRUE),
    dplyr::n_distinct(male_data$pid, na.rm = TRUE),
    length(genotype_ids),
    length(intersect(unique(as.character(id_key$bfile_id)), genotype_ids)),
    length(intersect(male_bfile_ids, genotype_ids)),
    length(phs_bfile_ids),
    length(intersect(unique(as.character(id_key$bfile_id)), phs_bfile_ids)),
    dplyr::n_distinct(phs_data$bfile_id, na.rm = TRUE),
    length(intersect(male_bfile_ids, unique(as.character(phs_data$bfile_id)))),
    nrow(combined_data),
    sum(combined_data$age >= 40),
    wgs_full_cohort_n,
    wgs_male_n,
    sum(combined_data$wgs == 1),
    sum(combined_data$wgs == 1 & combined_data$age >= 40)
  ),
  definition = c(
    "Unique participants in the NFCS/MoBa phenotype source file.",
    "Unique male participants in the NFCS/MoBa phenotype source file.",
    "Unique sample IDs in the configured post-QC genotype file; may include samples outside the NFCS phenotype cohort.",
    "Unique NFCS participant IDs present in the configured, post-QC genotype file.",
    "Male NFCS participant IDs present in the configured, post-QC genotype file.",
    "Unique sample IDs in the PHS-variant genotype file.",
    "Unique NFCS participant IDs present in the PHS-variant genotype file.",
    "Unique sample IDs with a score in the PHS profile file.",
    "Male NFCS participant IDs present in the PHS score file.",
    "Male participants with a PHS score and non-missing age and PCa status, after duplicate-ID resolution.",
    "The preceding group restricted to age 40 years or older; input to the primary model before KING filtering.",
    "Unique NFCS phenotype-source participants flagged as having WGS.",
    "Unique male NFCS participants flagged as having WGS.",
    "Male participants flagged for WGS with a PHS score and non-missing age and PCa status.",
    "The preceding WGS group restricted to age 40 years or older; input to the WGS analysis."
  )
)
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

summarize_age_vars <- function(df, birth_year_df) {
  df %>%
    left_join(birth_year_df, by = "pid") %>%
    mutate(
      age_of_metastasis     = if_else(metastasis == 1, age_of_metastasis, NA_real_),
      age_lastobs_censored  = if_else(status == 0, age_lastobs, NA_real_),
      follow_up_years       = NA_real_,
      age_first_follow_up   = NA_real_
    ) %>%
    summarize(
      across(
        c(birth_year, age_first_follow_up, age_diagnosis, age_of_metastasis,
          age_lastobs_censored, follow_up_years, phs),
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
}

proportion_summary <- function(df) {
  n_total <- nrow(df)
  tibble(
    variable = c("pct_pca", "pct_aggressive", "pct_censored", "pct_death_before_dx"),
    n = c(
      sum(df$status == 1),
      sum(df$metastasis == 1),
      sum(df$status == 0),
      sum(df$death == 1 & df$status == 0)
    )
  ) %>%
    mutate(pct = 100 * n / n_total) %>%
    add_row(variable = "n_total", n = n_total, pct = NA_real_, .before = 1) %>%
    mutate(median = NA_real_, q25 = NA_real_, q75 = NA_real_, iqr = NA_real_) %>%
    select(variable, median, q25, q75, iqr, n, pct)
}

analysis_data <- combined_data %>%
  filter(age >= 40)
wgs_analysis_data <- analysis_data %>%
  filter(wgs == 1)

age_summary     <- bind_rows(summarize_age_vars(analysis_data, birth_year_p212),
                             proportion_summary(analysis_data))
age_summary_wgs <- bind_rows(summarize_age_vars(wgs_analysis_data, birth_year_p212),
                             proportion_summary(wgs_analysis_data))

write_tsv(age_summary,     "age_summary.tsv")
write_tsv(age_summary_wgs, "age_summary_wgs.tsv")
record_sample_counts(sample_counts, sample_counts_file, replace_all = TRUE)