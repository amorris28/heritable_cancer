args <- commandArgs(trailingOnly = TRUE)

model_file <- args[1] # PHS file
inverse <- args[2]
genotype_bfile <- args[3]
phs_bfile <- args[4]
sample_counts_file <- args[5]

# initialize packages
library(tidyverse)
library(lubridate)
library(readxl)

# constants
source("src/config.R")
source("src/functions.R")

# import data
avt_data <- read_delim(
  "/ess/p2521/data/durable/HuskMent/phenotypes/NPRdata/Tilleggsbestilling_Aug_2024/AVT.csv",
  delim = ";"
)
som_data <- read_delim(
  "/ess/p2521/data/durable/HuskMent/phenotypes/NPRdata/Tilleggsbestilling_Aug_2024/SOM.csv",
  delim = ";"
)
phv_data <- read_delim(
  "/ess/p2521/data/durable/HuskMent/phenotypes/NPRdata/Tilleggsbestilling_Aug_2024/PHV.csv",
  delim = ";"
)
death_data <- read_excel(
  "/ess/p2521/data/durable/HuskMent/phenotypes/DAAR_P210891_KOBLET.xlsx"
)
kreft_data <- read_delim(
  "/ess/p2521/data/durable/HuskMent/phenotypes/Kreftregisteret/Utlevert_Haavik/Utlevert_Haavik.csv",
  delim = ";"
)
age_dob_data <- read_delim(
  "/ess/p2521/data/durable/HuskMent/phenotypes/data_2023-09-28/HUSKMENT_data_withkeys_2023-09-28.csv",
  delim = ","
)
birth_month_husk2 <- read_delim(
  "/ess/p2521/data/durable/HuskMent/phenotypes/HUSK_questionairre_data/WORK_QUERY_FOR_HUSK2_FMND_AAR.csv",
  delim = ","
)

# Check last death dates
# no deaths after 12.2020
death_data %>% arrange(desc(DAAR), desc(DMND))

# IDs for each wave
husk2ids <- birth_month_husk2$HUSKMentID
husk1ids <- age_dob_data %>% 
              filter(!HUSKMentID %in% birth_month_husk2$HUSKMentID) %>%
              pull(HUSKMentID)

## Step 1: All IDs
id_key <-
  age_dob_data %>%
  select(
    huskmentid = HUSKMentID,
    bfile_id = decode_subjectid
  ) %>%
  resolve_duplicate_ids(c("huskmentid", "bfile_id"), "p2521 id_key")
sum(!unique(avt_data$lopenr) %in% id_key$huskmentid)
sum(!unique(som_data$lopenr) %in% id_key$huskmentid)
sum(!unique(phv_data$lopenr) %in% id_key$huskmentid)
sum(!unique(death_data$HUSKMENTID) %in% id_key$huskmentid)
sum(!unique(kreft_data$huskmentid) %in% id_key$huskmentid)
sum(!unique(birth_month_husk2$HUSKMentID) %in% id_key$huskmentid)
# It doesn't matter that these are missing because if I don't have DECODE IDs
# for them then I have no genotype data (or no way of connecting genotype data
# to them)

## Step 2: male IDs
# It looks like only cancer registry has gender data
kreft_male_huskmentids <-
  kreft_data %>%
  select(huskmentid, kjoenn = KJOENN) %>%
  filter(kjoenn == "M") %>%
  pull(huskmentid)

male_huskmentids <-
  age_dob_data %>%
  select(
    huskmentid = HUSKMentID,
    sex1 = Gender,
    sex2 = kjona,
    sex3 = kjonn
  ) %>%
  mutate(sex1 = case_match(sex1, "M" ~ 1, "F" ~ 2)) %>%
  mutate(sex = coalesce(sex1, sex2, sex3)) %>%
  select(huskmentid, sex) %>%
  filter(sex == 1) %>%
  pull(huskmentid)

all_male_huskmentids <- unique(c(kreft_male_huskmentids, male_huskmentids))

#########################################################################
## Step 2: Age last obs
### Step 2a: birthdate

# ID DOB
birth_year <-
  age_dob_data %>%
  select(huskmentid = HUSKMentID, faar1, faar) %>%
  mutate(faar3 = if_else(is.na(faar1), faar, faar1)) %>%
  select(-faar, -faar1) %>%
  rename(faar1 = faar3)

birth_month <-
  birth_month_husk2 %>%
  rename(
    huskmentid = HUSKMentID,
    fmnd = invitasjonsregister_fmnd,
    faar2 = invitasjonsregister_faar
  )

dob_combined <-
  birth_year %>%
  full_join(birth_month, by = "huskmentid")

# 10k subjects are missing from birth_month_husk2 and therefore have no birth months
# but they do have IDs from decode so they have been genotype
fodselsdato <-
  dob_combined %>%
  mutate(faar = coalesce(faar2, faar1)) %>%
  select(huskmentid, faar, fmnd) %>%
  mutate(faar = faar + 1900) %>%
  mutate(fmnd = if_else(is.na(fmnd), 1, fmnd)) %>%
  mutate(fdato = ym(paste(faar, fmnd))) %>%
  select(huskmentid, fdato)

## 2b: Date last obs
### registry data
avt_data_proc <-
  avt_data %>%
  select(huskmentid = lopenr, dato = utdato)
som_data_proc <-
  som_data %>%
  select(huskmentid = lopenr, dato = utdato)
phv_data_proc <-
  phv_data %>%
  select(huskmentid = lopenr, dato = utdato)
npr_dates <-
  avt_data_proc %>%
  rbind(som_data_proc) %>%
  rbind(phv_data_proc)

### death data
death_dates <-
  death_data %>%
  select(huskmentid = HUSKMENTID, daar = DAAR, dmnd = DMND) %>%
  mutate(dato = ym(paste(daar, dmnd, sep = "-"))) %>%
  select(huskmentid, dato)

### prostate cancer death data
pca_death_data <-
  death_data %>%
  select(huskmentid = HUSKMENTID, daar = DAAR, dmnd = DMND, alle_koder = ALLE_KODER) %>%
  filter(grepl(icd10pca, alle_koder, fixed = TRUE)) %>%
  mutate(pca_death = 1L) %>%
  mutate(pca_death_date = ym(paste(daar, dmnd, sep = "-"))) %>%
  select(huskmentid, pca_death, pca_death_date)

### cancer data
kreft_dates <-
  kreft_data %>%
  select(huskmentid, diagnosedato = DIAGNOSEDATO) %>%
  mutate(dato = dmy(diagnosedato)) %>%
  select(huskmentid, dato)

censor_date <-
  npr_dates %>%
  rbind(death_dates) %>%
  rbind(kreft_dates) %>%
  group_by(huskmentid) %>%
  summarize(dato = max(dato)) %>%
  rename(date_lastobs = dato)

## 2c. Age last obs
age_lastobs <-
  fodselsdato %>%
  full_join(censor_date, by = "huskmentid") %>%
  mutate(age_lastobs = interval(fdato, date_lastobs) / years(1))

#########################################################################
## Step 3: ID cases (and diagnosis date)
# Process avt_data
avt_data_cases <-
  avt_data %>%
  select(huskmentid = lopenr, icd10, diagnose_dato = inndato) %>%
  mutate(case = grepl(icd10pca, icd10, fixed = TRUE)) %>%
  select(huskmentid, diagnose_dato, case)

# Process som_data
som_data_cases <-
  som_data %>%
  select(huskmentid = lopenr, icd10, diagnose_dato = inndato) %>%
  mutate(case = grepl(icd10pca, icd10, fixed = TRUE)) %>%
  select(huskmentid, diagnose_dato, case)

# Process phv_data
phv_data_cases <-
  phv_data %>%
  select(huskmentid = lopenr, icd10, diagnose_dato = inndato) %>%
  mutate(case = grepl(icd10pca, icd10, fixed = TRUE)) %>%
  select(huskmentid, diagnose_dato, case)

# Process death data
death_data_cases <-
  death_data %>%
  select(
    huskmentid = HUSKMENTID, daar = DAAR, dmnd = DMND, alle_koder = ALLE_KODER
  ) %>%
  mutate(diagnose_dato = ym(paste(daar, dmnd, sep = "-"))) %>%
  mutate(case = grepl(icd10pca, alle_koder, fixed = TRUE)) %>%
  select(huskmentid, diagnose_dato, case)

# Process kreft data
kreft_data_cases <-
  kreft_data %>%
  select(huskmentid, icd10 = ICD10, diagnose_dato = DIAGNOSEDATO) %>%
  mutate(case = grepl(icd10pca, icd10, fixed = TRUE)) %>%
  mutate(diagnose_dato = dmy(diagnose_dato)) %>%
  select(huskmentid, diagnose_dato, case)

diagnose_date <-
  avt_data_cases %>%
  rbind(som_data_cases) %>%
  rbind(phv_data_cases) %>%
  rbind(death_data_cases) %>%
  rbind(kreft_data_cases) %>%
  filter(case) %>%
  group_by(huskmentid) %>%
  mutate(first_diagnose_dato = min(diagnose_dato)) %>%
  ungroup() %>%
  filter(first_diagnose_dato == diagnose_dato) %>%
  select(-first_diagnose_dato)

phenotype_data <-
  age_lastobs %>%
  full_join(diagnose_date, by = "huskmentid") %>%
  mutate(status = if_else(is.na(case), 0, 1)) %>%
  mutate(age_diagnosis = interval(fdato, diagnose_dato) / years(1)) %>%
  mutate(age = if_else(status == 1, age_diagnosis, age_lastobs))

#########################################################################
## Step 4: metastasis and age of metastasis
metastasis_data <-
  kreft_data %>%
  select(huskmentid, icd10 = ICD10, diagnose_dato = DIAGNOSEDATO, metastasis_code = METASTASE) %>%
  mutate(case = grepl(icd10pca, icd10, fixed = TRUE)) %>%
  mutate(metastasis_date = dmy(diagnose_dato)) %>%
  mutate(metastasis = metastasis_code %in% metastase_strict) %>%
  filter(case & metastasis) %>%
  select(huskmentid, metastasis_code, metastasis_date, metastasis)

phenotype_data <-
  phenotype_data %>%
  full_join(metastasis_data, by = "huskmentid") %>%
  mutate(age_of_metastasis = interval(fdato, diagnose_dato) / years(1)) %>%
  mutate(metastasis = if_else(is.na(metastasis), 0, metastasis)) %>%
  mutate(
    age_of_metastasis = if_else(metastasis == 1, age_of_metastasis,
      if_else(is.na(age_lastobs), age, age_lastobs)
    )
  ) %>%
  mutate(conversion_to_metastasis_years = age_of_metastasis - age_diagnosis)

# PCa-specific mortality
phenotype_data <-
  phenotype_data %>%
  left_join(pca_death_data, by = "huskmentid") %>%
  mutate(pca_death = if_else(is.na(pca_death), 0L, pca_death)) %>%
  mutate(
    age_at_pca_death = if_else(
      pca_death == 1L,
      interval(fdato, pca_death_date) / years(1),
      age_lastobs
    )
  )


#########################################################################
## Step 5: date of first follow-up
date_first_follow_up <-
  age_dob_data %>%
  select(huskmentid = HUSKMentID, usdag, usmnd, usaar) %>%
  mutate(date_first_follow_up = ymd(paste(usaar, usmnd, usdag, sep = "-"))) %>%
  group_by(huskmentid) %>%
  summarize(date_first_follow_up = min(date_first_follow_up, na.rm = TRUE))

phenotype_data <-
  phenotype_data %>%
  full_join(date_first_follow_up, by = "huskmentid")


# Filter for male subjects
phenotype_data <-
  phenotype_data %>%
  full_join(id_key, by = "huskmentid") %>%
  filter(huskmentid %in% all_male_huskmentids)

# Add PHS
phs_data <- import_phs_data(model_file, inverse)

combined_data <- process_combined_data(phenotype_data, phs_data)

combined_data <- resolve_duplicate_ids(combined_data, c("huskmentid", "bfile_id"), "p2521 phenotypes")

genotype_ids <- read_fam_iids(genotype_bfile)
phs_bfile_ids <- read_fam_iids(phs_bfile)
male_bfile_ids <- unique(as.character(phenotype_data$bfile_id[!is.na(phenotype_data$bfile_id)]))
sample_counts <- sample_count_rows(
  "p2521",
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
    "male_with_phs_and_complete_phenotype_age40plus"
  ),
  n = c(
    dplyr::n_distinct(age_dob_data$HUSKMentID, na.rm = TRUE),
    length(unique(stats::na.omit(all_male_huskmentids))),
    length(genotype_ids),
    length(intersect(unique(as.character(id_key$bfile_id)), genotype_ids)),
    length(intersect(male_bfile_ids, genotype_ids)),
    length(phs_bfile_ids),
    length(intersect(unique(as.character(id_key$bfile_id)), phs_bfile_ids)),
    dplyr::n_distinct(phs_data$bfile_id, na.rm = TRUE),
    length(intersect(male_bfile_ids, unique(as.character(phs_data$bfile_id)))),
    nrow(combined_data),
    sum(combined_data$age >= 40)
  ),
  definition = c(
    "Unique participants in the HUSK age/DOB source file.",
    "Unique participants classified as male from HUSK phenotype/registry sources.",
    "Unique sample IDs in the configured post-QC genotype file; includes samples outside the HUSK phenotype cohort if present.",
    "Unique HUSK participant IDs present in the configured, post-QC genotype file.",
    "Male HUSK participant IDs present in the configured, post-QC genotype file.",
    "Unique sample IDs in the PHS-variant genotype file.",
    "Unique HUSK participant IDs present in the PHS-variant genotype file.",
    "Unique sample IDs with a score in the PHS profile file.",
    "Male HUSK participant IDs present in the PHS score file.",
    "Male participants with a PHS score and non-missing age and PCa status, after duplicate-ID resolution.",
    "The preceding group restricted to age 40 years or older; input to the primary model before KING filtering."
  )
)
# Save output
write_tsv(combined_data, "phenotypes.tsv")

#########################################################################
## Summary statistics: median and IQR
descriptive_data <- combined_data %>%
  filter(age >= 40)

age_summary <- descriptive_data %>%
  mutate(
    birth_year = year(fdato),
    age_first_follow_up = interval(fdato, date_first_follow_up) / years(1),
    age_of_metastasis = if_else(metastasis == 1, age_of_metastasis, NA_real_),
    age_lastobs_censored = if_else(status == 0, age_lastobs, NA_real_),
    # age_lastobs is last observed age for everyone in the registry, cases included
    end_of_follow_up = if_else(status == 0, age_lastobs, age_diagnosis),
    follow_up_years = end_of_follow_up - age_first_follow_up
  ) %>%
  summarize(
    across(
      c(birth_year, age_first_follow_up, age_diagnosis, age_of_metastasis,
        age_lastobs_censored, follow_up_years),
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

# PHS summary from the same age-eligible, genotyped subjects
phs_summary <- descriptive_data %>%
  summarize(
    median = median(phs, na.rm = TRUE),
    q25    = quantile(phs, 0.25, na.rm = TRUE),
    q75    = quantile(phs, 0.75, na.rm = TRUE)
  ) %>%
  mutate(variable = "phs", iqr = q75 - q25) %>%
  select(variable, median, q25, q75, iqr)

# died_before_dx: censored and appears in death register
phenotype_data_dx <- descriptive_data %>%
  mutate(died_before_dx = (huskmentid %in% death_dates$huskmentid) & (status == 0))

n_total <- nrow(phenotype_data_dx)
proportion_summary_p2521 <- tibble(
  variable = c("pct_pca", "pct_aggressive", "pct_censored", "pct_death_before_dx"),
  n = c(
    sum(phenotype_data_dx$status == 1),
    sum(phenotype_data_dx$metastasis == 1),
    sum(phenotype_data_dx$status == 0),
    sum(phenotype_data_dx$died_before_dx)
  )
) %>%
  mutate(pct = 100 * n / n_total) %>%
  add_row(variable = "n_total", n = n_total, pct = NA_real_, .before = 1) %>%
  mutate(median = NA_real_, q25 = NA_real_, q75 = NA_real_, iqr = NA_real_) %>%
  select(variable, median, q25, q75, iqr, n, pct)

age_summary <- bind_rows(age_summary, phs_summary, proportion_summary_p2521)

write_tsv(age_summary, "age_summary.tsv")
record_sample_counts(sample_counts, sample_counts_file, replace_all = TRUE)
