library(tidyverse)
library(survival)
library(lubridate)
source("src/00_config.R")
phenotypes = read_tsv(paste0(PROCESSED_DATA, "husk_phenotypes.tsv"))
phses = read_tsv(paste0(PROCESSED_DATA, "p2521_PHS601_phs"), col_names = c("subjectid", "phs"))
death_data = read_tsv(paste0(PROCESSED_DATA, "death_data_proc.tsv"))
diagnosedato = read_tsv(paste0(PROCESSED_DATA, "diagnosedato.tsv"))
id_key = read_tsv(paste0(PROCESSED_DATA, "id_key.tsv"))
sensurdato = read_tsv(paste0(PROCESSED_DATA, "husk_sensurdato.tsv"))

phenotypes %>%
    left_join(phses, by = "subjectid") ->
    combined_data

# Full survival analysis (diagnosis as event)

fit = coxph(Surv(age, status) ~ phs, data = combined_data)
summary(fit)

# Case-only survival analysis (diagnosis as event)
combined_data %>% 
    filter(status == 1) ->
    cases_only
casesfit = coxph(Surv(age, status) ~ phs, data = cases_only)
summary(casesfit)

# Case-only survival anlysis (all-cause mortality within 5 years of diagnosis as event)
# The final structure I want:
# Subjectid, Age at diagnosis, status (1 if dead at 5 years), phs

## Identify death date
death_data %>%
	mutate(dato = ym(paste(daar, dmnd, sep = "-"))) %>%
    select(subjectid, case, dato) %>%
    rename(died_of_pca = case) ->
    death_data_proc
combined_data

diagnosedato %>%
	left_join(id_key, by = "subjectid") %>%
	mutate(subjectid = decode_subjectid) %>%
	select(-decode_subjectid) ->
	diagnosedato

death_data_proc %>%
    rename(death_date = dato) %>%
    filter(!is.na(subjectid)) ->
    death_date
diagnosedato %>%
    rename(diagnose_date = dato) %>%
    select(-status) %>%
    filter(!is.na(subjectid)) ->
    diagnose_date

combined_data %>%
    filter(!is.na(subjectid)) %>%
    left_join(diagnose_date, by = "subjectid") %>%
    left_join(death_date, by = "subjectid") ->
    alldata

alldata %>%
    mutate(diag_dth_interval = interval(diagnose_date, death_date) / years(1)) %>%
    mutate(diag_dth_interval = if_else(is.na(diag_dth_interval), Inf, diag_dth_interval)) %>%
    mutate(dead_at_5 = as.numeric(diag_dth_interval <= 5)) ->
    all_data_acm 
all_data_acm %>% print(width = Inf)
# Censoring date is 12.2020 (last deaths in registry)
sensurdato %>%
    left_join(id_key, by = "subjectid") %>%
	mutate(subjectid = decode_subjectid) %>%
	select(-decode_subjectid) %>%
    rename(censoring_date = dato) ->
    censoring_dates

# comparing the last observed date to the last date in death registry (12/2020)
# some people were last observed 14+ years ago but are still absent from death
# registry. Therefore, I should use last observed unless the death registy date
# is later
censoring_dates %>%
    mutate(new_interval = interval(my("12-2020"), censoring_date) / years(1))

censoring_dates %>%
    mutate(new_date = if_else(censoring_date < my("12-2020"), my("12-2020"), censoring_date)) %>%
    mutate(censoring_date = new_date) %>%
    select(-new_date) ->
    censoring_dates

# censoring_date = my("12-2020")
all_data_acm %>%
    left_join(censoring_dates, by = "subjectid") %>%
    # mutate(censoring_date = censoring_date) %>%
    mutate(diag_cens_interval = interval(diagnose_date, censoring_date) / years(1)) %>%
    mutate(years = if_else(diag_dth_interval <= 5, diag_dth_interval, if_else(diag_cens_interval > 0, diag_cens_interval, 0))) ->
    all_data_acm_plus_censoring
all_data_acm_plus_censoring %>% print(width = Inf, n = 100)

casesfit = coxph(Surv(years, dead_at_5) ~ phs + age, data = all_data_acm_plus_censoring)
summary(casesfit)
casesfit = coxph(Surv(years, dead_at_5) ~ phs, data = all_data_acm_plus_censoring)
summary(casesfit)

# Case-only survival anlysis (death FROM PCA within 5 years of diagnosis as event)
all_data_acm_plus_censoring %>%
    mutate(died_of_pca_at_5 = if_else(dead_at_5 == 1 & died_of_pca, 1, 0)) ->
    all_data_pca_mortality

casesfit = coxph(Surv(years, died_of_pca_at_5) ~ phs + age, data = all_data_pca_mortality)
summary(casesfit)
casesfit = coxph(Surv(years, died_of_pca_at_5) ~ phs, data = all_data_pca_mortality)
summary(casesfit)
# Full survival analysis (metastatic disease as event with metastatic disease
# defined as long-term immunotherapy from drug registry)