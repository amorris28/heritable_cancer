library(tidyverse)
library(readxl)

source("src/00_config.R")

allData = read_excel("/ess/p212/home/p212-andrewmo/FW_Poly_tables/poly_final_220425_estimated_lastobs.xls", sheet = 1)
caseData = read_excel("/ess/p212/home/p212-andrewmo/FW_Poly_tables/poly_final_220425_estimated_lastobs.xls", sheet = 2)
rareVariants = read_excel("/ess/p212/home/p212-andrewmo/FW_Poly_tables/poly_final_220425_estimated_lastobs.xls", sheet = 3)

# Get case IDs and diagnosis date
caseData %>%
    filter(icd9 == 174) %>%
    select(pid, age_diag = age_ca) %>%
    mutate(status = 1) ->
    bca_cases

# Get control IDs and date last obs
allData %>%
    filter(SEX == "F") %>%
    select(pid = PID, age_lob = AGELASTOBS) ->
    female_subjects

# link cases and controls
female_subjects %>%
    full_join(bca_cases, by = "pid") %>%
    mutate(age = if_else(is.na(age_diag), age_lob, age_diag)) %>%
    mutate(status = if_else(is.na(status), 0, status)) %>%
    select(-age_lob, -age_diag) ->
    combined_data

# Pull out BRCA2 carriers
rareVariants %>%
    filter(GENE %in% c("BRCA1", "BRCA2")) %>%
    select(pid = PID) %>%
    mutate(carrier = 1) ->
    BRCA_carriers

combined_data %>%
    full_join(BRCA2_carriers, by = "pid") %>%
    mutate(carrier = if_else(is.na(carrier), 0 , carrier)) ->
    combined_data

combined_data %>%
    select(carrier, status) %>%
    table() ->
    any_BRCA_table

fisher.test(any_BRCA_table)

