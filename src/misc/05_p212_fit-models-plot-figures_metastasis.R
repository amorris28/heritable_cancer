library(tidyverse)
library(readxl)
library(survival)

phs = read_tsv("data/processed/PHS601_phs", col_names = c("bfile_id", "phs"))
phenotypes = read_tsv("data/processed/heritable_cancer_phenotypes.tsv")


xl = read_excel("/ess/p212/home/p212-andrewmo/FW_Poly_tables/EDTA_REGNR_06.05.2022.xlsx")
xl %>%
    mutate(bfile_id = paste0("Eivind_", EDTANR)) %>%
    select(bfile_id, regnr = PTN_REGNR) ->
    id_key

metastasis %>%
    left_join(id_key, by = "regnr") ->
    metastasis

phenotypes %>%
    full_join(phs) %>%
    full_join(metastasis) ->
    full_data

full_data %>%
    filter(age != ALDER)

full_data %>%
    filter(!is.na(bfile_id)) ->
    full_data

phs_fit = coxph(Surv(age, status) ~ phs, full_data)
summary(phs_fit)

meta_fit = coxph(Surv(age, metastase_strict) ~ phs, full_data)
summary(meta_fit)

meta_fit = coxph(Surv(age, metastase_loose) ~ phs, full_data)
summary(meta_fit)
