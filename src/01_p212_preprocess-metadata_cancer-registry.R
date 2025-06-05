library(tidyverse)
library(readxl)
library(survival)

cancer_reg = read_delim("/ess/p212/data/durable/brca2_cancer_registry/cancer_registry_20180628.tsv", delim = ";")

cancer_reg %>% 
    filter(ICD10_gr == "C61") %>%
    print(width = 1000, n = 10)

cancer_reg %>% 
    filter(ICD10_gr == "C61") %>%
    select(METASTASE) %>%
    print(n = 1000, width = 1000)

metastase_strict = c(1,2,4,6,7)
metastase_loose = c(1,2,4,5,6,7)


cancer_reg %>% 
    filter(KJOENN == "M") %>%
    filter(ICD10_gr == "C61") %>%
    select(regnr, icd10 = ICD10_gr, DOED, DAAR_CANCER, DIAGNOSEAAR, ALDER, METASTASE) %>%
    mutate(regnr = as.character(regnr)) %>%
    mutate(metastase_strict = if_else(METASTASE %in% metastase_strict, 1, 0)) %>%
    mutate(metastase_loose = if_else(METASTASE %in% metastase_loose, 1, 0)) ->
    metastasis

# No observations repeated for the same
metastasis %>%
    group_by(regnr) %>%
    summarize(n = n()) %>%
    arrange(desc(n))

metastasis %>%
    filter(DOED == "Død av cancer" & DAAR_CANCER == "C61")

sum(metastasis$metastase_loose)
sum(metastasis$metastase_strict)

nrow(metastasis)

write_tsv(metastasis, "data/interim/metastasis.tsv")

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
