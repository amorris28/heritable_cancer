library(tidyverse)
library(readxl)

source("src/00_config.R")

allData = read_excel("/tsd/p212/home/p212-andrewmo/FW__Poly_tables/poly_final_220425_estimated_lastobs.xls", sheet = 1)
caseData = read_excel("/tsd/p212/home/p212-andrewmo/FW__Poly_tables/poly_final_220425_estimated_lastobs.xls", sheet = 2)

# Grab all male IDs
allData %>% 
  filter(SEX == "M") %>% 
  pull(PID) %>% 
  unique() ->
  allMaleIDs

# Grab IDs of prostate cancer cases
caseData %>% 
  filter(pid %in% allMaleIDs) %>% 
  filter(icd9 == 185) %>% 
  pull(pid) %>% 
  unique ->
  prostateIDs

# Grab IDs of all cancer cases
caseData %>% 
  filter(pid %in% allMaleIDs) %>% 
  pull(pid) %>% 
  unique() ->
  caseIDs

# Grab IDs for all cancer-free males
allData %>% 
  filter(PID %in% allMaleIDs) %>% 
  filter(!(PID %in% caseIDs)) %>% 
  pull(PID) %>% 
  unique() ->
  controlIDs

# Double check no cancer cases or prostate cases are in control IDs
try(for (i in caseIDs) (if (i %in% controlIDs) stop("Found case ID in controls.")))
try(for (i in prostateIDs) (if (i %in% controlIDs) stop("Found prostate ID in controls.")))

case = tibble(PTN_REGNR = as.character(prostateIDs))
control = tibble(PTN_REGNR = as.character(controlIDs))

xl = read_excel("/tsd/p212/home/p212-andrewmo/FW__Poly_tables/EDTA_REGNR_06.05.2022.xlsx")

bind_rows(case, control) %>% 
  left_join(xl, by = "PTN_REGNR") %>% 
  mutate(bfile_id = paste("Eivind", EDTANR, sep = "_")) ->
  samples

samples %>% 
  mutate(bfile_id2 = bfile_id) %>% 
  select(bfile_id, bfile_id2) ->
  bfile_id_samples
bfile_id_samples$bfile_id = 0
write_tsv(bfile_id_samples, paste0(INTERIM_DATA, "bfile_id_samples.txt"), col_names = FALSE)

########################################################
# bfile_id, age, status
# Eivind_00001, 56, 0/1
caseData %>%
  filter(icd9 == 185) %>%
  rename(PTN_REGNR = pid) %>%
  rename(age = age_ca) %>%
  mutate(status = 1) %>%
  select(PTN_REGNR, age, status) %>%
  group_by(PTN_REGNR) %>%
  filter(age == min(age)) %>%
  ungroup() -> cases

allData %>%
  filter(PID %in% controlIDs) %>%
  rename(PTN_REGNR = PID, age = AGELASTOBS) %>%
  select(PTN_REGNR, age) %>%
  mutate(status = 0) %>%
  group_by(PTN_REGNR) -> controls

bind_rows(cases, controls) %>%
  mutate(PTN_REGNR = as.character(PTN_REGNR)) ->
  combined_data

samples %>%
  select(bfile_id, PTN_REGNR) %>%
  right_join(combined_data, by = "PTN_REGNR") %>%
  select(-PTN_REGNR) -> phenotypes

write_tsv(phenotypes, paste0(PROCESSED_DATA, "familial_cancer_phenotypes.tsv"))
