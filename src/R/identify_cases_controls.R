library(tidyverse)
library(readxl)

allData = read_excel('FW__Poly_tables/poly_final_220425_estimated_lastobs.xls', sheet = 1)
caseData = read_excel('FW__Poly_tables/poly_final_220425_estimated_lastobs.xls', sheet = 2)

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

write_lines(prostateIDs, 'case_pids.txt')
write_lines(controlIDs, 'control_pids.txt')

xl = read_excel("FW__Poly_tables/EDTA_REGNR_06.05.2022.xlsx")
case = read_csv("case_pids.txt", col_names = "PTN_REGNR", col_types = "c")
control = read_csv("control_pids.txt", col_names = "PTN_REGNR", col_types = "c")
bind_rows(case, control) %>% 
  left_join(xl, by = "PTN_REGNR") %>% 
  mutate(bfile_id = paste("Eivind", EDTANR, sep = "_")) ->
  samples
write_tsv(samples, "sample_key.txt")

samples %>% 
  mutate(bfile_id2 = bfile_id) %>% 
  select(bfile_id, bfile_id2) ->
  bfile_id_samples
write_tsv(bfile_id_samples, "bfile_id_samples.txt", col_names = FALSE)
