library(tidyverse)
library(readxl)

xl = read_excel("FW__Poly_tables/EDTA_REGNR_06.05.2022.xlsx")
case = read_csv("case_pids.txt", col_names = "PTN_REGNR", col_types = "c")
control = read_csv("control_pids.txt", col_names = "PTN_REGNR", col_types = "c")
bind_rows(case, control) %>% 
  inner_join(xl, by = "PTN_REGNR") %>% 
  mutate(Eivind = paste("Eivind", EDTANR, sep = "_")) %>% 
  mutate(Eivind2 = Eivind) ->
  samples
write_tsv(samples, "sample_key.txt", col_names = FALSE)

samples %>% 
  select(Eivind, Eivind2) ->
  eivind_samples
write_tsv(eivind_samples, "eivind_samples.txt", col_names = FALSE)