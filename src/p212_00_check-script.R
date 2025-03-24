library(tidyverse)
library(readxl)
library(genio)

source("src/00_config.R")

allData = read_excel("/tsd/p212/home/p212-andrewmo/FW_Poly_tables/poly_final_220425_estimated_lastobs.xls", sheet = 1)
caseData = read_excel("/tsd/p212/home/p212-andrewmo/FW_Poly_tables/poly_final_220425_estimated_lastobs.xls", sheet = 2)
rareVariants = read_excel("/tsd/p212/home/p212-andrewmo/FW_Poly_tables/poly_final_220425_estimated_lastobs.xls", sheet = 3)
fam_file = read_table("data/raw/moba_imputed_all.fam", col_names = FALSE) %>% as_tibble()


allData %>%
	pull("PID") ->
	all_subject_ids



