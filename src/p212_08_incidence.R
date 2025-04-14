library(tidyverse)
library(readxl)

source("src/00_config.R")

# p212 heritable cancer data
allData = read_excel("/ess/p212/home/p212-andrewmo/FW_Poly_tables/poly_final_220425_estimated_lastobs.xls", sheet = 1)
caseData = read_excel("/ess/p212/home/p212-andrewmo/FW_Poly_tables/poly_final_220425_estimated_lastobs.xls", sheet = 2)
rareVariants = read_excel("/ess/p212/home/p212-andrewmo/FW_Poly_tables/poly_final_220425_estimated_lastobs.xls", sheet = 3)
fam_file = read_table("data/raw/moba_imputed_all.fam", col_names = FALSE) %>% as_tibble()

# kreftregisteret statistics
incidence = read_delim("data/external/Incidence_20250320_150713.csv", delim = ";", col_select = 1:11, n_max = 9)
# prev_norge = read_delim("data/external/Prevalence_20250320_163654.csv", delim = ";")
# prev_rogaland = read_delim("data/external/Prevalence_20250320_163603_Rogaland.csv", delim = ";")

# prev_rogaland$cases = as.numeric(gsub(" ", "", prev_rogaland$cases))
# prev_norge$cases = as.numeric(gsub(" ", "", prev_norge$cases))

# Try to run survival analysis on kreftregisteret data

incidence %>%
    mutate(ASR_cumsum = cumsum(Aldersspesifikkrate)) %>%
    mutate(Age_group = as.numeric(substr(Alder, 0, 2))) ->
    new_incidence
t = new_incidence$Age_group
y = new_incidence$ASR_cumsum
    
add_nls <- nls(y ~ a*exp(r*t), 
               start = list(a = 0.5, r = 0.2))
nls_coef = coef(add_nls)
y = nls_coef[1]*exp(nls_coef[2]*t)
spline_int <- as.data.frame(spline(t, y))


########
# At least a plot of the cumulative incidence from cancer registry
# With age-specific rate per 100,000 on the y-axis
ggplot(spline_int, aes(x, y)) +
           geom_line() +
        geom_point(aes(x = Age_group, y = ASR_cumsum), data = new_incidence, color = "red")



#     # mutate(new_cases = cases - lag(cases))
# prev_norge %>% 
#     mutate(cases = as.numeric(gsub(" ", "", cases))) %>% 
#     mutate(pop = as.numeric(gsub(" ", "", pop))) %>%
#     mutate(aldersspesifikkrate = cases / pop * 100000) %>%
#     mutate(cum_case_rate = cumsum(aldersspesifikkrate)) %>% 
#     mutate(age_group = as.numeric(substr(Age, 1, 2))) %>%
#     print(width = 100) ->
#     prev_norge

# prev_rogaland %>% 
#     mutate(cases = as.numeric(gsub(" ", "", cases))) %>% 
#     mutate(pop = as.numeric(gsub(" ", "", pop))) %>% 
#     mutate(case_rate = cases / pop) %>%
#     mutate(cum_case_rate = cumsum(case_rate)) %>% 
#     mutate(age_group = as.numeric(substr(Age, 1, 2))) ->
#     prev_rogaland

# ggplot(prev_norge, aes(x = age_group, y = cum_case_rate)) +
#     geom_point() ->
#     p
# ggsave("prevalence_test.png", p)

# Incidence, the Aldersspesifikkrate is the rate of new diagnoses per
# 100 000 person-years. (i.e., Aldersspesifikkrate / 100000 * Personår = Tilfeller
incidence %>% 
    select(Alder, Personår:Tilfeller) ->
    incidence
incidence %>%
    mutate(csum_tilfeller = cumsum(Tilfeller)) %>% 
    mutate(csum_Aldersspesifikkrater = cumsum(Aldersspesifikkrate) / 100000) ->
    incidence

incidence %>% 
    mutate(Tilfeller2 = Aldersspesifikkrate / 10 / 10000 * Personår) %>%
    print(width = 10000)

ggplot(incidence, aes(x = Alder, y = csum_Aldersspesifikkrater)) +
    geom_point() ->
    p1
ggsave("incidence_test.png", p1)


####################################################################################
# Reshape p212 data into this format:
# Age group (tiaar), Pop, Cases
allData %>%
    filter(SEX == "M") %>%
    select(pid = PID, year_birth = YEAR_BIRTH, agelastobs = AGELASTOBS) ->
    all_data_cleaned

caseData %>%
    filter(icd9 == "185") %>%
    group_by(pid) %>%
    filter(age_ca == min(age_ca)) %>%
    ungroup() %>%
    select(pid, yr_diagn, icd9, age_ca) ->
    case_data_cleaned
    
all_data_cleaned %>%
    full_join(case_data_cleaned, by = "pid") ->
    full_data

## Check for multiple diagnosis dates for the same individual
# full_data %>% group_by(pid) %>% summarize(n = n()) %>% arrange(desc(n)) %>%
# filter(n == 2) %>% pull(pid) -> repeat_pids

# full_data %>%
#     filter(pid %in% repeat_pids)

full_data %>%
    select(icd9, age_ca, agelastobs) %>%
    mutate(case = if_else(is.na(icd9), 0, 1)) %>%
    mutate(age = if_else(case == 1, age_ca, agelastobs)) %>%
    select(age, case) %>%
    mutate(decade = substr(age, 1, 1)) %>%
    select(decade, case) %>%
    group_by(decade) %>%
    summarize(pop = n(),
              cases = sum(case)) %>%
    mutate(case_rate = cases / pop)

# Incidience, i.e., year of diagnosis 