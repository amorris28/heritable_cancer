library(tidyverse)

incidence = read_delim("data/external/Incidence_20250320_150713.csv", delim = ";", col_select = 1:11, n_max = 9)
prev_norge = read_delim("data/external/Prevalence_20250320_163654.csv", delim = ";")
prev_rogaland = read_delim("data/external/Prevalence_20250320_163603_Rogaland.csv", delim = ";")

prev_rogaland$cases = as.numeric(gsub(" ", "", prev_rogaland$cases))
prev_norge$cases = as.numeric(gsub(" ", "", prev_norge$cases))

prev_norge %>% 
    mutate(cases = as.numeric(gsub(" ", "", cases))) %>% 
    mutate(pop = as.numeric(gsub(" ", "", pop))) %>%
    mutate(new_cases = cases - lag(cases))
    mutate(aldersspesifikkrte = cases / pop * 100000) %>%
    mutate(cum_case_rate = cumsum(case_rate)) %>% 
    mutate(age_group = as.numeric(substr(Age, 1, 2))) ->
    prev_norge

prev_rogaland %>% 
    mutate(cases = as.numeric(gsub(" ", "", cases))) %>% 
    mutate(pop = as.numeric(gsub(" ", "", pop))) %>% 
    mutate(case_rate = cases / pop) %>%
    mutate(cum_case_rate = cumsum(case_rate)) %>% 
    mutate(age_group = as.numeric(substr(Age, 1, 2))) ->
    prev_rogaland

ggplot(prev_norge, aes(x = age_group, y = cum_case_rate)) +
    geom_point()

# Incidence, the Aldersspesifikkrate is the rate of new diagnoses per
# 100 000 person-years. (i.e., Aldersspesifikkrate / 100000 * Personår = Tilfeller
incidence %>% 
    mutate(csum_tilfeller = cumsum(Tilfeller)) %>% 
    mutate(csum_Aldersspesifikkrater = cumsum(Aldersspesifikkrate) / 100000) ->
    incidence

ggplot(incidence, aes(x = Alder, y = csum_Aldersspesifikkrater)) +
    geom_point()
