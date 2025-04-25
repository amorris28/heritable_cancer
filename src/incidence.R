library(tidyverse)
library(readxl)
library(survival)

source("src/00_config.R")

# p212 heritable cancer data
# allData = read_excel("/ess/p212/home/p212-andrewmo/FW_Poly_tables/poly_final_220425_estimated_lastobs.xls", sheet = 1)
# caseData = read_excel("/ess/p212/home/p212-andrewmo/FW_Poly_tables/poly_final_220425_estimated_lastobs.xls", sheet = 2)
# rareVariants = read_excel("/ess/p212/home/p212-andrewmo/FW_Poly_tables/poly_final_220425_estimated_lastobs.xls", sheet = 3)
# fam_file = read_table("data/raw/moba_imputed_all.fam", col_names = FALSE) %>% as_tibble()

# kreftregisteret statistics
incidence = read_delim("data/external/Incidence_20250320_150713.csv", delim = ";", col_select = 1:11, n_max = 9)
incidence = read_delim("data/external/incidence_norway_agegroups.csv", delim = ";", col_select = 1:11, n_max = 16)

# p212 phenotype data
p212_data = read_tsv(paste0(PROCESSED_DATA, "heritable_cancer_phenotypes.tsv"))
p212_phs = read_tsv(paste0(PROCESSED_DATA, "p212_PHS601_phs"), col_names = c("bfile_id", "phs"))

p212_data %>%
    full_join(phs_p212, by = "bfile_id") %>%
    filter(!is.na(phs)) %>%
    filter(!is.na(age)) ->
    full_data

# p2521 data
p2521_data = read_tsv(paste0(PROCESSED_DATA, "husk_phenotypes.tsv"))

###############################################################################
# plot incidence in Norway
###############################################################################
# According to SSB, the average life expectancy for males in Norway at age 85 is 6 years
# So I will use 85-91 as the top interval.
# Then, I calculate the median age for each interval
incidence$Alder[incidence$Alder == "85+"] = "85-91"
incidence %>% 
    rowwise() %>%
    mutate(median_age = median(seq(as.numeric(substr(Alder, 1, 2)), as.numeric(substr(Alder, 4, 5))))  ) ->
    incidence
# calculate natural log of incidence rate
incidence %>%
    mutate(log_incidencerate = log(Aldersspesifikkrate)) ->
    incidence
# First cases are in age group 35-40
incidence %>%
    filter(median_age >= 35) %>%
    mutate(median_age_minus35 = median_age - min(median_age)) ->
    incidence_over35

# Visualize age-specific incidence rate on log scale
# with polynomial fit
ggplot(incidence_over35, aes(x = median_age_minus35, y = log_incidencerate)) +
    geom_point() +
    stat_smooth(method='lm', formula = y ~ poly(x,2)) +
    theme_classic() +
    scale_x_continuous(breaks = c(5, 15, 25, 35, 45), labels = c(40, 50, 60, 70, 80)) +
    labs(x = "Age", y = "log(Age-specific Incidence Rate)") ->
    p
ggsave("reports/figures/log_incidence_rate.png", p, width = 4.5, height = 3.5, dpi = 300)

# Compare polynomial fit between fitting a linear model to the log-transformed
# data or fit a non-linear model to the original data
lm_fit = lm(log_incidencerate ~ poly(median_age_minus35, 2, raw = TRUE), data = incidence_over35)
# plot(incidence_over35$agegroup_minus35, resid(lm_fit))
# abline(h = 0, lty = 2)

nls_fit <- nls(Aldersspesifikkrate ~ b0 * exp(b1 * median_age_minus35 + b2 * median_age_minus35 ^ 2), 
               data = incidence_over35,
               start = list(b0 = 0.5, b1 = 0.4, b2 = -0.005))
coef(nls_fit)
# plot(incidence_over35$agegroup_minus35, resid(nls_fit))
# abline(h = 0, lty = 2)

# The linear model is more uniformly distributed residuals than the 
# non-linear model, likely because variance increases with the 


B0 = lm_fit$coefficients[1]
B1 = lm_fit$coefficients[2]
B2 = lm_fit$coefficients[3]

x = incidence_over35$median_age_minus35

pca_incidence_poly = function(x) {
    B0 + B1 * x + B2 * x ^ 2
}

new_data = incidence_over35 %>% select(median_age_minus35)

y_equation = pca_incidence_poly(x)
y_predict = predict(lm_fit, new_data)

y = y_equation

ggplot(mapping = aes(x, y)) + 
    geom_smooth() + #loess
    geom_point(aes(x = median_age_minus35, y = log_incidencerate), data = incidence_over35) +
    theme_classic() +
    scale_x_continuous(breaks = c(5, 15, 25, 35, 45), labels = c(40, 50, 60, 70, 80)) +
    labs(x = "Age", y = "log(Age-specific Incidence Rate)") ->
    p
ggsave("reports/figures/log_incidence_rate_lmfit.png", p, width = 4.5, height = 3.5, dpi = 300)

###########################
# if I take the log of y and fit a second-order polynomial model I get this:
# log(y) = log(b0) + b1 * x + b2 * x ^ 2 
# To exponentiate it this back to linear space, we take the exp of both sizes
# the exp(log()) cancels out for y and b0 and we are left with:
# y = b0 + e ^ {b1 * x + b2 * x ^ 2}

b0 = exp(B0)
b1 = B1
b2 = B2

pca_incidence_poly_linearized = function(x) {
    b0 + exp(b1 * x + b2 * x ^ 2)
}

x_new = seq(from = 0, to = 50, by = 0.01)
y_new = pca_incidence_poly_linearized(x_new)

# Exponentiate the model prediction to return it to the original units
ggplot() +
    geom_point(aes(x = median_age_minus35, y = Aldersspesifikkrate), incidence_over35) +
    geom_smooth(aes(x_new, y_new)) +
    scale_x_continuous(breaks = c(5, 15, 25, 35, 45), labels = c(40, 50, 60, 70, 80)) +
    labs(x = "Age", y = "Age-specific Incidence Rate (per 100,000)") +
    theme_classic() ->
    p
ggsave("reports/figures/log_incidence_rate_lmfit_exp.png", p, width = 4.5, height = 3.5, dpi = 300)

# The model doesn't match up very well with the original in this case
# Try the nls model:
new_data = data.frame(median_age_minus35 = x_new)
y_predict = predict(nls_fit, new_data)
y = y_predict
ggplot()  +
    geom_point(aes(x = median_age_minus35, y = Aldersspesifikkrate), incidence_over35) +
    geom_smooth(aes(x_new, y)) +
    scale_x_continuous(breaks = c(5, 15, 25, 35, 45), labels = c(40, 50, 60, 70, 80)) +
    labs(x = "Age", y = "Age-specific Incidence Rate (per 100,000)") +
    theme_classic() ->
    p1
ggsave("reports/figures/norway_incidence.png", p1, width = 4.5, height = 3.5, dpi = 300)

# The non-linear model seems to fit better

fit = coxph(Surv(age, status) ~ 1, data = p212_data)
p212_incidence = basehaz(fit)
p212_incidence %>%
    mutate(cum_inc_per_100k = hazard * 100000) %>%
    mutate(age_minus_35 = time - 37) ->
    p212_incidence
p1 = ggplot(p212_incidence, aes(x = time, y = cum_inc_per_100k)) +
        geom_point()
ggsave("reports/figures/p212_incidence.png", p1, width = 4.5, height = 3.5, dpi = 300)


p3 = ggplot() +
    geom_smooth(aes(x_new, y)) +
    geom_point(aes(age_minus_35, cum_inc_per_100k), data = p212_incidence) +
    scale_x_continuous(breaks = c(5, 15, 25, 35, 45), labels = c(40, 50, 60, 70, 80)) +
    labs(x = "Age", y = "Age-specific Incidence Rate (per 100,000)") +
    theme_classic()
ggsave("reports/figures/combined_incidence_p212.png", p3, width = 4.5, height = 3.5, dpi = 300)

#################################################################################################
# Adjust curves by phs groups

fit = coxph(Surv(age, status) ~ phs, data = full_data)
summary(fit)
new_data = data.frame(phs = quantile(full_data$phs, 0.80))
percentile_80 = predict(fit, new_data)


# prev_norge = read_delim("data/external/Prevalence_20250320_163654.csv", delim = ";")
# prev_rogaland = read_delim("data/external/Prevalence_20250320_163603_Rogaland.csv", delim = ";")

# prev_rogaland$cases = as.numeric(gsub(" ", "", prev_rogaland$cases))
# prev_norge$cases = as.numeric(gsub(" ", "", prev_norge$cases))

# Try to plot cumulatice incidence with p212 data using coxph 

fit = coxph(Surv(age, status) ~ 1, data = p212_data)
p212_incidence = basehaz(fit)
p212_incidence %>%
    mutate(cum_inc_per_100k = hazard * 100000) ->
    p212_incidence
p1 = ggplot(p212_incidence, aes(x = time, y = cum_inc_per_100k)) +
        geom_point()
ggsave("reports/figures/p212_incidence.png", p1)

# Try to plot the p2521 data

p2521_fit = coxph(Surv(age, status) ~ 1, data = p2521_data)
summary(p2521_fit)
p2521_survfit = survfit(p2521_fit)
summary(p2521_survfit)

p2521_incidence = basehaz(p2521_fit)
p2521_incidence %>%
    mutate(cum_inc_per_100k = hazard * 100000) ->
    p2521_incidence
p1 = ggplot(p2521_incidence, aes(x = time, y = cum_inc_per_100k)) +
        geom_point()
ggsave("reports/figures/p2521_incidence.png", p1)

# Plot cumulative incidence based on kreftregisteret data

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
p2 = ggplot(spline_int, aes(x, y)) +
           geom_line() +
        geom_point(aes(x = Age_group, y = ASR_cumsum), data = new_incidence, color = "red")

ggsave("reports/figures/kreftregisteret_incidence.png", p2)


# Incidence combined

p3 = ggplot() +
        geom_line(aes(x, y), data = spline_int) +
        geom_point(aes(time, cum_inc_per_100k), data = p212_incidence)
ggsave("reports/figures/combined_incidence_p212.png", p3)
# Incidence combined

p4 = ggplot() +
        geom_line(aes(x, y), data = spline_int) +
        geom_point(aes(time, cum_inc_per_100k), data = p2521_incidence)
ggsave("reports/figures/combined_incidence_p2521.png", p4)

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