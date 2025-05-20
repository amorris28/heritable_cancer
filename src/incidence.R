library(tidyverse)
library(readxl)
library(survival)

source("src/00_config.R")

# cancer registry incidence data
incidence = read_delim("data/external/incidence_norway_agegroups.csv", delim = ";", col_select = 1:11, n_max = 16)

# p212 phenotype data
p212_data = read_tsv(paste0(PROCESSED_DATA, "heritable_cancer_phenotypes.tsv"))
p212_phs = read_tsv(paste0(PROCESSED_DATA, "p212_PHS601_phs"), col_names = c("bfile_id", "phs"))

p212_data %>%
    full_join(p212_phs, by = "bfile_id") %>%
    filter(!is.na(phs)) %>%
    filter(!is.na(age)) ->
    full_data

# p2521 data
p2521_data = read_tsv(paste0(PROCESSED_DATA, "husk_phenotypes.tsv"))
p2521_phs = read_tsv(paste0(PROCESSED_DATA, "p2521_PHS601_phs"), col_names = c("subjectid", "phs"))

p2521_data %>%
    full_join(p2521_phs, by = "subjectid") %>%
    filter(!is.na(phs)) %>%
    filter(!is.na(age)) ->
    full_data

################################################################################
# 1. Plot incidence curve for Norway
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
    mutate(median_age_minus35 = median_age - 35) ->
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
#ggsave("reports/figures/log_incidence_rate.png", p, width = 4.5, height = 3.5, dpi = 300)

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
pca_incidence_poly = function(x) {
    B0 + B1 * x + B2 * x ^ 2
}

x = incidence_over35$median_age_minus35

new_data = incidence_over35 %>% select(median_age_minus35)

y_equation = pca_incidence_poly(x)
y_predict = predict(lm_fit, new_data) # should return the same

y = y_equation

ggplot(mapping = aes(x, y)) + 
    geom_smooth() + #loess
    geom_point(aes(x = median_age_minus35, y = log_incidencerate), data = incidence_over35) +
    theme_classic() +
    scale_x_continuous(breaks = c(5, 15, 25, 35, 45), labels = c(40, 50, 60, 70, 80)) +
    labs(x = "Age", y = "log(Age-specific Incidence Rate)") ->
    p
#ggsave("reports/figures/log_incidence_rate_lmfit.png", p, width = 4.5, height = 3.5, dpi = 300)

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

x_new = seq(from = min(incidence_over35$median_age_minus35), to = max(incidence_over35$median_age_minus35), by = 0.01)
y_new = pca_incidence_poly_linearized(x_new)

# Exponentiate the model prediction to return it to the original units
ggplot() +
    geom_point(aes(x = median_age_minus35, y = Aldersspesifikkrate), incidence_over35) +
    geom_smooth(aes(x_new, y_new)) +
    scale_x_continuous(breaks = c(5, 15, 25, 35, 45), labels = c(40, 50, 60, 70, 80)) +
    labs(x = "Age", y = "Age-specific Incidence Rate (per 100,000)") +
    theme_classic() ->
    p
#ggsave("reports/figures/log_incidence_rate_lmfit_exp.png", p, width = 4.5, height = 3.5, dpi = 300)

# The model doesn't match up very well with the original in this case

# Try the nls model:
x_new = seq(from = min(incidence_over35$median_age_minus35), to = max(incidence_over35$median_age_minus35), by = 0.01)

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
#ggsave("reports/figures/norway_incidence.png", p1, width = 4.5, height = 3.5, dpi = 300)

# The non-linear model seems to fit better

# Predict the incidence rate at the midpoint of each 1-year age group
yearly_data = data.frame(median_age_minus35 = seq(from = min(incidence_over35$median_age_minus35) + 0.5, to = max(incidence_over35$median_age_minus35) + 0.5, by = 1))
yearly_data$y_predict = predict(nls_fit, yearly_data)
ggplot()  +
    geom_point(aes(x = median_age_minus35, y = Aldersspesifikkrate), incidence_over35) +
    geom_smooth(aes(x_new, y_predict)) +
    geom_point(aes(median_age_minus35, y_predict), yearly_data, color = "red") +
    scale_x_continuous(breaks = c(5, 15, 25, 35, 45), labels = c(40, 50, 60, 70, 80)) +
    labs(x = "Age", y = "Age-specific Incidence Rate (per 100,000)") +
    theme_classic() ->
    p1
ggsave("reports/figures/norway_incidence_yearly.png", p1, width = 4.5, height = 3.5, dpi = 300)

# Now, take the cumulative incidence for later comparison to p212 and p2521
yearly_data %>%
    mutate(cumulative_incidence = cumsum(y_predict)) ->
    yearly_data

ggplot(yearly_data, aes(x = median_age_minus35, y = cumulative_incidence)) +
    geom_point() ->
    p

ggsave("reports/figures/norway_cumulative_incidence.png", p, width = 4.5, height = 3.5, dpi = 300)

###################################################################################
# Plot the time-to-event data from the two cohorts next to the cancer registry data

fit = coxph(Surv(age, status) ~ phs, data = full_data)
H = data.frame(basehaz(fit)) # This is the cumulative baseline hazard
H %>%
    mutate(S = exp(-hazard)) %>% # This is the cumulative surival
    mutate(I = 1-S) %>%
    mutate(I_100000 = I * 100000) %>%
    mutate(time_minus_35 = time - 35) -> # This is the cumulative incidence
    H


ggplot() +
    geom_step(aes(time_minus_35, I_100000), H, color = "red") +
    geom_line(aes(x = median_age_minus35, y = cumulative_incidence), yearly_data) +
    scale_x_continuous(breaks = c(5, 15, 25, 35, 45, 55), labels = c(40, 50, 60, 70, 80, 90), limits = c(0, 55)) +
    labs(x = "Age", y = "Cumulative Incidence Rate (per 100,000)") +
    theme_classic() ->
    p

ggsave("reports/figures/cumulative_comparison.png", p, width = 4.5, height = 3.5, dpi = 300)

#############################################################################
# Next, stratify the cancer registry by percentiles

quantiles = quantile(full_data$phs, probs = c(0, 0.2, 0.5, 0.8, 1))

as.numeric(quantiles)
new_data = data.frame(phs = quantiles)
phs_prediction = predict(fit, new_data)
HR_80_50 = phs_prediction["80%"] - phs_prediction["50%"]
HR_20_50 = phs_prediction["20%"] - phs_prediction["50%"]

##############################################################################
# 3. Stratify incidence curve by HRs
#x_new
H0 = yearly_data$y_predict
#H0 = yearly_data$cumulative_incidence
H_80 = H0*exp(HR_80_50)
#S_80 = exp(-H_80)

H_20 = H0*exp(HR_20_50)
#S_20 = exp(-H_20)
adj_curve = data.frame(age = yearly_data$median_age_minus35,
                        H0,
                        H_20,
                        H_80) %>%
                pivot_longer(!age)
adj_curve$name <- factor(adj_curve$name, levels = c("H_80", "H0", "H_20"))
# ggplot(mapping = aes(x = x_new)) +
#     geom_line(aes(y = H0)) +
#     geom_line(aes(y = H_20)) +
#     geom_line(aes(y = H_80)) +
#     scale_x_continuous(breaks = c(5, 15, 25, 35, 45), labels = c(40, 50, 60, 70, 80)) +
#     labs(x = "Age", y = "Age-specific Incidence Rate (per 100,000)") +
#     theme_classic() ->
#     p
#ratio = with(adj_curve, ggthemes::bank_slopes(age, value, method = "as"))
ggplot(adj_curve, aes(age, value, color = name)) +
    geom_line() +
    scale_x_continuous(breaks = c(5, 15, 25, 35, 45, 55), 
                       labels = c(40, 50, 60, 70, 80, 90)) +
    labs(x = "Age", y = expression(atop("Age-specific Incidence Rate", "(per 100,000)"))) +
    scale_color_brewer(name = "", 
                         labels = c("PHS 80%", "Population Baseline", "PHS 20%"),
                         palette = "Set1") +
    lims(y = c(0, 2000)) +
    coord_fixed(ratio = 1/50) +
    theme_classic() ->
    p
ggsave("reports/figures/stratified_incidence_agespecific.png", p, width = 5, height = 3, dpi = 300)


#H0 = yearly_data$y_predict
H0 = yearly_data$cumulative_incidence
H_80 = H0*exp(HR_80_50)
#S_80 = exp(-H_80)

H_20 = H0*exp(HR_20_50)
#S_20 = exp(-H_20)
adj_curve = data.frame(age = yearly_data$median_age_minus35,
                        H0,
                        H_20,
                        H_80) %>%
                pivot_longer(!age)
adj_curve$name <- factor(adj_curve$name, levels = c("H_80", "H0", "H_20"))

ggplot(adj_curve, aes(age, value, color = name)) +
    geom_line() +
    scale_x_continuous(breaks = c(5, 15, 25, 35, 45, 55), 
                       labels = c(40, 50, 60, 70, 80, 90)) +
    labs(x = "Age", y = expression(atop("Cumulative Incidence Rate", "(per 100,000)"))) +
    scale_color_brewer(name = "", 
                         labels = c("PHS 80%", "Population Baseline", "PHS 20%"),
                         palette = "Set1") +
    lims(y = c(0, 40000)) +
    coord_fixed(ratio = 1/1000) +
    theme_classic() ->
    p
ggsave("reports/figures/stratified_incidence_cumulative.png", p, width = 5, height = 3, dpi = 300)

####################################################################################################

ggplot(mapping = aes(x = x_new)) +
    geom_line(aes(y = S_20)) +
    geom_line(aes(y = S_50)) +
    geom_line(aes(y = S_80)) +
    scale_x_continuous(breaks = c(5, 15, 25, 35, 45), labels = c(40, 50, 60, 70, 80)) +
    labs(x = "Age", y = "Survival Rate") +
    theme_classic() ->
    p
ggsave("reports/figures/stratified_survival.png", p, width = 4.5, height = 3.5, dpi = 300)




full_data %>% 

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