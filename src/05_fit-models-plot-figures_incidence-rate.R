library(tidyverse)
library(readxl)
library(survival)

source("src/00_config.R")


#########################################
# First argument: PHS file
# tab delimited file, no header, first column ID, second column PHS score
# Second argument: metadata file
# column names. Missing is okay, different order is okay, must match string exactly: 
# pid	death	wgs	age_lastobs	age_diagnosis	status	age	bfile_id	age_of_metastasis	
# metastasis_code	metastasis	conversion_to_metastasis_years	doed_av_pca	age_at_death	
# diagnosis_to_death_years	brca2_carrier	hoxb13_carrier	hoxb13_snp
# tab delimited file with column names bfile_id, age, status. 
# Additional columns will be ignored.
# Third argument: 
# boolean indicating whether to inverse (x * -1) the PHS scores
# to reverse the direction of effect.

args = commandArgs(trailingOnly=TRUE)

model_file = args[1] # PHS file
metadata_file = args[2] # Phenotype data
inverse = args[3]
tsd_proj = args[4]

model_file = "data/processed/PHS601_phs" # PHS file
metadata_file = "data/processed/heritable_cancer_phenotypes.tsv" # Phenotype data
# metadata_file = "data/processed/husk_phenotypes.tsv"
inverse = FALSE
tsd_proj = "p212"
# tsd_proj = "p2521"

model = gsub("(.+?)(\\_.*)", "\\1", basename(model_file))

#########################################
# Import data

# Huyhn-LE 2020 ASIS data
column_names = c("age", "low_risk", "fav_interm", "unfav_interm", "high_risk", "regional", "metastatic")
huynhle_asis = read_excel("data/external/ASIS_per100k_Huynh-Le_2020.xlsx", skip = 2, col_names = column_names)

huynhle_asis %>%
    rowwise() %>%
    mutate(total = sum(across(low_risk:metastatic))) %>%
    mutate(age_nice = if_else(age == "<50", "45-49", if_else(age == "≥85", "85-89", age))) %>%
    mutate(age_y = as.numeric(substr(age_nice, 1, 2))) %>%
    mutate(median_age = median(seq(as.numeric(substr(age_nice, 1, 2)), as.numeric(substr(age_nice, 4, 5))))  ) ->
    incidence

# import phs data
phs_data = import_phs_data(model_file, inverse)

# import phenotype data
metadata = read_tsv(metadata_file)

# combine phs and metadata
full_data = process_combined_data(metadata, phs_data)

####################################################################################
# Any PCA
####################################################################################

################################################################################
# 1. Plot incidence curve for Norway

# calculate natural log of incidence rate
incidence %>%
    mutate(log_incidencerate = log(total)) ->
    incidence


# Visualize age-specific incidence rate on log scale
# with polynomial fit
ggplot(incidence, aes(x = median_age, y = log_incidencerate)) +
    geom_point() +
    stat_smooth(method='lm', formula = y ~ poly(x,2)) +
    theme_classic() +
    # scale_x_continuous(breaks = c(5, 15, 25, 35, 45), labels = c(40, 50, 60, 70, 80)) +
    labs(x = "Age", y = "log(Age-specific Incidence Rate)")
#ggsave(paste0(results, "log_incidence_rate.png"), p, width = 4.5, height = 3.5, dpi = 300)

# Compare polynomial fit between fitting a linear model to the log-transformed
# data or fit a non-linear model to the original data
lm_fit = lm(log_incidencerate ~ poly(median_age, 2, raw = TRUE), data = incidence)
plot(incidence$median_age, resid(lm_fit))
abline(h = 0, lty = 2)

nls_fit <- nls(total ~ b0 * exp(b1 * median_age + b2 * median_age ^ 2), 
               data = incidence,
               start = list(b0 = -17, b1 = 0.6, b2 = -0.004))
coef(nls_fit)
plot(incidence$median_age, resid(nls_fit))
abline(h = 0, lty = 2)

B0 = lm_fit$coefficients[1]
B1 = lm_fit$coefficients[2]
B2 = lm_fit$coefficients[3]
pca_incidence_poly = function(x) {
    B0 + B1 * x + B2 * x ^ 2
}

x = incidence$median_age

new_data = incidence %>% select(median_age)

y_equation = pca_incidence_poly(x)
y_predict = predict(lm_fit, new_data) # should return the same

y = y_equation

ggplot(mapping = aes(x, y)) + 
    geom_smooth() + #loess
    geom_point(aes(x = median_age, y = log_incidencerate), data = incidence) +
    theme_classic() +
    # scale_x_continuous(breaks = c(5, 15, 25, 35, 45), labels = c(40, 50, 60, 70, 80)) +
    labs(x = "Age", y = "log(Age-specific Incidence Rate)")
#ggsave(paste0(results, "log_incidence_rate_lmfit.png"), p, width = 4.5, height = 3.5, dpi = 300)

###########################
# if I take the log of y and fit a second-order polynomial model I get this:
# log(y) = b0 + b1 * x + b2 * x ^ 2 
# To exponentiate it this back to linear space, we take the exp of both sizes
# the exp(log()) cancels out for y and we are left with:
# y = e ^ {b0 + b1 * x + b2 * x ^ 2}

b0 = B0
b1 = B1
b2 = B2

pca_incidence_poly_linearized = function(x) {
    exp(b0 + b1 * x + b2 * x ^ 2)
}

x_new = seq(from = min(incidence$median_age), to = max(incidence$median_age), by = 0.01)
y_lm = pca_incidence_poly_linearized(x_new)

lp <- predict(lm_fit, newdata = new_data)      # linear predictor on log scale
sigma2 <- summary(lm_fit)$sigma^2
y_mean <- exp(lp + 0.5 * sigma2)             # mean on original scale

# Exponentiate the model prediction to return it to the original units
ggplot() +
    geom_point(aes(x = median_age, y = total), incidence) +
    geom_smooth(aes(x_new, y_lm)) +
    # geom_point(aes(x = new_data$median_age, y = y_mean), color = "red") +
    # scale_x_continuous(breaks = c(5, 15, 25, 35, 45), labels = c(40, 50, 60, 70, 80)) +
    labs(x = "Age", y = "Age-specific Incidence Rate (per 100,000)") +
    theme_classic()
#ggsave(paste0(results, "log_incidence_rate_lmfit_exp.png"), p, width = 4.5, height = 3.5, dpi = 300)

# The model doesn't match up very well with the original in this case

# Try the nls model:
x_new = seq(from = min(incidence$median_age), to = max(incidence$median_age), by = 0.01)

new_data = data.frame(median_age = x_new)
y_predict = predict(nls_fit, new_data)
y_nls = y_predict
ggplot()  +
    geom_point(aes(x = median_age, y = total), incidence) +
    geom_smooth(aes(x = median_age, y = total), incidence, method = "loess", color = "green", se = F, span = 0.75) +
    geom_smooth(aes(x_new, y_nls)) +
    geom_smooth(aes(x_new, y_lm), color = "red") +
    # scale_x_continuous(breaks = c(5, 15, 25, 35, 45), labels = c(40, 50, 60, 70, 80)) +
    labs(x = "Age", y = "Age-specific Incidence Rate (per 100,000)") +
    theme_classic()
#ggsave(paste0(results, "norway_incidence.png"), p1, width = 4.5, height = 3.5, dpi = 300)

# The non-linear model seems to fit better

# Predict the incidence rate at the midpoint of each 1-year age group
yearly_data = data.frame(median_age = seq(from = min(incidence$median_age) + 0.5, to = max(incidence$median_age) - 0.5, by = 1))
yearly_data$y_predict = predict(nls_fit, yearly_data)
ggplot()  +
    geom_point(aes(x = median_age, y = total), incidence) +
    geom_smooth(aes(x_new, y_predict)) +
    geom_point(aes(median_age, y_predict), yearly_data, color = "red") +
    # scale_x_continuous(breaks = c(5, 15, 25, 35, 45), labels = c(40, 50, 60, 70, 80)) +
    labs(x = "Age", y = "Age-specific Incidence Rate (per 100,000)") +
    theme_classic()
# ggsave(paste0(results, "norway_incidence_yearly.png"), p1,
#            height = 8 * 0.618, width = 8)

# Now, take the cumulative incidence for later comparison to cohort data
yearly_data %>%
    mutate(cumulative_incidence = cumsum(y_predict)) ->
    yearly_data

ggplot(yearly_data, aes(x = median_age, y = cumulative_incidence)) +
    geom_point()

# ggsave(paste0(results, "norway_cumulative_incidence.png"), p,
#            height = 8 * 0.618, width = 8)

###################################################################################
# Plot the time-to-event data from the two cohorts next to the cancer registry data

fit = coxph(Surv(age, status) ~ phs, data = full_data)
H = data.frame(basehaz(fit)) # This is the cumulative baseline hazard
H %>%
    mutate(S = exp(-hazard)) %>% # This is the cumulative surival
    mutate(I = 1-S) %>%
    mutate(I_100000 = I * 100000) %>%
    mutate(time = time) -> # This is the cumulative incidence
    H


(ggplot() +
    geom_step(aes(time, I_100000), H, color = "red") +
    geom_line(aes(x = median_age, y = cumulative_incidence), yearly_data) +
    # scale_x_continuous(breaks = c(5, 15, 25, 35, 45, 55), labels = c(40, 50, 60, 70, 80, 90), limits = c(0, 55)) +
    labs(x = "Age", y = "Cumulative Incidence Rate (per 100,000)") +
    lims(y = c(0, 30000), x = c(40, 90)) +
    theme_classic() +
    theme(aspect.ratio = 1) ->
    p)

ggsave(paste0(results, paste(tsd_proj, model, "cumulative_comparison.png", sep = "_")), p,
           height = 6 * 0.618, width = 6)
ggsave(paste0(results, paste(tsd_proj, model, "cumulative_comparison.pdf", sep = "_")), p,
           height = 6 * 0.618, width = 6)
saveRDS(p, paste0(results, paste(tsd_proj, model, "cumulative_comparison.rds", sep = "_")))

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
adj_curve = data.frame(age = yearly_data$median_age,
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
ggsave(paste0(results, paste(tsd_proj, model, "stratified_incidence_agespecific.png", sep = "_")), 
       p, height = 8 * 0.618, width = 8)
ggsave(paste0(results, paste(tsd_proj, model, "stratified_incidence_agespecific.pdf", sep = "_")), 
       p, height = 8 * 0.618, width = 8)

#H0 = yearly_data$y_predict
H0 = yearly_data$cumulative_incidence
H_80 = H0*exp(HR_80_50)
#S_80 = exp(-H_80)

H_20 = H0*exp(HR_20_50)
#S_20 = exp(-H_20)
adj_curve = data.frame(age = yearly_data$median_age,
                        H0,
                        H_20,
                        H_80) %>%
                pivot_longer(!age)
adj_curve$name <- factor(adj_curve$name, levels = c("H_80", "H0", "H_20"))

ggplot(adj_curve, aes(age, value, color = name)) +
    geom_line() +
    # scale_x_continuous(breaks = c(5, 15, 25, 35, 45, 55), 
    #                    labels = c(40, 50, 60, 70, 80, 90)) +
    labs(x = "Age", y = expression(atop("Cumulative Incidence Rate", "(per 100,000)"))) +
    scale_color_brewer(name = "", 
                         labels = c("PHS 80%", "Population Baseline", "PHS 20%"),
                         palette = "Set1") +
    lims(y = c(0, 40000)) +
    coord_fixed(ratio = 1/1000) +
    theme_classic() ->
    p
# ggsave(paste0(results, paste(tsd_proj, model, "stratified_incidence_cumulative.png", sep = "_")), p,
#            height = 8 * 0.618, width = 8)
# ggsave(paste0(results, paste(tsd_proj, model, "stratified_incidence_cumulative.pdf", sep = "_")), p,
#            height = 8 * 0.618, width = 8)

####################################################################################################

# ggplot(mapping = aes(x = x_new)) +
#     geom_line(aes(y = S_20)) +
#     geom_line(aes(y = S_50)) +
#     geom_line(aes(y = S_80)) +
#     scale_x_continuous(breaks = c(5, 15, 25, 35, 45), labels = c(40, 50, 60, 70, 80)) +
#     labs(x = "Age", y = "Survival Rate") +
#     theme_classic() ->
#     p
# ggsave(paste0(results, paste(tsd_proj, model, "stratified_survival.png", sep = "_")), p, width = 4.5, height = 3.5, dpi = 300)




# full_data %>% 

# fit = coxph(Surv(age, status) ~ phs, data = full_data)
# summary(fit)
# new_data = data.frame(phs = quantile(full_data$phs, 0.80))
# percentile_80 = predict(fit, new_data)


# prev_norge = read_delim("data/external/Prevalence_20250320_163654.csv", delim = ";")
# prev_rogaland = read_delim("data/external/Prevalence_20250320_163603_Rogaland.csv", delim = ";")

# prev_rogaland$cases = as.numeric(gsub(" ", "", prev_rogaland$cases))
# prev_norge$cases = as.numeric(gsub(" ", "", prev_norge$cases))

# Try to plot cumulatice incidence with cohort data using coxph 

fit = coxph(Surv(age, status) ~ 1, data = full_data)
cohort_incidence = basehaz(fit)
cohort_incidence %>%
    mutate(cum_inc_per_100k = hazard * 100000) ->
    cohort_incidence
p1 = ggplot(cohort_incidence, aes(x = time, y = cum_inc_per_100k)) +
        geom_point()
# ggsave(paste0(results, paste(tsd_proj, model, "incidence.png", sep = "_")), p1)

# Plot cumulative incidence based on kreftregisteret data

incidence %>%
    mutate(ASR_cumsum = cumsum(total)) %>%
    mutate(Age_group = age_y) ->
    new_incidence
t = new_incidence$Age_group
y = new_incidence$ASR_cumsum

lm_fit = lm(log(ASR_cumsum) ~ poly(Age_group, 1, raw = TRUE), data = new_incidence)

add_nls <- nls(y ~ a*exp(r*t), 
               start = list(a = 1.3, r = 0.07))
nls_coef = coef(add_nls)
y = nls_coef[1]*exp(nls_coef[2]*t)
spline_int <- as.data.frame(spline(t, y))


########
# At least a plot of the cumulative incidence from cancer registry
# With age-specific rate per 100,000 on the y-axis
p2 = ggplot(spline_int, aes(x, y)) +
           geom_line() +
        geom_point(aes(x = Age_group, y = ASR_cumsum), data = new_incidence, color = "red")

# ggsave(paste0(results, "kreftregisteret_incidence.png"), p2)


# Incidence combined

p3 = ggplot() +
        geom_line(aes(x, y), data = spline_int) +
        geom_point(aes(time, cum_inc_per_100k), data = cohort_incidence)
# ggsave(paste0(results, paste(tsd_proj, model, "combined_incidence.png", sep = "_")), p3)


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
# ggsave(paste0(results, "incidence_test.png"), p1)


# ####################################################################################
# # Reshape p212 data into this format:
# # Age group (tiaar), Pop, Cases
# allData %>%
#     filter(SEX == "M") %>%
#     select(pid = PID, year_birth = YEAR_BIRTH, agelastobs = AGELASTOBS) ->
#     all_data_cleaned

# caseData %>%
#     filter(icd9 == "185") %>%
#     group_by(pid) %>%
#     filter(age_ca == min(age_ca)) %>%
#     ungroup() %>%
#     select(pid, yr_diagn, icd9, age_ca) ->
#     case_data_cleaned
    
# all_data_cleaned %>%
#     full_join(case_data_cleaned, by = "pid") ->
#     full_data

# ## Check for multiple diagnosis dates for the same individual
# # full_data %>% group_by(pid) %>% summarize(n = n()) %>% arrange(desc(n)) %>%
# # filter(n == 2) %>% pull(pid) -> repeat_pids

# # full_data %>%
# #     filter(pid %in% repeat_pids)

# full_data %>%
#     select(icd9, age_ca, agelastobs) %>%
#     mutate(case = if_else(is.na(icd9), 0, 1)) %>%
#     mutate(age = if_else(case == 1, age_ca, agelastobs)) %>%
#     select(age, case) %>%
#     mutate(decade = substr(age, 1, 1)) %>%
#     select(decade, case) %>%
#     group_by(decade) %>%
#     summarize(pop = n(),
#               cases = sum(case)) %>%
#     mutate(case_rate = cases / pop)

# # Incidience, i.e., year of diagnosis 

# Huyhn-LE 2020 ASIS data
column_names = c("age", "low_risk", "fav_interm", "unfav_interm", "high_risk", "regional", "metastatic")

huynhle_asis = read_excel("data/external/ASIS_per100k_Huynh-Le_2020.xlsx", skip = 2, col_names = column_names)

huynhle_asis %>%
    rowwise() %>%
    mutate(total = sum(across(low_risk:metastatic))) %>%
    mutate(age_nice = age)  ->
    huynhle_asis

huynhle_asis$age_nice[huynhle_asis$age == "<50"] = "45-49"
huynhle_asis$age_nice[huynhle_asis$age == "≥85"] = "85-89"

huynhle_asis %>%
    mutate(age_y = as.numeric(substr(age_nice, 1, 2))) ->
    huynhle_asis

huynhle_asis$cum_inc = cumsum(huynhle_asis$total)

ggplot(huynhle_asis, aes(x = age_y, y = total)) + 
    geom_smooth() + 
    theme_classic()

# 1. Fit curve to metastatic data
fm_lm <- lm(log(metastatic) ~ age_y, huynhle_asis)
st <- list(a = exp(coef(fm_lm)[1]), b = coef(fm_lm)[2])

#  a*exp(b*age_y)
nls_fit = nls(metastatic ~ a * exp(b * age_y), huynhle_asis, start = st)
coef(nls_fit)

x_new = seq(from = min(huynhle_asis$age_y), to = max(huynhle_asis$age_y), by = 0.01)

new_data = data.frame(age_y = x_new)
y_predict = predict(nls_fit, new_data)
y = y_predict
ggplot()  +
    geom_point(aes(x = age_y, y = metastatic), huynhle_asis) +
    geom_smooth(aes(x_new, y)) +
    scale_x_continuous(breaks = c(5, 15, 25, 35, 45), labels = c(40, 50, 60, 70, 80)) +
    labs(x = "Age", y = "Age-specific Incidence Rate (per 100,000)") +
    theme_classic() ->
    p1

###################################################################################
# Plot the time-to-event data from the two cohorts next to the cancer registry data

fit = coxph(Surv(age, status) ~ phs, data = full_data)
H = data.frame(basehaz(fit)) # This is the cumulative baseline hazard
H %>%
    mutate(S = exp(-hazard)) %>% # This is the cumulative surival
    mutate(I = 1-S) %>%
    mutate(I_100000 = I * 100000) %>%
    mutate(age_y = time) -> # This is the cumulative incidence
    H


ggplot() +
    geom_step(aes(age_y, I_100000), H, color = "red") +
    geom_line(aes(x = age_y, y = metastatic), huynhle_asis) +
    # scale_x_continuous(breaks = c(5, 15, 25, 35, 45, 55), labels = c(40, 50, 60, 70, 80, 90), limits = c(0, 55)) +
    labs(x = "Age", y = "Age-specific Incidence Rate (per 100,000)") +
    # lims(y = c(0, 30000)) +
    theme_classic() ->
    p

# ggsave(paste0(results, paste(tsd_proj, model, "cumulative_comparison.png", sep = "_")), p,
#            height = 6 * 0.618, width = 6)
# ggsave(paste0(results, paste(tsd_proj, model, "cumulative_comparison.pdf", sep = "_")), p,
#            height = 6 * 0.618, width = 6)

# 2. Multiply curve by beta from model. (metastatic model??)

#############################################################################
# Next, stratify the cancer registry by percentiles
yearly_data = data.frame(age_y = seq(from = min(huynhle_asis$age_y) + 0.5, to = max(huynhle_asis$age_y) + 0.5, by = 1))
yearly_data$y_predict = predict(nls_fit, yearly_data)

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
adj_curve = data.frame(age = yearly_data$age_y,
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

# 3. Plot curves and clean up design and save


ggplot(adj_curve, aes(age, value, color = name)) +
    geom_line() +
    labs(x = "Age", y = expression(atop("Age-specific Incidence Rate", "(per 100,000)"))) +
    scale_color_brewer(name = "", 
                         labels = c("PHS 80%", "Population Baseline", "PHS 20%"),
                         palette = "Set1") +
    theme_classic() ->
    p
ggsave(paste0(results, paste(tsd_proj, model, "stratified_incidence_agespecific.png", sep = "_")), 
       p, height = 6 * 0.618, width = 6)
ggsave(paste0(results, paste(tsd_proj, model, "stratified_incidence_agespecific.pdf", sep = "_")), 
       p, height = 6 * 0.618, width = 6)
saveRDS(p, paste0(results, paste(tsd_proj, model, "stratified_incidence_agespecific.rds", sep = "_")))
