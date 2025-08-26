library(tidyverse)
library(readxl)
library(survival)

source("src/00_config/config.R")
source("src/00_config/functions.R")

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
# metadata_file = "data/processed/heritable_cancer_phenotypes.tsv" # Phenotype data
metadata_file = "data/processed/husk_phenotypes.tsv"
inverse = FALSE
# tsd_proj = "p212"
tsd_proj = "p2521"

model = gsub("(.+?)(\\_.*)", "\\1", basename(model_file))

#########################################
# Import data

# Huyhn-LE 2020 ASIS data
column_names = c("age", "low_risk", "fav_interm", "unfav_interm", "high_risk", "regional", "metastatic")
huynhle_asis = read_excel("data/external/ASIS_per100k_Huynh-Le_2020.xlsx", skip = 2, col_names = column_names)

huynhle_asis %>%
    rowwise() %>%
    mutate(total = sum(across(low_risk:metastatic))) %>%
    mutate(reg_met = sum(across(regional:metastatic))) %>%
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

full_data = full_data %>% filter(fdato >= "1950-01-01")
# # Troubleshooting

# response =  paste0("log(", incidence_variable, ")")
#     predictor = "poly(median_age, 2, raw = TRUE)"
#     predictor = "median_age"
# lm_fit = lm(log(total) ~ poly(median_age, 2, raw = TRUE), data = incidence)
# lm_fit = lm(log(total) ~ poly(median_age, 2, raw = TRUE), data = incidence)
# ggplot()

    ## Fit the nls model from the starting parameters of the log model
    predictor = "b0 * exp(b1 * median_age + b2 * median_age ^ 2)"
    predictor = "b0 * exp(b1 * median_age)"
    nls_fit <- nls(as.formula(sprintf('%s~%s', incidence_variable, predictor)),
               data = incidence,
            #    start = list(b0 = coef(lm_fit)[1], b1 = coef(lm_fit)[2], b2 = coef(lm_fit)[3]))
               start = list(b0 = coef(lm_fit)[1], b1 = coef(lm_fit)[2]))
    # Create a smooth curve by predicting very fine time steps (here: 1% of a year)
    x_new = seq(from = min(incidence$median_age), to = max(incidence$median_age), by = 0.01)
    new_data = data.frame(median_age = x_new)   # new data frame for prediction
    y_nls = predict(nls_fit, new_data)          # prediction from the nls model
    # to calculate the annualized incidence rate, 
    # predict at a 1 year time step centered on the mid-point of the year
    yearly_data = data.frame(median_age = seq(from = min(incidence$median_age) + 0.5, to = max(incidence$median_age) - 0.5, by = 1))
    yearly_data %>%
        mutate(y_predict = predict(nls_fit, yearly_data)) %>%
        mutate(cumulative_incidence = cumsum(y_predict)) ->
        yearly_data
    return(yearly_data)

####################################################################################
# Function definitions
####################################################################################

fit_incidence_model = function(incidence, incidence_variable = "total") {
    response =  paste0("log(", incidence_variable, ")")
    predictor = "poly(median_age, 2, raw = TRUE)"
    predictor = "median_age"
    lm_fit = lm(as.formula(sprintf('%s~%s', response, predictor)), data = incidence)
    ## Fit the nls model from the starting parameters of the log model
    predictor = "b0 * exp(b1 * median_age + b2 * median_age ^ 2)"
    predictor = "b0 * exp(b1 * median_age)"
    nls_fit <- nls(as.formula(sprintf('%s~%s', incidence_variable, predictor)),
               data = incidence,
            #    start = list(b0 = coef(lm_fit)[1], b1 = coef(lm_fit)[2], b2 = coef(lm_fit)[3]))
               start = list(b0 = coef(lm_fit)[1], b1 = coef(lm_fit)[2]))
    # Create a smooth curve by predicting very fine time steps (here: 1% of a year)
    x_new = seq(from = min(incidence$median_age), to = max(incidence$median_age), by = 0.01)
    new_data = data.frame(median_age = x_new)   # new data frame for prediction
    y_nls = predict(nls_fit, new_data)          # prediction from the nls model
    # to calculate the annualized incidence rate, 
    # predict at a 1 year time step centered on the mid-point of the year
    yearly_data = data.frame(median_age = seq(from = min(incidence$median_age) + 0.5, to = max(incidence$median_age) - 0.5, by = 1))
    yearly_data %>%
        mutate(y_predict = predict(nls_fit, yearly_data)) %>%
        mutate(cumulative_incidence = cumsum(y_predict)) ->
        yearly_data
    return(yearly_data)
}

plot_comparison_figure = function(yearly_data, cohort_data, age = "age", status = "status") {
    response = paste0("Surv(", age, ", ", status, ")")
    pca_fit = coxph(as.formula(sprintf('%s~%s', response, "phs")), data = cohort_data)
    H = data.frame(basehaz(pca_fit))    # This is the cumulative baseline hazard
    H %>%
        mutate(S = exp(-hazard)) %>%        # This is the cumulative survival
        mutate(I = 1-S) %>%                 # The inverse of survival is incidence
        mutate(I_1e5 = I * 1e5) ->    # Convert to units of per 100k
        H
    ggplot() +
        geom_step(aes(time, I_1e5), H, color = "red") +
        geom_line(aes(x = median_age, y = cumulative_incidence), yearly_data) +
        # scale_x_continuous(breaks = c(5, 15, 25, 35, 45, 55), labels = c(40, 50, 60, 70, 80, 90), limits = c(0, 55)) +
        labs(x = "Age", y = "Cumulative Incidence Rate (per 100,000)") +
        # lims(y = c(0, 30000), x = c(40, 90)) +
        lims(x = c(40, 90)) +
        theme_classic() +
        theme(aspect.ratio = 1)
}

plot_stratified_figure = function(yearly_data, cohort_data, age = "age", status = "status") {
    response = paste0("Surv(", age, ", ", status, ")")
    pca_fit = coxph(as.formula(sprintf('%s~%s', response, "phs")), data = cohort_data)
    quantiles = quantile(cohort_data$phs, probs = c(0, 0.2, 0.5, 0.8, 1))
    new_data = data.frame(phs = quantiles)
    phs_prediction = predict(pca_fit, new_data)
    HR_80_50 = phs_prediction["80%"] - phs_prediction["50%"]
    HR_20_50 = phs_prediction["20%"] - phs_prediction["50%"]
    H0 = yearly_data$y_predict
    H_80 = H0*exp(HR_80_50)
    H_20 = H0*exp(HR_20_50)
    adj_curve = data.frame(age = yearly_data$median_age,
                        H0,
                        H_20,
                        H_80) %>%
                pivot_longer(!age)
    adj_curve$name <- factor(adj_curve$name, levels = c("H_80", "H0", "H_20"))
    ggplot(adj_curve, aes(age, value, color = name)) +
    geom_line() +
    labs(x = "Age", y = expression(atop("Age-specific Incidence Rate", "(per 100,000)"))) +
    scale_color_brewer(name = "", 
                         labels = c("PHS 80%", "Population Baseline", "PHS 20%"),
                         palette = "Set1") +
    theme_classic()
}

######################################################################################
# Any PCa
yearly_data = fit_incidence_model(incidence, "total")

p = plot_comparison_figure(yearly_data, full_data)

ggsave(paste0(results, paste(tsd_proj, model, "cumulative_comparison.png", sep = "_")), p,
           height = 6 * 0.618, width = 6)
ggsave(paste0(results, paste(tsd_proj, model, "cumulative_comparison.pdf", sep = "_")), p,
           height = 6 * 0.618, width = 6)
saveRDS(p, paste0(results, paste(tsd_proj, model, "cumulative_comparison.rds", sep = "_")))

p = plot_stratified_figure(yearly_data, full_data)

ggsave(paste0(results, paste(tsd_proj, model, "stratified_incidence_agespecific.png", sep = "_")), 
       p, height = 6 * 0.618, width = 6)
ggsave(paste0(results, paste(tsd_proj, model, "stratified_incidence_agespecific.pdf", sep = "_")), 
       p, height = 6 * 0.618, width = 6)
saveRDS(p, paste0(results, paste(tsd_proj, model, "stratified_incidence_agespecific.rds", sep = "_")))

######################################################################################
# Metastatic PCa
yearly_data = fit_incidence_model(incidence, "reg_met")

p = plot_comparison_figure(yearly_data, full_data, age = "age_of_metastasis", status = "metastasis")

ggsave("reg_met_comparison_loose.png", p, height = 6 * 0.618, width = 6)

ggsave(paste0(results, paste(tsd_proj, model, "cumulative_comparison_metastatic.png", sep = "_")), p,
           height = 6 * 0.618, width = 6)
ggsave(paste0(results, paste(tsd_proj, model, "cumulative_comparison_metastatic.pdf", sep = "_")), p,
           height = 6 * 0.618, width = 6)
saveRDS(p, paste0(results, paste(tsd_proj, model, "cumulative_comparison_metastatic.rds", sep = "_")))

p = plot_stratified_figure(yearly_data, full_data, age = "age_of_metastasis", status = "metastasis")

ggsave(paste0(results, paste(tsd_proj, model, "stratified_incidence_agespecific_metastatic.png", sep = "_")), 
       p, height = 6 * 0.618, width = 6)
ggsave(paste0(results, paste(tsd_proj, model, "stratified_incidence_agespecific_metastatic.pdf", sep = "_")), 
       p, height = 6 * 0.618, width = 6)
saveRDS(p, paste0(results, paste(tsd_proj, model, "stratified_incidence_agespecific_metastatic.rds", sep = "_")))
