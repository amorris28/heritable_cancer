library(tidyverse)
library(readxl)
library(survival)

source("src/config.R")
source("src/functions.R")

args = commandArgs(trailingOnly=TRUE)

model_file = args[1] # PHS file
metadata_file = args[2] # Phenotype data
inverse = args[3]

model = tools::file_path_sans_ext(basename(model_file))
column_names = c("age", "low_risk", "fav_interm", "unfav_interm", "high_risk", "regional", "metastatic")
huynhle_asis = read_excel("data/external/ASIRs_per100k_Huynh-Le_2020.xlsx", skip = 2, col_names = column_names)

huynhle_asis %>%
    rowwise() %>%
    mutate(total = sum(across(low_risk:metastatic))) %>%
    mutate(reg_met = sum(across(regional:metastatic))) %>%
    mutate(age_nice = if_else(age == "<50", "45-49", if_else(age == "≥85", "85-89", age))) %>%
    mutate(age_y = as.numeric(substr(age_nice, 1, 2))) %>%
    mutate(median_age = median(seq(as.numeric(substr(age_nice, 1, 2)), as.numeric(substr(age_nice, 4, 5))))) %>%
    mutate(rate_per_year = total / 1e5) %>%
    mutate(meta_rate_per_year = reg_met / 1e5) ->
    incidence

# import phs data
phs_data = import_phs_data(model_file, inverse)

# import phenotype data
metadata = read_tsv(metadata_file)

# combine phs and metadata
full_data = process_combined_data(metadata, phs_data)

# Compute the various curves in a function that takes arguments for
# registry data, age column, response column, linear model or polynomial
return_cum_risk = function(incidence_data, response, order = 1) {
    if (order == 1) {
        fit_reg <- lm(as.formula(sprintf('log(%s)~median_age', response)), data = incidence_data)
    } else if (order == 2){
        fit_reg <- lm(as.formula(sprintf('log(%s)~poly(median_age, 2, raw = TRUE)', response)), data = incidence_data)
    } else {
        stop("Invalid value for order. Must be 1 or 2")
    }
    # mid-year ages
    tibble(age = seq(min(incidence$median_age) + 0.5, max(incidence$median_age) - 0.5, by = 1)) %>%
        mutate(haz_reg = exp(predict(fit_reg, newdata = data.frame(median_age = age))),
            annual_risk = 1 - exp(-haz_reg),
            cum_H = cumsum(haz_reg),
            cum_risk = 1 - exp(-cum_H)) ->
            output_df

    return(output_df)
}

pca_risk = return_cum_risk(incidence, "rate_per_year", order = 2)
meta_risk = return_cum_risk(incidence, "meta_rate_per_year", order = 1)
cox_fit_pca <- coxph(Surv(age, status) ~ phs, data = full_data)
cox_fit_meta <- coxph(Surv(age_of_metastasis, metastasis) ~ phs, data = full_data)

# cohort data, pca_risk data
stratify_risk = function(full_data, pca_risk, fit_cox) {
    # Fit cox model of cohort data

    # Get predicted hazard at each quantile
    q <- quantile(full_data$phs, probs = c(0.2, 0.5, 0.8))
    eta <- predict(fit_cox, newdata = data.frame(phs = q), type = "lp")
    # HR between top20/bottom20 and median
    HR80 <- exp(eta[3] - eta[2])
    HR20 <- exp(eta[1] - eta[2])

    # stratify registry data (haz_reg) by high and low HR
    pca_risk %>%
        mutate(haz_reg_80 = haz_reg * HR80,
            annual_risk_80 = 1 - exp(-haz_reg_80),
            cum_risk_80 = 1 - exp(-cumsum(haz_reg_80))) %>%
        mutate(haz_reg_20 = pca_risk$haz_reg * HR20,
            annual_risk_20 = 1 - exp(-haz_reg_20),
            cum_risk_20 = 1 - exp(-cumsum(haz_reg_20))) ->
        pca_risk
    pca_risk %>%
        select(age, starts_with("cum_risk")) %>%
        pivot_longer(starts_with("cum_risk"), names_to = "curve", values_to = "cum_risk") %>%
        mutate(curve = factor(curve, levels = c("cum_risk_80", "cum_risk", "cum_risk_20"))) ->
        stratified_df_long
    return(stratified_df_long)
}

stratified_df_pca = stratify_risk(full_data, pca_risk, cox_fit_pca)
stratified_df_meta = stratify_risk(full_data, meta_risk, cox_fit_meta)

group_labels <- c(
    cum_risk_80 = "PHS 80%",
    cum_risk = "Population baseline",
    cum_risk_20 = "PHS 20%"
)

# Plot stratified (figure 4)
p = ggplot(stratified_df_pca, aes(x = age, y = cum_risk, color = curve)) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A"), labels = group_labels) +
    labs(x = "Age", y = "Cumulative Incidence", color = "", 
         title = "Absolute risk (cumulative incidence) by risk group") +
    theme_classic(base_size = 13) +
    theme(aspect.ratio = 1)
saveRDS(p, paste(model, "stratified_incidence_agespecific.rds", sep = "_"))

# Plot stratified (figure 4)
p = ggplot(stratified_df_meta, aes(x = age, y = cum_risk, color = curve)) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A"), labels = group_labels) +
    labs(x = "Age", y = "Cumulative Incidence", color = "", 
         title = "Absolute risk (cumulative incidence) by risk group") +
    theme_classic(base_size = 13) +
    theme(aspect.ratio = 1)
saveRDS(p, paste(model, "stratified_incidence_agespecific_metastatic.rds", sep = "_"))

# cumulative incidence function
return_cohort_ci = function(cox_fit) {
    bh <- basehaz(cox_fit, centered = FALSE)  # time column is age
    cohort_ci <- tibble(
        age   = bh$time,
        H0    = bh$hazard,
        F_cohort = 1 - exp(-H0)    # cumulative risk by age
    )
}
cohort_ci_pca = return_cohort_ci(cox_fit_pca)
p = ggplot() +
    geom_step(aes(age, F_cohort), cohort_ci_pca, color = "red", linewidth = 1.1) +
    geom_line(aes(age, cum_risk), pca_risk, color = "black", linewidth = 1) +
    labs(x = "Age", y = "Cumulative Incidence (per 100,000)",
            color = "", title = "Absolute risk (cumulative incidence) by age") +
    theme_classic(base_size = 13) +
    theme(aspect.ratio = 1)
saveRDS(p, paste(model, "cumulative_comparison.rds", sep = "_"))

cohort_ci_meta = return_cohort_ci(cox_fit_meta)
p = ggplot() +
    geom_step(aes(age, F_cohort), cohort_ci_meta, color = "red", linewidth = 1.1) +
    geom_line(aes(age, cum_risk), meta_risk, color = "black", linewidth = 1) +
    labs(x = "Age", y = "Cumulative Incidence (per 100,000)",
            color = "", title = "Absolute risk (cumulative incidence) by age") +
    theme_classic(base_size = 13) +
    theme(aspect.ratio = 1)
saveRDS(p, paste(model, "cumulative_comparison_metastatic.rds", sep = "_"))

