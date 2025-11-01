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
        mutate(hazard = exp(predict(fit_reg, newdata = data.frame(median_age = age))),
            annual_risk = 1 - exp(-hazard),
            cum_H = cumsum(hazard),
            cum_risk = 1 - exp(-cum_H)) ->
            output_df

    return(output_df)
}

pca_risk = return_cum_risk(incidence, "rate_per_year", order = 2)
meta_risk = return_cum_risk(incidence, "meta_rate_per_year", order = 1)
cox_fit_pca <- coxph(Surv(age, status) ~ phs, data = full_data)
cox_fit_meta <- coxph(Surv(age_of_metastasis, metastasis) ~ phs, data = full_data)

# cohort data, pca_risk data
stratify_curve = function(full_data, pca_risk, fit_cox, output = "hazard") {
    # output must be "hazard" or "cum_risk"

    # Get predicted hazard at each quantile
    q <- quantile(full_data$phs, probs = c(0.2, 0.5, 0.8))
    eta <- predict(fit_cox, newdata = data.frame(phs = q), type = "lp")
    # HR between top20/bottom20 and median
    HR80 <- exp(eta[3] - eta[2])
    HR20 <- exp(eta[1] - eta[2])

    # stratify registry data (hazard) by high and low HR
    pca_risk %>%
        mutate(hazard_80 = hazard * HR80,
            annual_risk_80 = 1 - exp(-hazard_80),
            cum_risk_80 = 1 - exp(-cumsum(hazard_80))) %>%
        mutate(hazard_20 = pca_risk$hazard * HR20,
            annual_risk_20 = 1 - exp(-hazard_20),
            cum_risk_20 = 1 - exp(-cumsum(hazard_20))) ->
        pca_risk

    pca_risk %>%
        select(age, starts_with(output)) %>%
        pivot_longer(starts_with(output), names_to = "curve", values_to = output) %>%
        mutate(curve = factor(curve, levels = paste0(output, c("_80", "", "_20")))) ->
        stratified_df_long
    return(stratified_df_long)
}

stratified_hazard_pca = stratify_curve(full_data, pca_risk, cox_fit_pca)
stratified_hazard_meta = stratify_curve(full_data, meta_risk, cox_fit_meta)
stratified_risk_pca = stratify_curve(full_data, pca_risk, cox_fit_pca, output = "cum_risk")
stratified_risk_meta = stratify_curve(full_data, meta_risk, cox_fit_meta, output = "cum_risk")

group_labels <- c(
    cum_risk_80 = "PHS 80%",
    cum_risk = "Population baseline",
    cum_risk_20 = "PHS 20%"
)
group_labels <- c(
    "PHS 80%",
    "Population baseline",
    "PHS 20%"
)
plot_stratified = function(stratified_df, group_labels, y_var = "hazard") {
  # y_var must be "cum_risk" or "hazard"
  if (y_var == "hazard") {
    y_label = "Hazard (per 100,000)"
    stratified_df[, y_var] = stratified_df[, y_var] * 1e5
  } else if (y_var == "cum_risk") {
    y_label = "Cumulative Incidence (% Risk)"
    stratified_df[, y_var] = stratified_df[, y_var] * 1e2
  }
  y_var = sym(y_var)
  ggplot(stratified_df, aes(x = age, y = !!y_var, color = curve)) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A"), labels = group_labels) +
    labs(x = "Age", y = y_label, color = "") +
    theme_classic(base_size = 13) +
    theme(aspect.ratio = 1)
}
# Plot stratified (figure 4)
# p_stratified_hazard_pca = plot_stratified(stratified_hazard_pca, group_labels)
# p_stratified_hazard_meta = plot_stratified(stratified_hazard_meta, group_labels)
p_stratified_risk_pca = plot_stratified(stratified_risk_pca, group_labels, "cum_risk")
p_stratified_risk_meta = plot_stratified(stratified_risk_meta, group_labels, "cum_risk")

# saveRDS(p_stratified_hazard_pca, paste(model, "stratified_hazard_pca.rds", sep = "_"))
# saveRDS(p_stratified_hazard_meta, paste(model, "stratified_hazard_meta.rds", sep = "_"))
saveRDS(p_stratified_risk_pca, paste(model, "stratified_incidence.rds", sep = "_"))
saveRDS(p_stratified_risk_meta, paste(model, "stratified_incidence_metastatic.rds", sep = "_"))

# cumulative incidence function
return_cohort_ci = function(cox_fit) {
    bh <- basehaz(cox_fit, centered = FALSE)  # time column is age
    cohort_ci <- tibble(
        age   = bh$time,
        H0    = bh$hazard,
        F_cohort = 1 - exp(-H0)    # cumulative risk by age
    )
    return(cohort_ci)
}

cohort_ci_pca = return_cohort_ci(cox_fit_pca)
cohort_ci_meta = return_cohort_ci(cox_fit_meta)

# Two switches: pca/meta and hazard/risk
plot_comparison = function(pca_risk, cohort_ci, y_var = "hazard") {
  # y_var must be either "hazard" or "cum_risk"
  if (y_var == "hazard") {
    y_cohort = "H0"
    y_label = "Hazard (per 100,000)"
    pca_risk[, y_var] = pca_risk[, y_var] * 1e5
    cohort_ci[, y_cohort] = cohort_ci[, y_cohort] * 1e5
  } else if (y_var == "cum_risk") {
    y_cohort = "F_cohort"
    y_label = "Cumulative Incidence (% risk)"
    pca_risk[, y_var] = pca_risk[, y_var] * 1e2
    cohort_ci[, y_cohort] = cohort_ci[, y_cohort] * 1e2
  }
  y_var = sym(y_var)
  y_cohort = sym(y_cohort)
  ggplot() +
    geom_line(aes(age, !!y_var), pca_risk, color = "black", linewidth = 1) +
    geom_step(aes(age, !!y_cohort), cohort_ci, color = "red", linewidth = 1.1) +
    labs(x = "Age", y = y_label, color = "") +
    theme_classic(base_size = 13) +
    theme(aspect.ratio = 1)
}
# p_cumulative_comparison_pca = plot_comparison(pca_risk, cohort_ci_pca)
# p_cumulative_comparison_metastatic = plot_comparison(meta_risk, cohort_ci_meta)
p_cumulative_comparison_pca = plot_comparison(pca_risk, cohort_ci_pca, "cum_risk")
p_cumulative_comparison_metastatic = plot_comparison(meta_risk, cohort_ci_meta, "cum_risk")

saveRDS(p_cumulative_comparison_pca, paste(model, "cumulative_comparison.rds", sep = "_"))
saveRDS(p_cumulative_comparison_metastatic, paste(model, "cumulative_comparison_metastatic.rds", sep = "_"))
