library(tidyverse)
library(readxl)
library(survival)
library(hazrd)
source("src/config.R")

args <- commandArgs(trailingOnly = TRUE)

metadata_file <- args[1] # Phenotype data
output_prefix <- args[2]

column_names <- c(
  "age", "low_risk", "fav_interm", "unfav_interm",
  "high_risk", "regional", "metastatic"
)
huynhle_asis <- read_excel("data/external/ASIRs_per100k_Huynh-Le_2020.xlsx",
  skip = 2, col_names = column_names
)

incidence <-
  huynhle_asis %>%
  rowwise() %>%
  mutate(total = sum(across(low_risk:metastatic))) %>%
  mutate(reg_met = sum(across(regional:metastatic))) %>%
  mutate(age_nice = if_else(age == "<50", "45-49",
    if_else(age == "≥85", "85-89", age)
  )) %>%
  mutate(age_y = as.numeric(substr(age_nice, 1, 2))) %>%
  mutate(median_age = median(seq(
    as.numeric(substr(age_nice, 1, 2)),
    as.numeric(substr(age_nice, 4, 5))
  ))) %>%
  mutate(rate_per_year = total / 1e5) %>%
  mutate(meta_rate_per_year = reg_met / 1e5)

# import data
full_data <- read_tsv(metadata_file)

# Compute cumulative risk for the first- or second-order polynomial curves
# Takes registry data, response column, and linear or poly
return_cum_risk <- function(incidence_data, response, order = 1) {
  if (order == 1) {
    fit_reg <- lm(as.formula(sprintf("log(%s)~median_age", response)),
      data = incidence_data
    )
  } else if (order == 2) {
    fit_reg <- lm(
      as.formula(sprintf(
        "log(%s)~poly(median_age, 2, raw = TRUE)",
        response
      )),
      data = incidence_data
    )
  } else {
    stop("Invalid value for order. Must be 1 or 2")
  }
  # mid-year ages
  output_df <-
    tibble(age = seq(min(incidence$median_age) + 0.5, max(incidence$median_age) - 0.5, by = 1)) %>%
    dplyr::mutate(
      hazard = exp(predict(fit_reg, newdata = data.frame(median_age = age))),
      annual_risk = 1 - exp(-hazard),
      cum_H = cumsum(hazard),
      cum_risk = 1 - exp(-cum_H)
    )

  return(output_df)
}

pca_risk <- return_cum_risk(incidence, "rate_per_year", order = 2)
meta_risk <- return_cum_risk(incidence, "meta_rate_per_year", order = 1)
cox_fit_pca <- coxph(Surv(age, status) ~ phs, data = full_data)
cox_fit_meta <- coxph(Surv(age_of_metastasis, metastasis) ~ phs,
  data = full_data
)

# cohort data, pca_risk data
# Registry hazard is treated as fixed/known; only the cohort's bootstrapped HR
# estimate contributes uncertainty, so the CI ribbon reflects HR uncertainty only.
stratify_curve <- function(full_data, pca_risk, fit_cox, output = "hazard") {
  # output must be "hazard" or "cum_risk"
  hr_bands <- list(
    "95" = c(0.95, 1.00),
    "80" = c(0.80, 1.00),
    "20" = c(0.00, 0.20)
  )

  hr_summary <- purrr::map_dfr(names(hr_bands), function(band) {
    res <- phs_metrics(
      full_data,
      metrics   = "HR",
      hr_pairs  = list(list(numerator = hr_bands[[band]], denominator = c(0.30, 0.70))),
      bootstrap = TRUE,
      n_boot    = bootstrap_iters,
      seed      = 42
    )
    tibble(
      band = band, estimate = res$estimate,
      conf_low = res$conf_low, conf_high = res$conf_high
    )
  })

  # stratify registry data (hazard) by cohort HR point estimate and its bootstrap CI
  for (band in names(hr_bands)) {
    hr <- hr_summary[hr_summary$band == band, ]
    pca_risk[[paste0("hazard_", band)]] <- pca_risk$hazard * hr$estimate
    pca_risk[[paste0("cum_risk_", band)]] <-
      1 - exp(-cumsum(pca_risk[[paste0("hazard_", band)]]))
    pca_risk[[paste0("cum_risk_", band, "_lo")]] <-
      1 - exp(-cumsum(pca_risk$hazard * hr$conf_low))
    pca_risk[[paste0("cum_risk_", band, "_hi")]] <-
      1 - exp(-cumsum(pca_risk$hazard * hr$conf_high))
  }
  # Population baseline (HR = 1) has no bootstrapped CI; no ribbon is drawn for it.
  pca_risk$cum_risk_lo <- NA_real_
  pca_risk$cum_risk_hi <- NA_real_

  estimate_cols <- c(paste0(output, "_95"), paste0(output, "_80"), output, paste0(output, "_20"))
  if (output == "cum_risk") {
    lo_cols <- paste0(estimate_cols, "_lo")
    hi_cols <- paste0(estimate_cols, "_hi")
    # baseline column has no "_lo"/"_hi" suffix variants; use the NA placeholders added above
    lo_cols[estimate_cols == output] <- "cum_risk_lo"
    hi_cols[estimate_cols == output] <- "cum_risk_hi"
  } else {
    # No bootstrapped CI computed for the raw (non-cumulative) hazard curves.
    pca_risk[paste0(estimate_cols, "_lo")] <- NA_real_
    pca_risk[paste0(estimate_cols, "_hi")] <- NA_real_
    lo_cols <- paste0(estimate_cols, "_lo")
    hi_cols <- paste0(estimate_cols, "_hi")
  }

  stratified_df_long <-
    pca_risk %>%
    select(age, all_of(estimate_cols)) %>%
    pivot_longer(all_of(estimate_cols), names_to = "curve", values_to = output) %>%
    mutate(curve = factor(curve, levels = estimate_cols))

  ci_long <-
    map_dfr(seq_along(estimate_cols), function(i) {
      tibble(
        age       = pca_risk$age,
        curve     = factor(estimate_cols[i], levels = estimate_cols),
        conf_low  = pca_risk[[lo_cols[i]]],
        conf_high = pca_risk[[hi_cols[i]]]
      )
    })

  left_join(stratified_df_long, ci_long, by = c("age", "curve"))
}

stratified_hazard_pca <- stratify_curve(full_data, pca_risk, cox_fit_pca)
stratified_hazard_meta <- stratify_curve(full_data, meta_risk, cox_fit_meta)
stratified_risk_pca <- stratify_curve(full_data, pca_risk, cox_fit_pca, output = "cum_risk")
stratified_risk_meta <- stratify_curve(full_data, meta_risk, cox_fit_meta, output = "cum_risk")

group_labels <- c(
  "PHS 95-100th",
  "PHS 80-100th",
  "Population baseline",
  "PHS 0-20th"
)
plot_stratified <- function(stratified_df, group_labels, y_var = "hazard") {
  # y_var must be "cum_risk" or "hazard"
  if (y_var == "hazard") {
    y_label <- "Hazard (per 100,000)"
    stratified_df[, y_var] <- stratified_df[, y_var] * 1e5
    stratified_df[, c("conf_low", "conf_high")] <- stratified_df[, c("conf_low", "conf_high")] * 1e5
  } else if (y_var == "cum_risk") {
    y_label <- "Cumulative incidence, %\n(Absolute risk)"
    stratified_df[, y_var] <- stratified_df[, y_var] * 1e2
    stratified_df[, c("conf_low", "conf_high")] <- stratified_df[, c("conf_low", "conf_high")] * 1e2
  }
  y_var <- sym(y_var)
  ggplot(stratified_df, aes(x = age, y = !!y_var, color = curve, fill = curve)) +
    geom_ribbon(aes(ymin = conf_low, ymax = conf_high), alpha = 0.15, color = NA, na.rm = TRUE) +
    geom_line(linewidth = 0.8) +
    ggplot2::scale_color_brewer(palette = "Set1", labels = group_labels) +
    ggplot2::scale_fill_brewer(palette = "Set1", labels = group_labels) +
    labs(x = "Age", y = y_label, color = "", fill = "") +
    theme_classic(base_size = 12) +
    theme(aspect.ratio = 1)
}
# Plot stratified (figure 4)
p_stratified_risk_pca <- plot_stratified(stratified_risk_pca, group_labels, "cum_risk")
p_stratified_risk_meta <- plot_stratified(stratified_risk_meta, group_labels, "cum_risk")

saveRDS(p_stratified_risk_pca, paste(output_prefix, "stratified_incidence.rds", sep = "_"))
saveRDS(p_stratified_risk_meta, paste(output_prefix, "stratified_incidence_metastatic.rds", sep = "_"))

# cumulative incidence function
return_cohort_ci <- function(cox_fit) {
  bh <- basehaz(cox_fit, centered = FALSE) # time column is age
  cohort_ci <- tibble(
    age = bh$time,
    H0 = bh$hazard,
    F_cohort = 1 - exp(-H0) # cumulative risk by age
  )
  return(cohort_ci)
}

cohort_ci_pca <- return_cohort_ci(cox_fit_pca)
cohort_ci_meta <- return_cohort_ci(cox_fit_meta)

# Two switches: pca/meta and hazard/risk
plot_comparison <- function(pca_risk, cohort_ci, y_var = "hazard") {
  # y_var must be either "hazard" or "cum_risk"
  if (y_var == "hazard") {
    y_cohort <- "H0"
    y_label <- "Hazard (per 100,000)"
    pca_risk[, y_var] <- pca_risk[, y_var] * 1e5
    cohort_ci[, y_cohort] <- cohort_ci[, y_cohort] * 1e5
  } else if (y_var == "cum_risk") {
    y_cohort <- "F_cohort"
    y_label <- "Cumulative Incidence (% risk)"
    pca_risk[, y_var] <- pca_risk[, y_var] * 1e2
    cohort_ci[, y_cohort] <- cohort_ci[, y_cohort] * 1e2
  }
  y_var <- sym(y_var)
  y_cohort <- sym(y_cohort)
  ggplot() +
    geom_line(aes(age, !!y_var), pca_risk, color = "black", linewidth = 1) +
    geom_step(aes(age, !!y_cohort), cohort_ci, color = "red", linewidth = 1) +
    labs(x = "Age", y = y_label, color = "") +
    theme_classic(base_size = 12) #+
    # theme(aspect.ratio = 1)
}
p_cumulative_comparison_pca <- plot_comparison(pca_risk, cohort_ci_pca, "cum_risk")
p_cumulative_comparison_metastatic <- plot_comparison(meta_risk, cohort_ci_meta, "cum_risk")

saveRDS(p_cumulative_comparison_pca, paste(output_prefix, "cumulative_comparison.rds", sep = "_"))
saveRDS(p_cumulative_comparison_metastatic, paste(output_prefix, "cumulative_comparison_metastatic.rds", sep = "_"))
