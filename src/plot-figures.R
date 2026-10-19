library(tidyverse)
library(stringr)
library(survival)
library(patchwork)
source("src/config.R")
source("src/functions.R")

args <- commandArgs(trailingOnly = TRUE)

metadata_file <- args[1] # Phenotype data
output_prefix <- args[2]

#########################################
# Import PHS data and metadata
combined_data <- read_tsv(metadata_file)

#########################################
combined_data <-
  combined_data %>%
  filter(age >= 40) %>%
  mutate(age_minus_40 = age - 40) %>%
  mutate(age_of_metastasis_minus_40 = age_of_metastasis - 40)

# ZPH TT
# model
fit <- coxph(Surv(age, status) ~ phs, data = combined_data)

# Residuals from cox zph fit
zph_fit <- cox.zph(fit)

# Create data frames for plotting
resid_data <- data.frame(
  time = zph_fit$time,
  beta = zph_fit$y[, "phs"] # Schoenfeld residuals
)

# Get the coefficient from the original model
coef_phs <- coef(fit)["phs"]

# Fit the time-transformed model
phs_fit_tt <- coxph(Surv(age_minus_40, status) ~ phs + tt(phs), combined_data,
  tt = function(x, t, ...) x * t
)

# Extract coefficients and variance-covariance matrix
coefs <- coef(phs_fit_tt)
vcov_matrix <- vcov(phs_fit_tt)

# Coefficients
beta_phs <- coefs[1] # Main effect
beta_tt <- coefs[2] # Time interaction effect

# Standard errors
se_phs <- sqrt(vcov_matrix[1, 1])
se_tt <- sqrt(vcov_matrix[2, 2])
cov_phs_tt <- vcov_matrix[1, 2] # Covariance between the two

# Create sequence of time points
time_seq <- seq(min(combined_data$age_minus_40, na.rm = TRUE),
  max(combined_data$age_minus_40, na.rm = TRUE),
  length.out = 100
)

# Calculate time-varying coefficient: β(t) = β_phs + β_tt * t
beta_t <- beta_phs + beta_tt * time_seq

# Calculate standard error of β(t)
# SE[β(t)] = sqrt(Var(β_phs) + t^2 * Var(β_tt) + 2*t*Cov(β_phs, β_tt))
se_beta_t <- sqrt(
  vcov_matrix[1, 1] + time_seq^2 * vcov_matrix[2, 2] + 2 * time_seq * cov_phs_tt
)

# Calculate 95% confidence intervals
z_value <- qnorm(0.975) # Cox models use normal distribution
ci_lower <- beta_t - z_value * se_beta_t
ci_upper <- beta_t + z_value * se_beta_t

# Create plotting data frame
plot_data <- data.frame(
  time = time_seq + 40,
  beta_t = beta_t,
  ci_lower = ci_lower,
  ci_upper = ci_upper
)

# Create the plot
zph_tt <- ggplot(plot_data, aes(x = time, y = beta_t)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),
    fill = "steelblue", alpha = 0.3
  ) +
  geom_line(color = "steelblue", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_point(data = resid_data, aes(x = time, y = beta)) +
  stat_smooth() +
  geom_hline(
    yintercept = coef_phs,
    linetype = "dashed",
    color = "blue"
  ) +
  labs(
    title = "Schoenfeld Residuals for phs",
    subtitle = paste0(
      "Test for PH: p = ",
      round(zph_fit$table["phs", "p"], 4)
    ),
    x = "Time (age)",
    y = "Beta(t) for phs"
  ) +
  theme_classic()
saveRDS(zph_tt, paste(output_prefix, "zph_tt.rds", sep = "_"))

#########################################################
# Plot Cumulative Incidence Curves
any_pca_cuminc <- phs_cuminc_strata(
  combined_data,
  intervals = list(c(0.95, 1), c(0.80, 1), c(0.30, 0.70), c(0, 0.20))
)
any_pca_km_curve <- plot_cuminc_strata(any_pca_cuminc)

######################################################
# metastatic
#####################################################
# ZPH TT
# model
fit <- coxph(Surv(age_of_metastasis, metastasis) ~ phs, data = combined_data)

# Residuals from cox zph fit
zph_fit <- cox.zph(fit)

# Create data frames for plotting
resid_data <- data.frame(
  time = zph_fit$time,
  beta = zph_fit$y[, "phs"] # Schoenfeld residuals
)

# Get the coefficient from the original model
coef_phs <- coef(fit)["phs"]

# Fit the time-transformed model
phs_fit_tt <- coxph(
  Surv(age_of_metastasis_minus_40, metastasis) ~ phs + tt(phs), combined_data,
  tt = function(x, t, ...) x * t
)

# Extract coefficients and variance-covariance matrix
coefs <- coef(phs_fit_tt)
vcov_matrix <- vcov(phs_fit_tt)

# Coefficients
beta_phs <- coefs[1] # Main effect
beta_tt <- coefs[2] # Time interaction effect

# Standard errors
se_phs <- sqrt(vcov_matrix[1, 1])
se_tt <- sqrt(vcov_matrix[2, 2])
cov_phs_tt <- vcov_matrix[1, 2] # Covariance between the two

# Create sequence of time points
time_seq <- seq(
  min(combined_data$age_of_metastasis_minus_40, na.rm = TRUE),
  max(combined_data$age_of_metastasis_minus_40, na.rm = TRUE),
  length.out = 100
)

# Calculate time-varying coefficient: β(t) = β_phs + β_tt * t
beta_t <- beta_phs + beta_tt * time_seq

# Calculate standard error of β(t)
# SE[β(t)] = sqrt(Var(β_phs) + t^2 * Var(β_tt) + 2*t*Cov(β_phs, β_tt))
se_beta_t <- sqrt(
  vcov_matrix[1, 1] + time_seq^2 * vcov_matrix[2, 2] + 2 * time_seq * cov_phs_tt
)

# Calculate 95% confidence intervals
z_value <- qnorm(0.975) # Cox models use normal distribution
ci_lower <- beta_t - z_value * se_beta_t
ci_upper <- beta_t + z_value * se_beta_t

# Create plotting data frame
plot_data <- data.frame(
  time = time_seq + 40,
  beta_t = beta_t,
  ci_lower = ci_lower,
  ci_upper = ci_upper
)

# Create the plot
zph_tt <- ggplot(plot_data, aes(x = time, y = beta_t)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),
    fill = "steelblue", alpha = 0.3
  ) +
  geom_line(color = "steelblue", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_point(data = resid_data, aes(x = time, y = beta)) +
  stat_smooth() +
  geom_hline(
    yintercept = coef_phs,
    linetype = "dashed",
    color = "blue"
  ) +
  labs(
    title = "Schoenfeld Residuals for phs",
    subtitle = paste0(
      "Test for PH: p = ",
      round(zph_fit$table["phs", "p"], 4)
    ),
    x = "Time (age)",
    y = "Beta(t) for phs"
  ) +
  theme_classic()
saveRDS(zph_tt, paste(output_prefix, "zph_tt_aggr.rds", sep = "_"))

#########################################################

# KM Curves for metastatic disease
metastatic_cuminc <- phs_cuminc_strata(
  combined_data,
  time      = "age_of_metastasis",
  event     = "metastasis",
  intervals = list(c(0.95, 1), c(0.80, 1), c(0.30, 0.70), c(0, 0.20))
)
metastatic_km_curve <- plot_cuminc_strata(metastatic_cuminc)

# Build numbers-at-risk tables
any_pca_risk_table <- plot_risk_table(
  phs_risk_at_times(any_pca_cuminc, c(50, 60, 70, 80, 90)),
  x_range = range(any_pca_cuminc$time)
)
metastatic_risk_table <- plot_risk_table(
  phs_risk_at_times(metastatic_cuminc, c(50, 60, 70, 80, 90)),
  x_range = range(metastatic_cuminc$time)
)

# Cumulative incidence (1 - KM) by a predefined age for each PHS stratum, for Table 2
risk_age <- 70
cuminc_at_age <- function(cuminc_df, age, model) {
  cuminc_df %>%
    group_by(stratum) %>%
    summarise(
      last = max(which(time <= age)),
      estimate = risk_estimate[last], conf_low = conf_low[last],
      conf_high = conf_high[last], .groups = "drop"
    ) %>%
    transmute(cohort = output_prefix, model = model, age = age,
              stratum = as.character(stratum), estimate, conf_low, conf_high)
}
cuminc_age_table <- bind_rows(
  cuminc_at_age(any_pca_cuminc, risk_age, "any_pca"),
  cuminc_at_age(metastatic_cuminc, risk_age, "metastatic")
)
write_tsv(cuminc_age_table, paste(output_prefix, "cuminc_at_age.tsv", sep = "_"))

# Save Output
saveRDS(any_pca_km_curve,       paste(output_prefix, "any_pca_abs_risk.rds",       sep = "_"))
saveRDS(any_pca_risk_table,     paste(output_prefix, "any_pca_risk_table.rds",     sep = "_"))
saveRDS(metastatic_km_curve,    paste(output_prefix, "metastatic_abs_risk.rds",    sep = "_"))
saveRDS(metastatic_risk_table,  paste(output_prefix, "metastatic_risk_table.rds",  sep = "_"))
