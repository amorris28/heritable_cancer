library(tidyverse)
library(stringr)
library(hazrd)
library(survival)
library(patchwork)
source("src/config.R")
source("src/functions.R")

args = commandArgs(trailingOnly=TRUE)

metadata_file = args[1] # Phenotype data
output_prefix = args[2]

#########################################
# Import PHS data and metadata
combined_data = read_tsv(metadata_file)

#########################################
# Filter over 40 
combined_data %>%
  filter(status == 1) %>%
  summarize(min_age = min(age)) # 40ish

combined_data %>%
  filter(age >= 40) %>%
  mutate(age_minus_40 = age - 40) %>%
  mutate(age_of_metastasis_minus_40 = age_of_metastasis - 40) ->
  combined_data

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
phs_fit_tt = coxph(Surv(age_minus_40, status) ~ phs + tt(phs), combined_data,
                  tt = function(x, t, ...) x * t)


# Extract coefficients and variance-covariance matrix
coefs <- coef(phs_fit_tt)
vcov_matrix <- vcov(phs_fit_tt)

# Coefficients
beta_phs <- coefs[1]        # Main effect
beta_tt <- coefs[2]         # Time interaction effect

# Standard errors
se_phs <- sqrt(vcov_matrix[1, 1])
se_tt <- sqrt(vcov_matrix[2, 2])
cov_phs_tt <- vcov_matrix[1, 2]  # Covariance between the two

# Create sequence of time points
time_seq <- seq(min(combined_data$age_minus_40, na.rm = TRUE), 
                max(combined_data$age_minus_40, na.rm = TRUE), 
                length.out = 100)

# Calculate time-varying coefficient: β(t) = β_phs + β_tt * t
beta_t <- beta_phs + beta_tt * time_seq

# Calculate standard error of β(t)
# SE[β(t)] = sqrt(Var(β_phs) + t^2 * Var(β_tt) + 2*t*Cov(β_phs, β_tt))
se_beta_t <- sqrt(vcov_matrix[1, 1] + 
                  time_seq^2 * vcov_matrix[2, 2] + 
                  2 * time_seq * cov_phs_tt)

# Calculate 95% confidence intervals
z_value <- qnorm(0.975)  # Cox models use normal distribution
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
zph_tt = ggplot(plot_data, aes(x = time, y = beta_t)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),
              fill = "steelblue", alpha = 0.3) +
  geom_line(color = "steelblue", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  geom_point(data = resid_data, aes(x = time, y = beta)) +
  stat_smooth() +
  geom_hline(yintercept = coef_phs, 
             linetype = "dashed", 
             color = "blue") +
  labs(title = "Schoenfeld Residuals for phs",
       subtitle = paste0("Test for PH: p = ", 
                        round(zph_fit$table["phs", "p"], 4)),
       x = "Time (age)",
       y = "Beta(t) for phs") +
  theme_classic()
saveRDS(zph_tt, paste(output_prefix, "zph_tt.rds", sep = "_"))

#########################################################

# Plot KM Curves
any_pca_km_curve = plot_km_curves(combined_data)
any_pca_inc_curve = plot_km_curves(combined_data, style = "risk")


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
phs_fit_tt = coxph(Surv(age_of_metastasis_minus_40, metastasis) ~ phs + tt(phs), combined_data,
                  tt = function(x, t, ...) x * t)


# Extract coefficients and variance-covariance matrix
coefs <- coef(phs_fit_tt)
vcov_matrix <- vcov(phs_fit_tt)

# Coefficients
beta_phs <- coefs[1]        # Main effect
beta_tt <- coefs[2]         # Time interaction effect

# Standard errors
se_phs <- sqrt(vcov_matrix[1, 1])
se_tt <- sqrt(vcov_matrix[2, 2])
cov_phs_tt <- vcov_matrix[1, 2]  # Covariance between the two

# Create sequence of time points
time_seq <- seq(min(combined_data$age_of_metastasis_minus_40, na.rm = TRUE), 
                max(combined_data$age_of_metastasis_minus_40, na.rm = TRUE), 
                length.out = 100)

# Calculate time-varying coefficient: β(t) = β_phs + β_tt * t
beta_t <- beta_phs + beta_tt * time_seq

# Calculate standard error of β(t)
# SE[β(t)] = sqrt(Var(β_phs) + t^2 * Var(β_tt) + 2*t*Cov(β_phs, β_tt))
se_beta_t <- sqrt(vcov_matrix[1, 1] + 
                  time_seq^2 * vcov_matrix[2, 2] + 
                  2 * time_seq * cov_phs_tt)

# Calculate 95% confidence intervals
z_value <- qnorm(0.975)  # Cox models use normal distribution
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
zph_tt = ggplot(plot_data, aes(x = time, y = beta_t)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),
              fill = "steelblue", alpha = 0.3) +
  geom_line(color = "steelblue", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  geom_point(data = resid_data, aes(x = time, y = beta)) +
  stat_smooth() +
  geom_hline(yintercept = coef_phs, 
             linetype = "dashed", 
             color = "blue") +
  labs(title = "Schoenfeld Residuals for phs",
       subtitle = paste0("Test for PH: p = ", 
                        round(zph_fit$table["phs", "p"], 4)),
       x = "Time (age)",
       y = "Beta(t) for phs") +
  theme_classic()
saveRDS(zph_tt, paste(output_prefix, "zph_tt_aggr.rds", sep = "_"))

#########################################################

# KM Curves for metastatic disease
metastatic_km_curve = plot_km_curves(combined_data, age = "age_of_metastasis", status = "metastasis")
metastatic_inc_curve = plot_km_curves(combined_data, age = "age_of_metastasis", status = "metastasis", style = "risk")

# Save Output
saveRDS(any_pca_km_curve, paste(output_prefix, "any_pca_km_curve.rds", sep = "_"))
saveRDS(metastatic_km_curve, paste(output_prefix, "metastatic_km_curve.rds", sep = "_"))
saveRDS(any_pca_inc_curve, paste(output_prefix, "any_pca_km_curve.rds", sep = "_"))
saveRDS(metastatic_inc_curve, paste(output_prefix, "metastatic_km_curve.rds", sep = "_"))
