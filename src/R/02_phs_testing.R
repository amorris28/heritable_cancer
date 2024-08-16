library(tidyverse)
library(survival)
library(survminer)
library(broom)

#########################################
# First argument: sumstats file
# Second argument: metadata file

args = commandArgs(trailingOnly=TRUE)

model = args[1] # Sumstats file

metadata = read_tsv(args[2]) # PHS data


#########################################
# Import PHS data

phs = read_tsv(model, col_names = c("bfile_id", "phs"))
phs$phs = phs$phs * -1
metadata %>%
  left_join(phs, by = "bfile_id") %>%
  mutate(status = factor) %>%
  filter(!(is.na(phs))) ->
  combined_data

phs_hist = ggplot(combined_data, aes(x = phs, group = status, fill = status)) +
  geom_histogram()

ggsave(paste0("figs/", model, "_hist.png"), phs_hist)

####################################
# Generate data
set.seed(70454416)
n = 2000
combined_data = tibble(
  age = round(runif(n, 40, 100)),
  phs = rnorm(n),
  status = rbinom(n, 1, 0.4 + (0.1 * phs))
  )

ggplot(combined_data, aes(x = phs, group = status, fill = status)) +
  geom_histogram()

######################
# Fit model

agevals = 40:100
reference_PHS = combined_data$phs
prctile_list = c(0, 0.2, 0.7, 0.8, 0.98, 1)
prctile_labels = c('0-20th', '20-70th', '70-80th', '80-98th', '98-100th')
combined_data$percentile = 
  with(combined_data,
       cut(phs, breaks = quantile(phs, probs = prctile_list),
                                      labels = prctile_labels,
                                      include.lowest = TRUE))

cox_model = coxph(Surv(age, status) ~ phs, data = combined_data)

bt = cox_model$coefficients[["phs"]]

combined_data %>% 
  filter(percentile %in% c('0-20th', '20-70th', '80-98th', '98-100th')) ->
  combined_data

combined_data %>% 
  group_by(percentile) %>% 
  summarize(median_phs = median(phs)) ->
  median_tibble


median_phs_individual = quantile(reference_PHS, 0.50)
critvalvec = median_tibble$median_phs

#########################
# Calculate idealized curves
#
# Basically copied this from:
# https://github.com/cmig-research-group/phs/blob/master/phs_model_example/phs_generate_new_demo.m
a = 0.0700 / 100;
b = 0.0753
pcainc = a * exp(b * (agevals - 40)) # this is the population incidence rate for PCa
H0 = a / b * exp(b * (agevals - 40)) # compute baseline hazard
H0 = H0 - H0[1] # assume cumulative hazard == 0 before age 40

H = matrix(nrow = length(H0), ncol = length(critvalvec))
S = matrix(nrow = length(H0), ncol = length(critvalvec))
for (cnt in 1:length(critvalvec)) {
  H[, cnt] = H0 * exp(bt * critvalvec[cnt] - bt * median_phs_individual[[1]])
  S[, cnt] = exp(-H[, cnt])
}

S_tibble = as_tibble(S, .name_repair = ~ as.character(median_tibble$percentile))
S_tibble$agevals = agevals
S_tibble %>%
  pivot_longer(-agevals, names_to = "percentile") ->
  S_tibble
##########################

# Plot K-M curves for centiles
mod = survfit(Surv(age, status) ~ percentile, data = combined_data)
mod_data = tidy(mod)

tibble(strata = unique(mod_data$strata),
       percentile = unique(S_tibble$percentile)) ->
  key

mod_data %>% 
  left_join(key, by = 'strata') ->
  mod_data

kmcurve = ggplot(mod_data, aes(x = time, y = estimate, color = percentile)) +
  geom_step() +
  theme_minimal() +
  xlim(40, 100) + 
  geom_line(aes(x = agevals, y = value, color = percentile), S_tibble) + 
  scale_color_brewer(palette = "Set1") +
  labs(x = "Age", y = "PCa-free Survival")

ggsave(paste0("figs/", model, "_km_curve.png"),
       kmcurve)

combined_data %>% 
  group_by(status) %>% 
  summarize(n = n()) ->
  cases_controls

## Hazard ratios
source("R/RK_get_perf.R")
perform = with(combined_data, RK_get_perf(phs, age, status, ref = phs, swc.popnumcases = cases_controls[[2, 2]], swc.popnumcontrols = cases_controls[[1, 2]]))
perform = stack(perform)[, c(2, 1)]
write_tsv(perform, paste0(model, "_performance.txt"))
