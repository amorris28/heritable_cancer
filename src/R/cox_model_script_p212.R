library(tidyverse)
library(survival)
library(survminer)
library(broom)

model = "PRACTICAL_PCa_Huynh-Le2022_phs_imputed_antegenes"

metadata = read_tsv('combined_data_simplified.tsv')

phs = read_tsv(model, col_names = c("bfile_id", "phs"))
phs$phs = phs$phs * -1
metadata %>% 
  left_join(phs, by = "bfile_id") %>% 
  filter(!(is.na(phs))) ->
  combined_data

p1 = ggplot(combined_data, aes(x = phs, group = status, fill = status)) +
  geom_histogram()

ggsave(paste0(model, "hist.png"), p1)

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

mod = survfit(Surv(age, status) ~ percentile, data = combined_data)

mod_data = tidy(mod)

tibble(strata = unique(mod_data$strata),
       percentile = unique(S_tibble$percentile)) ->
  key

mod_data %>% 
  left_join(key, by = 'strata') ->
  mod_data

km_curve = ggplot(mod_data, aes(x = time, y = estimate, color = percentile)) +
  geom_step() +
  theme_minimal() +
  xlim(40, 100) + 
  geom_line(aes(x = agevals, y = value, color = percentile), S_tibble) + 
  scale_color_brewer(palette = "Set1") +
  labs(x = "Age", y = "PCa-free Survival")

ggsave(paste0(model, "_km_curve.png"), km_curve)

## Hazard ratios
source('RK_get_perf.R')
perform = with(combined_data, RK_get_perf(phs, age, status, ref = phs))
