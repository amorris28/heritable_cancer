library(tidyverse)
library(survival)
library(survminer)
library(broom)

set.seed(1984)
n = 2000

my_tibble = tibble(
  age = round(runif(n, 40, 100)),
  phs = rnorm(n),
  status = rbinom(n, 1, 0.4 + (0.1 * phs))
  )
# 
# ####
# 
# RK_get_regPHS(my_tibble$age, my_tibble$status, my_tibble$phs)
# 
# ####

ggplot(my_tibble, aes(x = phs, group = status, fill = status)) +
  geom_histogram()

agevals = 40:100
reference_PHS = my_tibble$phs
prctile_list = c(0, 0.2, 0.7, 0.8, 0.98, 1)
prctile_labels = c('0-20th', '20-70th', '70-80th', '80-98th', '98-100th')
my_tibble$percentile = 
  with(my_tibble,
       cut(phs, breaks = quantile(phs, probs = prctile_list),
                                      labels = prctile_labels,
                                      include.lowest = TRUE))


cox_model = coxph(Surv(age, status) ~ phs, data = my_tibble)
bt = cox_model$coefficients[["phs"]]

my_tibble %>% 
  filter(percentile %in% c('0-20th', '20-70th', '80-98th', '98-100th')) ->
  my_tibble

my_tibble %>% 
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

mod = survfit(Surv(age, status) ~ percentile, data = my_tibble)

mod_data = tidy(mod)

tibble(strata = unique(mod_data$strata),
       percentile = unique(S_tibble$percentile)) ->
  key

mod_data %>% 
  left_join(key, by = 'strata') ->
  mod_data

ggplot(mod_data, aes(x = time, y = estimate, color = percentile)) +
  geom_step() +
  theme_minimal() +
  xlim(40, 100) + 
  geom_line(aes(x = agevals, y = value, color = percentile), S_tibble) + 
  scale_color_brewer(palette = "Set1") +
  labs(x = "Age", y = "PCa-free Survival")

## Hazard ratios
 
perform = with(my_tibble, RK_get_perf(phs, age, status, ref = phs))

