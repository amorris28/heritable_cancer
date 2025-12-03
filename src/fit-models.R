library(tidyverse)
library(stringr)
library(hazrd)
library(survival)
library(patchwork)
source("src/config.R")
source("src/functions.R")

args = commandArgs(trailingOnly=TRUE)

metadata_file = args[1] # Phenotype data
perform_bool = args[2]
output_prefix = args[3]
#########################################
# Import PHS data and metadata
combined_data = read_tsv(metadata_file)

#########################################
# Analyze PHS data

# Coxph
combined_data %>%
  filter(status == 1) %>%
  summarize(min_age = min(age)) # 40ish

combined_data %>%
  filter(age >= 40) %>%
  mutate(age_minus_40 = age - 40) ->
  combined_data

fit = coxph(Surv(age, status) ~ phs, data = combined_data)
zph_fit = cox.zph(fit)

#plot(zph_fit)
#summary(fit)

zph_fit %>%
  print() %>%
  as_tibble %>%
  mutate(cohort = output_prefix, PCa = "any") ->
phs_zph_tidy
phs_zph_tidy = phs_zph_tidy[1,]

any_pca_coxph = broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)

# Time transformation
phs_fit_tt = coxph(Surv(age_minus_40, status) ~ phs + tt(phs), combined_data,
                  tt = function(x, t, ...) x * t)
phs_tt_tidy = broom::tidy(phs_fit_tt, exponentiate = TRUE, conf.int = TRUE)

# Performance metrics

if (perform_bool) {
  any_pca_perform = get_performance_metrics(combined_data, bootstrap_iters, output_prefix)
}

# Time to metastatic diseaes
metastasis_fit = coxph(Surv(age_of_metastasis, metastasis) ~ phs, combined_data)
metastatic_coxph = broom::tidy(metastasis_fit, exponentiate = TRUE, conf.int = TRUE)
meta_zph_fit = cox.zph(metastasis_fit)


meta_zph_fit %>%
  print() %>%
  as_tibble %>%
  mutate(cohort = output_prefix, PCa = "aggr") ->
meta_phs_zph_tidy
meta_phs_zph_tidy = meta_phs_zph_tidy[1,]
# Performance metrics

if (perform_bool) {
  metastatic_perform = get_performance_metrics(combined_data, bootstrap_iters, output_prefix, age = "age_of_metastasis", status = "metastasis")
}


# Gather sample sizes
combined_data %>%
  select(phs, age, status) ->
  any_pca_sample_size
any_pca_sample_size_cases = sum(any_pca_sample_size$status)
any_pca_sample_size_total = nrow(any_pca_sample_size)
any_pca_sample_size_controls = any_pca_sample_size_total - any_pca_sample_size_cases


combined_data %>%
  select(phs, age_of_metastasis, metastasis) -> #%>%
  # drop_na() ->
  metastatic_sample_size
metastatic_sample_size_cases = sum(metastatic_sample_size$metastasis)
metastatic_sample_size_total = nrow(metastatic_sample_size)
metastatic_sample_size_controls = 
  metastatic_sample_size_total - 
  metastatic_sample_size_cases

sample_size = tibble(group = c("any_pca", "metastatic"),
       total = c(any_pca_sample_size_total, metastatic_sample_size_total),
       cases = c(any_pca_sample_size_cases, metastatic_sample_size_cases),
       controls = c(any_pca_sample_size_controls, metastatic_sample_size_controls))

# Combine Coxph results
any_pca_coxph$model = "any_pca"
metastatic_coxph$model = "metastatic"
phs_tt_tidy$model = "phs_tt"

rbind(any_pca_coxph, 
      metastatic_coxph, 
      phs_tt_tidy) ->
      coxph_results

# Combine zph results
phs_zph_tidy %>%
  bind_rows(meta_phs_zph_tidy) ->
  phs_zph_tidy

# Save Output
write_tsv(phs_zph_tidy, paste(output_prefix, "phs_zph_tidy.tsv", sep = "_"))
write_tsv(any_pca_perform, paste(output_prefix, "any_pca_performance.tsv", sep = "_"))
write_tsv(metastatic_perform, paste(output_prefix, "metastatic_performance.tsv", sep = "_"))
write_tsv(sample_size, paste(output_prefix, "sample_size.tsv", sep = "_"))
write_tsv(coxph_results, paste(output_prefix, "coxph.tsv", sep = "_"))

##################################################################
# case-only local to met analysis
combined_data %>%
  filter(status == 1) ->
  local_to_met_data

fit = coxph(Surv(conversion_to_metastasis_years, metastasis) ~ phs, data = local_to_met_data)
summary(fit)
##################################################################
# case-only birth to met analysis
combined_data %>%
  filter(status == 1) ->
  birth_to_met_data

fit = coxph(Surv(age_of_metastasis, metastasis) ~ phs, data = birth_to_met_data)
summary(fit)
