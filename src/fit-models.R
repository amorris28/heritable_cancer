library(tidyverse)
library(stringr)
library(hazrd)
library(survival)
library(patchwork)
source("src/functions.R")

args <- commandArgs(trailingOnly = TRUE)

metadata_file <- args[1] # Phenotype data
perform_bool <- args[2]
output_prefix <- args[3]
bootstrap_iters <- args[4]


#########################################
# Import PHS data and metadata
combined_data <- read_tsv(metadata_file)

#########################################
# Analyze PHS data

# Coxph, min case age is 40ish
combined_data <-
  combined_data %>%
  filter(age >= 40) %>%
  mutate(age_minus_40 = age - 40)

fit <- coxph(Surv(age, status) ~ phs, data = combined_data)
zph_fit <- cox.zph(fit)

phs_zph_tidy <-
  zph_fit %>%
  print() %>%
  as_tibble() %>%
  mutate(cohort = output_prefix, PCa = "any")

phs_zph_tidy <- phs_zph_tidy[1, ]

any_pca_coxph <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)

# Time transformation
phs_fit_tt <- coxph(Surv(age_minus_40, status) ~ phs + tt(phs), combined_data,
  tt = function(x, t, ...) x * t
)
phs_tt_tidy <- broom::tidy(phs_fit_tt, exponentiate = TRUE, conf.int = TRUE)

# Performance metrics
if (perform_bool) {
  any_pca_perform <-
    performance_metrics(combined_data, output_prefix)
}

# Time to metastatic diseaes
metastasis_fit <- coxph(
  Surv(age_of_metastasis, metastasis) ~ phs,
  combined_data
)
metastatic_coxph <- broom::tidy(metastasis_fit,
  exponentiate = TRUE,
  conf.int = TRUE
)
meta_zph_fit <- cox.zph(metastasis_fit)

meta_phs_zph_tidy <-
  meta_zph_fit %>%
  print() %>%
  as_tibble() %>%
  mutate(cohort = output_prefix, PCa = "aggr")
meta_phs_zph_tidy <- meta_phs_zph_tidy[1, ]

# Performance metrics
if (perform_bool) {
  metastatic_perform <-
      performance_metrics(
        combined_data,
        output_prefix,
        age = "age_of_metastasis",
        status = "metastasis")
}

# Gather sample sizes
any_pca_sample_size <-
  combined_data %>%
  select(phs, age, status)
any_pca_sample_size_cases <- sum(any_pca_sample_size$status)
any_pca_sample_size_total <- nrow(any_pca_sample_size)
any_pca_sample_size_controls <-
  any_pca_sample_size_total - any_pca_sample_size_cases

metastatic_sample_size <-
  combined_data %>%
  select(phs, age_of_metastasis, metastasis)
metastatic_sample_size_cases <- sum(metastatic_sample_size$metastasis)
metastatic_sample_size_total <- nrow(metastatic_sample_size)
metastatic_sample_size_controls <-
  metastatic_sample_size_total -
  metastatic_sample_size_cases

sample_size <- tibble(
  group = c("any_pca", "metastatic"),
  total = c(
    any_pca_sample_size_total,
    metastatic_sample_size_total
  ),
  cases = c(
    any_pca_sample_size_cases,
    metastatic_sample_size_cases
  ),
  controls = c(
    any_pca_sample_size_controls,
    metastatic_sample_size_controls
  )
)

# Combine Coxph results
any_pca_coxph$model <- "any_pca"
metastatic_coxph$model <- "metastatic"
phs_tt_tidy$model <- "phs_tt"

coxph_results <-
  rbind(
    any_pca_coxph,
    metastatic_coxph,
    phs_tt_tidy
  )

# Combine zph results
phs_zph_tidy <-
  phs_zph_tidy %>%
  bind_rows(meta_phs_zph_tidy)

# Save Output
write_tsv(phs_zph_tidy, paste(output_prefix, "phs_zph_tidy.tsv", sep = "_"))
write_tsv(
  any_pca_perform,
  paste(output_prefix, "any_pca_performance.tsv", sep = "_")
)
write_tsv(
  metastatic_perform,
  paste(output_prefix, "metastatic_performance.tsv", sep = "_")
)
write_tsv(
  sample_size,
  paste(output_prefix, "sample_size.tsv", sep = "_")
)
write_tsv(coxph_results, paste(output_prefix, "coxph.tsv", sep = "_"))
