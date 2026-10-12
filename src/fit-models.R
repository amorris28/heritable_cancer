library(tidyverse)
library(stringr)
library(survival)
library(patchwork)
source("src/functions.R")

args <- commandArgs(trailingOnly = TRUE)

metadata_file <- args[1] # Phenotype data
perform_bool <- as.logical(args[2])
output_prefix <- args[3]
bootstrap_iters <- as.integer(args[4])
relatedness_keep_files <- if (length(args) >= 5) args[5:length(args)] else character()

if (length(relatedness_keep_files) > 0 && length(relatedness_keep_files) != 3) {
  stop("Provide all three relatedness keep files: 0deg, 2deg, and 3deg")
}

read_relatedness_ids <- function(path, sample_ids) {
  sample_ids <- as.character(sample_ids)
  if (!file.exists(path)) {
    stop("Relatedness keep file does not exist: ", path)
  }
  id_table <- readr::read_table(
    path,
    col_names = FALSE,
    col_types = readr::cols(.default = readr::col_character()),
    comment = "#",
    progress = FALSE
  )
  if (nrow(id_table) == 0 || ncol(id_table) == 0) {
    stop("Relatedness keep file is empty: ", path)
  }
  id_table <- id_table %>%
    filter(!if_all(everything(), ~ toupper(.x) %in% c("FID", "IID")))
  if (nrow(id_table) == 0) {
    stop("Relatedness keep file contains no sample IDs: ", path)
  }
  overlaps <- vapply(id_table, function(column) {
    length(intersect(sample_ids, column))
  }, integer(1))
  if (max(overlaps) == 0) {
    stop("No IDs in ", path, " match phenotype bfile_id values")
  }
  ids <- unique(id_table[[which.max(overlaps)]])
  message(
    "Read ", length(ids), " retained IDs from ", path,
    "; ", max(overlaps), " match the phenotype cohort"
  )
  ids
}


#########################################
# Import PHS data and metadata
combined_data <- read_tsv(metadata_file)

#########################################
# Analyze PHS data

# Coxph, min case age is 40ish
analysis_data <-
  combined_data %>%
  filter(age >= 40) %>%
  mutate(age_minus_40 = age - 40)

# The cutoff-0 list removes duplicate/monozygotic-level samples while retaining
# one sample per pair. Use it for the primary p2521 analyses.
if (length(relatedness_keep_files) == 3) {
  primary_ids <- read_relatedness_ids(relatedness_keep_files[[1]], analysis_data$bfile_id)
  before_n <- nrow(analysis_data)
  combined_data <- analysis_data %>%
    filter(as.character(bfile_id) %in% primary_ids)
  if (nrow(combined_data) == 0) {
    stop("The 0-degree relatedness list removed every analysis participant")
  }
  message(
    "Primary 0-degree set: retained ", nrow(combined_data), " of ", before_n,
    " participants after age eligibility"
  )
} else {
  combined_data <- analysis_data
}

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
  any_pca_perform <- phs_performance_metrics(
    combined_data,
    n_boot = bootstrap_iters,
    seed   = 42
  )
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
  metastatic_perform <- phs_performance_metrics(
    combined_data,
    time   = "age_of_metastasis",
    event  = "metastasis",
    n_boot = bootstrap_iters,
    seed   = 42
  )
}

# PCa-specific mortality
pca_death_fit <- coxph(
  Surv(age_at_pca_death, pca_death) ~ phs,
  combined_data
)
pca_death_coxph <- broom::tidy(pca_death_fit,
  exponentiate = TRUE,
  conf.int = TRUE
)
pca_death_zph_fit <- cox.zph(pca_death_fit)

pca_death_zph_tidy <-
  pca_death_zph_fit %>%
  print() %>%
  as_tibble() %>%
  mutate(cohort = output_prefix, PCa = "pca_death")
pca_death_zph_tidy <- pca_death_zph_tidy[1, ]

# Performance metrics
if (perform_bool) {
  pca_death_perform <- phs_performance_metrics(
    combined_data,
    time   = "age_at_pca_death",
    event  = "pca_death",
    n_boot = bootstrap_iters,
    seed   = 42
  )
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

pca_death_sample_size <-
  combined_data %>%
  select(phs, age_at_pca_death, pca_death)
pca_death_sample_size_cases <- sum(pca_death_sample_size$pca_death)
pca_death_sample_size_total <- nrow(pca_death_sample_size)
pca_death_sample_size_controls <-
  pca_death_sample_size_total - pca_death_sample_size_cases

sample_size <- tibble(
  group = c("any_pca", "metastatic", "pca_death"),
  total = c(
    any_pca_sample_size_total,
    metastatic_sample_size_total,
    pca_death_sample_size_total
  ),
  cases = c(
    any_pca_sample_size_cases,
    metastatic_sample_size_cases,
    pca_death_sample_size_cases
  ),
  controls = c(
    any_pca_sample_size_controls,
    metastatic_sample_size_controls,
    pca_death_sample_size_controls
  )
)

# Estimate the same three PHS coefficients in increasingly unrelated subsets.
relatedness_sensitivity <- NULL
if (length(relatedness_keep_files) == 3) {
  sensitivity_specs <- list(
    any_pca = c(time = "age", event = "status"),
    metastatic = c(time = "age_of_metastasis", event = "metastasis"),
    pca_death = c(time = "age_at_pca_death", event = "pca_death")
  )
  sensitivity_labels <- c("0deg", "2deg", "3deg")
  sensitivity_cutoffs <- c("0deg" = 0.35333, "2deg" = 0.0883883, "3deg" = 0.0441942)
  relatedness_sensitivity <- purrr::map2_dfr(
    relatedness_keep_files, sensitivity_labels,
    function(keep_file, degree_label) {
      keep_ids <- read_relatedness_ids(keep_file, analysis_data$bfile_id)
      sensitivity_data <- analysis_data %>%
        filter(as.character(bfile_id) %in% keep_ids)
      purrr::imap_dfr(sensitivity_specs, function(spec, model_name) {
        model_data <- sensitivity_data %>%
          filter(!is.na(.data[[spec[["time"]]]]), !is.na(.data[[spec[["event"]]]]))
        if (nrow(model_data) == 0 || sum(model_data[[spec[["event"]]]]) == 0) {
          stop(
            "No analyzable events for ", model_name, " in ",
            degree_label, " relatedness subset"
          )
        }
        formula <- stats::as.formula(paste0(
          "survival::Surv(", spec[["time"]], ", ", spec[["event"]], ") ~ phs"
        ))
        fit <- survival::coxph(formula, data = model_data)
        broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
          mutate(
            model = paste0(model_name, "_", degree_label),
            relatedness_set = degree_label,
            kinship_cutoff = unname(sensitivity_cutoffs[degree_label]),
            n = stats::nobs(fit),
            events = fit$nevent
          )
      })
    }
  )
}

# Combine Coxph results
any_pca_coxph$model <- "any_pca"
metastatic_coxph$model <- "metastatic"
pca_death_coxph$model <- "pca_death"
phs_tt_tidy$model <- "phs_tt"

coxph_results <-
  rbind(
    any_pca_coxph,
    metastatic_coxph,
    pca_death_coxph,
    phs_tt_tidy
  )

# Combine zph results
phs_zph_tidy <-
  phs_zph_tidy %>%
  bind_rows(meta_phs_zph_tidy) %>%
  bind_rows(pca_death_zph_tidy)

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
  pca_death_perform,
  paste(output_prefix, "pca_death_performance.tsv", sep = "_")
)
write_tsv(
  sample_size,
  paste(output_prefix, "sample_size.tsv", sep = "_")
)
write_tsv(coxph_results, paste(output_prefix, "coxph.tsv", sep = "_"))
if (!is.null(relatedness_sensitivity)) {
  write_tsv(
    relatedness_sensitivity,
    paste(output_prefix, "relatedness_sensitivity.tsv", sep = "_")
  )
}
