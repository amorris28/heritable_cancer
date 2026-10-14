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
relatedness_keep_files <- if (length(args) >= 7) args[5:7] else character()
pca_file <- if (length(args) >= 8) args[8] else NULL

if (length(args) != 8) {
  stop(
    "Usage: fit-models.R <phenotypes.tsv> <perform> <output_prefix> ",
    "<bootstrap_iters> <0deg.keep> <2deg.keep> <3deg.keep> <pcs.txt>"
  )
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

read_pca_scores <- function(path, sample_ids) {
  if (!file.exists(path)) {
    stop("PCA scores file does not exist: ", path)
  }
  scores <- readr::read_table(
    path,
    col_types = readr::cols(.default = readr::col_character()),
    progress = FALSE
  )
  normalized_names <- toupper(sub("^#", "", names(scores)))
  id_index <- match("IID", normalized_names)
  if (is.na(id_index)) {
    id_index <- match("ID", normalized_names)
  }
  if (is.na(id_index)) {
    stop("Could not find an IID or ID column in PCA scores file: ", path)
  }
  pc_indices <- which(grepl("^PC[0-9]+$", normalized_names))
  if (length(pc_indices) < 10) {
    stop(
      "Expected at least 10 PC columns in ", path, "; found ",
      length(pc_indices)
    )
  }
  pc_numbers <- as.integer(sub("^PC", "", normalized_names[pc_indices]))
  pc_indices <- pc_indices[order(pc_numbers)]
  pc_numbers <- sort(pc_numbers)
  if (anyDuplicated(pc_numbers) || !all(seq_len(10) %in% pc_numbers)) {
    stop("PCA scores must contain PC1 through PC10: ", path)
  }
  pca_ids <- as.character(scores[[id_index]])
  if (anyNA(pca_ids) || anyDuplicated(pca_ids)) {
    stop("PCA scores have missing or duplicate sample IDs: ", path)
  }
  missing_ids <- setdiff(as.character(sample_ids), pca_ids)
  if (length(missing_ids) > 0) {
    stop(
      length(missing_ids), " phenotype sample IDs have no PCA scores in ", path,
      "; examples: ", paste(head(missing_ids, 5), collapse = ", ")
    )
  }
  relevant_rows <- match(as.character(sample_ids), pca_ids)
  pca_data <- tibble::tibble(bfile_id = as.character(sample_ids))
  for (pc_number in seq_len(10)) {
    column_index <- pc_indices[match(pc_number, pc_numbers)]
    pc_values <- suppressWarnings(as.numeric(scores[[column_index]][relevant_rows]))
    if (anyNA(pc_values)) {
      stop(
        "PCA column ", names(scores)[column_index],
        " contains missing or non-numeric values for phenotype samples"
      )
    }
    pca_data[[paste0("PC", pc_number)]] <- pc_values
  }
  pca_data
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

# Preserve the original all-subject, unadjusted models as the primary analysis.
combined_data <- analysis_data

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

# Compare PHS coefficients across all subjects and KING retained sets, with no
# PCs and with the first 5 or 10 PCs as covariates.
sensitivity_specs <- list(
  any_pca = c(time = "age", event = "status"),
  metastatic = c(time = "age_of_metastasis", event = "metastasis"),
  pca_death = c(time = "age_at_pca_death", event = "pca_death")
)
relatedness_labels <- c("0deg", "2deg", "3deg")
relatedness_cutoffs <- c("0deg" = 0.35333, "2deg" = 0.0883883, "3deg" = 0.0441942)
relatedness_ids <- purrr::map2(
  relatedness_keep_files, relatedness_labels,
  ~ read_relatedness_ids(.x, analysis_data$bfile_id)
)
names(relatedness_ids) <- relatedness_labels
pca_scores <- read_pca_scores(pca_file, analysis_data$bfile_id)
analysis_data <- analysis_data %>%
  mutate(bfile_id = as.character(bfile_id)) %>%
  left_join(pca_scores, by = "bfile_id")

relatedness_sets <- c(
  list(all = list(data = analysis_data, cutoff = NA_real_)),
  purrr::imap(relatedness_ids, function(ids, label) {
    list(
      data = analysis_data %>% filter(bfile_id %in% ids),
      cutoff = unname(relatedness_cutoffs[[label]])
    )
  })
)
names(relatedness_sets) <- c("all", relatedness_labels)

genetic_qc_sensitivity <- purrr::imap_dfr(
  relatedness_sets,
  function(set, relatedness_label) {
    purrr::imap_dfr(sensitivity_specs, function(spec, model_name) {
      purrr::map_dfr(c(0L, 5L, 10L), function(pc_count) {
        model_data <- set$data %>%
          filter(!is.na(.data[[spec[["time"]]]]), !is.na(.data[[spec[["event"]]]]))
        pc_terms <- if (pc_count == 0L) character() else paste0("PC", seq_len(pc_count))
        if (length(pc_terms) > 0) {
          model_data <- model_data %>%
            filter(if_all(all_of(pc_terms), ~ !is.na(.x)))
        }
        if (
          nrow(model_data) == 0 ||
          sum(model_data[[spec[["event"]]]], na.rm = TRUE) == 0
        ) {
          stop(
            "No analyzable events for ", model_name, " in ",
            relatedness_label, " relatedness set with ", pc_count, " PCs"
          )
        }
        formula <- stats::as.formula(paste(
          "survival::Surv(", spec[["time"]], ", ", spec[["event"]], ") ~ phs",
          if (length(pc_terms) > 0) paste("+", paste(pc_terms, collapse = " + ")) else ""
        ))
        fit <- survival::coxph(formula, data = model_data)
        broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
          filter(term == "phs") %>%
          mutate(
            model = model_name,
            relatedness_set = relatedness_label,
            kinship_cutoff = set$cutoff,
            pc_count = pc_count,
            n = fit$n,
            events = fit$nevent
          )
      })
    })
  }
)

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
write_tsv(
  genetic_qc_sensitivity,
  paste(output_prefix, "genetic_qc_sensitivity.tsv", sep = "_")
)
