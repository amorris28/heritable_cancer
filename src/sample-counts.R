library(tidyverse)
source("src/functions.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop(
    "Usage: sample-counts.R <phenotypes.tsv> <p2521|p212> ",
    "<0deg.keep> <2deg.keep> <3deg.keep> <sample_counts.tsv>"
  )
}

phenotypes_file <- args[1]
cohort <- args[2]
relatedness_files <- args[3:5]
sample_counts_file <- args[6]

if (!cohort %in% c("p2521", "p212")) {
  stop("Cohort must be p2521 or p212; got: ", cohort)
}

phenotypes <- read_tsv(phenotypes_file, show_col_types = FALSE)
required_columns <- c("bfile_id", "age", "status", "phs")
missing_columns <- setdiff(required_columns, names(phenotypes))
if (length(missing_columns) > 0) {
  stop("Missing required phenotype columns: ", paste(missing_columns, collapse = ", "))
}

analysis_data <- phenotypes %>%
  filter(age >= 40)
if (nrow(analysis_data) == 0) {
  stop("No age-40+ participants found in ", phenotypes_file)
}
if (anyDuplicated(as.character(analysis_data$bfile_id))) {
  stop("Duplicate bfile_id values found in the age-40+ phenotype data")
}

sample_ids <- as.character(analysis_data$bfile_id)
sample_count_file <- read_tsv(sample_counts_file, show_col_types = FALSE)
preprocessed_age40_n <- sample_count_file %>%
  filter(stage == "male_with_phs_and_complete_phenotype_age40plus") %>%
  pull(n)
if (length(preprocessed_age40_n) != 1 || preprocessed_age40_n != length(sample_ids)) {
  stop(
    "Age-40+ phenotype count does not match preprocessing report: ",
    length(sample_ids), " vs ", paste(preprocessed_age40_n, collapse = ", ")
  )
}
labels <- c("0deg", "2deg", "3deg")
retained_ids <- purrr::map(relatedness_files, read_keep_ids, sample_ids)
names(retained_ids) <- labels
retained_n <- purrr::map_int(retained_ids, ~ length(intersect(sample_ids, .x)))
pre_king_n <- length(sample_ids)

counts <- sample_count_rows(
  cohort,
  stage = c(
    "age40plus_before_king",
    paste0("age40plus_after_king_", labels),
    "age40plus_removed_by_king_0deg"
  ),
  n = c(
    pre_king_n,
    unname(retained_n),
    pre_king_n - unname(retained_n[["0deg"]])
  ),
  definition = c(
    "Male participants with PHS, age and PCa status, aged 40 years or older; primary model input before KING relatedness filtering.",
    paste0(
      "Age-40+ model-input participants retained by KING ", labels,
      " list; before model-specific outcome exclusions."
    ),
    "Participants removed from the age-40+ primary model input by the KING 0-degree list (duplicate/MZ-level); difference from before-KING count."
  )
)
record_sample_counts(counts, sample_counts_file)

message("Sample counts written to ", sample_counts_file)
