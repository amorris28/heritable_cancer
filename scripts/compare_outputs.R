#!/usr/bin/env Rscript
# Dev utility: compare a cohort's freshly regenerated data/public/ outputs
# against a baseline saved by scripts/backup_outputs_for_diff.sh.
#
# Usage:
#   Rscript scripts/compare_outputs.R <p212|p2521> <tag>
#   Rscript scripts/compare_outputs.R manuscript <tag>
#
# For each file present in the baseline, reports:
#   - IDENTICAL       byte-for-byte or numerically identical
#   - CLOSE (tol)     .tsv numeric columns differ by <= tol (default 1e-6);
#                     expected for bootstrap-derived columns (conf_low,
#                     conf_high, se) since the refactor replaced hazrd's
#                     bootstrap internals with a boot::boot-based one — same
#                     seed, different resampling code path, so cannot be
#                     expected to match exactly even when correct.
#   - DIFFERS         otherwise; printed with enough detail to inspect by hand
#   - MISSING (new)   file not present in the freshly regenerated output
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript scripts/compare_outputs.R <p212|p2521|manuscript> <tag>")
}
mode <- args[1]
tag <- args[2]

# Columns known to come from bootstrap resampling (not expected to match
# exactly across the refactor, only approximately).
bootstrap_cols <- c("conf_low", "conf_high", "se")

if (mode == "manuscript") {
  current_dir <- "manuscript/figures_and_tables"
  baseline_dir <- file.path(current_dir, "_baseline", tag)
} else {
  current_dir <- "data/public"
  baseline_dir <- file.path(current_dir, "_baseline", tag)
}

if (!dir.exists(baseline_dir)) {
  stop("No baseline found at ", baseline_dir, " — run scripts/backup_outputs_for_diff.sh first")
}

baseline_files <- list.files(baseline_dir)
if (length(baseline_files) == 0) {
  stop("Baseline directory ", baseline_dir, " is empty")
}

compare_tsv <- function(old_path, new_path) {
  old <- suppressMessages(read_tsv(old_path, show_col_types = FALSE))
  new <- suppressMessages(read_tsv(new_path, show_col_types = FALSE))

  if (!identical(dim(old), dim(new)) || !identical(sort(names(old)), sort(names(new)))) {
    return(list(status = "DIFFERS", detail = sprintf(
      "shape/columns differ: old %s [%s] vs new %s [%s]",
      paste(dim(old), collapse = "x"), paste(names(old), collapse = ","),
      paste(dim(new), collapse = "x"), paste(names(new), collapse = ",")
    )))
  }

  new <- new[names(old)] # align column order
  mismatches <- c()
  approx_only <- TRUE
  for (col in names(old)) {
    o <- old[[col]]
    n <- new[[col]]
    if (is.numeric(o) && is.numeric(n)) {
      tol <- if (col %in% bootstrap_cols) 0.25 else 1e-6
      diffs <- abs(o - n)
      bad <- which(diffs > tol & !(is.na(o) & is.na(n)))
      if (length(bad) > 0) {
        mismatches <- c(mismatches, sprintf(
          "  col '%s': max|diff|=%.6g over %d row(s) (tol=%.3g)%s",
          col, max(diffs[bad], na.rm = TRUE), length(bad), tol,
          if (col %in% bootstrap_cols) "  [bootstrap-derived, approx expected]" else ""
        ))
        if (!(col %in% bootstrap_cols)) approx_only <- FALSE
      }
    } else if (!identical(o, n)) {
      mismatches <- c(mismatches, sprintf("  col '%s': non-numeric mismatch", col))
      approx_only <- FALSE
    }
  }

  if (length(mismatches) == 0) {
    list(status = "IDENTICAL", detail = NULL)
  } else if (approx_only) {
    list(status = "CLOSE (bootstrap cols only)", detail = paste(mismatches, collapse = "\n"))
  } else {
    list(status = "DIFFERS", detail = paste(mismatches, collapse = "\n"))
  }
}

compare_rds <- function(old_path, new_path) {
  old <- readRDS(old_path)
  new <- readRDS(new_path)
  if (isTRUE(all.equal(old, new, tolerance = 1e-8))) {
    list(status = "IDENTICAL", detail = NULL)
  } else {
    list(status = "DIFFERS", detail = paste(utils::capture.output(all.equal(old, new)), collapse = "\n"))
  }
}

results <- map_dfr(baseline_files, function(fname) {
  old_path <- file.path(baseline_dir, fname)
  new_path <- file.path(current_dir, fname)

  if (!file.exists(new_path)) {
    return(tibble(file = fname, status = "MISSING (new)", detail = NA_character_))
  }

  ext <- tools::file_ext(fname)
  res <- tryCatch(
    if (ext == "tsv") {
      compare_tsv(old_path, new_path)
    } else if (ext == "rds") {
      compare_rds(old_path, new_path)
    } else {
      # Figures (pdf/png/tiff etc.) - compare raw bytes only.
      if (identical(readBin(old_path, "raw", file.info(old_path)$size),
                     readBin(new_path, "raw", file.info(new_path)$size))) {
        list(status = "IDENTICAL", detail = NULL)
      } else {
        list(status = "DIFFERS (binary)", detail = "byte content differs (expected for regenerated plots; inspect visually)")
      }
    },
    error = function(e) list(status = "ERROR", detail = conditionMessage(e))
  )
  tibble(file = fname, status = res$status, detail = res$detail %||% NA_character_)
})

cat(sprintf("\nComparing %s/ (new) vs %s/ (baseline '%s')\n\n", current_dir, baseline_dir, tag))
for (i in seq_len(nrow(results))) {
  cat(sprintf("[%s] %s\n", results$status[i], results$file[i]))
  if (!is.na(results$detail[i])) cat(results$detail[i], "\n")
}

n_bad <- sum(!results$status %in% c("IDENTICAL", "CLOSE (bootstrap cols only)"))
cat(sprintf(
  "\n%d file(s) compared: %d identical, %d close-only-in-bootstrap-cols, %d needing review.\n",
  nrow(results),
  sum(results$status == "IDENTICAL"),
  sum(results$status == "CLOSE (bootstrap cols only)"),
  n_bad
))
if (n_bad > 0) quit(status = 1)
