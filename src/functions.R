# Drop exact duplicate rows, then drop every row whose value in any of id_cols
# is shared with another row (NA ids are ignored). Conflicting duplicates are
# removed entirely because the correct record cannot be determined.
resolve_duplicate_ids <- function(data, id_cols, label = "data") {
  n_start <- nrow(data)
  data <- dplyr::distinct(data)
  message(label, ": removed ", n_start - nrow(data), " exact duplicate rows")

  dup <- rep(FALSE, nrow(data))
  for (col in id_cols) {
    ids <- data[[col]]
    d <- !is.na(ids) & (duplicated(ids) | duplicated(ids, fromLast = TRUE))
    message(label, ": ", length(unique(ids[d])), " conflicting duplicated ", col, " values")
    dup <- dup | d
  }
  if (any(dup)) {
    n_ids <- nrow(dplyr::distinct(data[dup, id_cols]))
    message(
      label, ": removing ", sum(dup), " rows (", n_ids,
      " unique ", paste(id_cols, collapse = "/"), " combinations) with conflicting duplicate ids"
    )
    print(data[dup, ], n = Inf)
  }
  data[!dup, ]
}

# Function definitons
##  Data processing
import_phs_data <- function(model_file, inverse = FALSE) {
  # phs_raw is the unstandardized plink weighted sum (log-hazard units if the
  # model Betas are log-HRs); phs is the cohort-specific z-score used for the
  # primary analyses.
  readr::read_table(model_file) %>%
    dplyr::select(bfile_id = IID, phs_raw = SCORESUM) %>%
    dplyr::mutate(phs = as.vector(scale(phs_raw, center = TRUE, scale = TRUE))) ->
  phs_data
  if (inverse) {
    dplyr::mutate(phs_data, phs_raw = phs_raw * -1, phs = phs * -1) -> phs_data
  }
  return(phs_data)
}

process_combined_data <- function(metadata, phs_data) {
  combined_data <- metadata %>%
    left_join(phs_data, by = "bfile_id") %>%
    drop_na(phs, age, status)

  return(combined_data)
}

##  Formatting
pvalr <- function(pvals, sig.limit = .001, digits = 3, html = FALSE) {
  roundr <- function(x, digits = 1) {
    res <- sprintf(paste0("%.", digits, "f"), x)
    zzz <- paste0("0.", paste(rep("0", digits), collapse = ""))
    res[res == paste0("-", zzz)] <- zzz
    res
  }

  sapply(pvals, function(x, sig.limit) {
    if (x < sig.limit) {
      if (html) {
        return(sprintf("&lt; %s", format(sig.limit)))
      } else {
        return(sprintf("< %s", format(sig.limit)))
      }
    }
    if (x > .1) {
      return(roundr(x, digits = 2))
    } else {
      return(roundr(x, digits = digits))
    }
  }, sig.limit = sig.limit)
}

##  PHS percentile-band / calibration helpers (replaces hazrd:: dependency)

# Cumulative-incidence plot for PHS percentile strata (replaces hazrd::phs_abs_risk plot output)
plot_cuminc_strata <- function(cuminc_df) {
  ggplot2::ggplot(cuminc_df, ggplot2::aes(x = time, y = risk_estimate, color = stratum)) +
    ggplot2::geom_step(linewidth = 0.8) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = conf_low, ymax = conf_high, fill = stratum), alpha = 0.15, color = NA) +
    ggplot2::scale_y_continuous(labels = scales::percent_format()) +
    ggplot2::scale_color_brewer(palette = "Set1") +
    ggplot2::scale_fill_brewer(palette = "Set1") +
    ggplot2::labs(x = "Time", y = "Cumulative Risk", color = "Stratum", fill = "Stratum") +
    ggplot2::theme_minimal()
}

# Numbers-at-risk strip below a cumulative-incidence plot (replaces hazrd::phs_risk_table plot output)
plot_risk_table <- function(risk_tbl, x_range = NULL) {
  strata_levels <- levels(risk_tbl$stratum)
  p <- ggplot2::ggplot(risk_tbl, ggplot2::aes(x = time, y = stratum, label = n_risk, color = stratum)) +
    ggplot2::geom_text(size = 3.2, show.legend = FALSE) +
    ggplot2::scale_y_discrete(limits = rev(strata_levels)) +
    ggplot2::scale_color_brewer(palette = "Set1", drop = FALSE) +
    ggplot2::labs(x = "Time", y = NULL, title = "Number at risk") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid   = ggplot2::element_blank(),
      plot.title   = ggplot2::element_text(size = 9, face = "plain"),
      plot.margin  = ggplot2::margin(0, 5.5, 5.5, 5.5)
    )
  if (!is.null(x_range)) p <- p + ggplot2::scale_x_continuous(limits = x_range)
  p
}

# Percentile-band membership, matching hazrd's ecdf-based stratum definition:
# lo/hi are percentiles in [0, 1]; upper bound is inclusive only for hi >= 1.
phs_band_mask <- function(phs_vals, lo, hi, ref_phs = phs_vals) {
  phs_pct <- as.numeric(ecdf(ref_phs)(phs_vals))
  upper_mask <- if (hi >= 1) phs_pct <= hi else phs_pct < hi
  phs_pct >= lo & upper_mask
}

# Cumulative incidence (1 - KM) per PHS percentile band, tidy long format.
phs_cuminc_strata <- function(data, phs = "phs", time = "age", event = "status",
                               intervals = list(c(0.95, 1), c(0.80, 1), c(0.30, 0.70), c(0, 0.20))) {
  labels <- sapply(intervals, function(iv) paste0(round(iv[1] * 100), "-", round(iv[2] * 100), "%"))
  purrr::map2_dfr(intervals, labels, function(iv, lbl) {
    sub <- data[phs_band_mask(data[[phs]], iv[1], iv[2]), , drop = FALSE]
    fit <- survival::survfit(survival::Surv(sub[[time]], sub[[event]]) ~ 1)
    tibble::tibble(
      time          = fit$time,
      risk_estimate = pmin(pmax(1 - fit$surv,  0), 1),
      conf_low      = pmin(pmax(1 - fit$upper, 0), 1),
      conf_high     = pmin(pmax(1 - fit$lower, 0), 1),
      n_risk        = fit$n.risk,
      n_event       = fit$n.event,
      stratum       = lbl
    )
  }) %>%
    dplyr::mutate(stratum = factor(stratum, levels = labels))
}

# Numbers-at-risk per stratum at requested time points (step-function lookback).
phs_risk_at_times <- function(cuminc_df, time_points) {
  strata_levels <- levels(cuminc_df$stratum)
  purrr::map_dfr(strata_levels, function(s) {
    sub <- dplyr::filter(cuminc_df, stratum == s)
    n_at_t <- sapply(time_points, function(t) {
      prior <- dplyr::filter(sub, time <= t)
      if (nrow(prior) == 0) max(sub$n_risk, na.rm = TRUE) else dplyr::last(prior$n_risk)
    })
    tibble::tibble(
      time    = time_points,
      stratum = factor(s, levels = strata_levels),
      n_risk  = as.integer(round(n_at_t))
    )
  })
}

# Cox-predicted survival curve for an individual at given PHS percentile(s).
phs_cox_curve_at_percentile <- function(cxph, ref_phs, percentiles) {
  phs_at_pct <- stats::quantile(ref_phs, probs = percentiles, type = 7)
  purrr::map2_dfr(percentiles, phs_at_pct, function(pct, phs_val) {
    fit <- survival::survfit(cxph, newdata = data.frame(phs = as.numeric(phs_val)))
    tibble::tibble(
      time              = fit$time,
      estimate          = pmin(pmax(fit$surv, 0), 1),
      percentile_value  = pct
    )
  })
}

# Group HR from the mean linear-predictor difference between two PHS percentile
# bands (hazrd's "continuous_group" method).
calc_hr_band <- function(cxph, data, phs, numerator, denominator) {
  beta <- as.numeric(cxph$coefficients["phs"])
  ix_num <- which(phs_band_mask(data[[phs]], numerator[1], numerator[2]))
  ix_den <- which(phs_band_mask(data[[phs]], denominator[1], denominator[2]))
  beta_phs <- beta * data[[phs]]
  tibble::tibble(
    metric = paste0(
      "HR[", round(numerator[1] * 100), "-", round(numerator[2] * 100),
      "]_[", round(denominator[1] * 100), "-", round(denominator[2] * 100), "]"
    ),
    estimate      = exp(mean(beta_phs[ix_num]) - mean(beta_phs[ix_den])),
    n_numerator   = length(ix_num),
    n_denominator = length(ix_den)
  )
}

# HR per 1-SD increase in PHS.
calc_hr_sd <- function(cxph, data, phs) {
  beta <- as.numeric(cxph$coefficients["phs"])
  exp(beta * sd(data[[phs]]))
}

# Generic bootstrap CI helper: stat_fn(data) -> named numeric vector of estimates.
boot_ci <- function(data, stat_fn, n_boot = 999, seed = NULL, ci_level = 0.95) {
  if (!is.null(seed)) set.seed(seed)
  boot_out <- boot::boot(
    data      = data,
    statistic = function(d, indices) stat_fn(d[indices, , drop = FALSE]),
    R         = n_boot
  )
  metric_names <- names(boot_out$t0)
  purrr::map_dfr(seq_along(metric_names), function(i) {
    ci <- tryCatch(
      boot::boot.ci(boot_out, index = i, type = "perc", conf = ci_level),
      error = function(e) NULL
    )
    tibble::tibble(
      metric    = metric_names[i],
      estimate  = boot_out$t0[i],
      conf_low  = if (!is.null(ci)) ci$percent[4] else NA_real_,
      conf_high = if (!is.null(ci)) ci$percent[5] else NA_real_,
      se        = sd(boot_out$t[, i], na.rm = TRUE)
    )
  })
}

# HR[band pairs] + C_index + HR_SD with bootstrap CIs, in the same tibble shape
# hazrd::phs_metrics() used to return (metric, estimate, conf_low, conf_high,
# se, n_numerator, n_denominator, method, adjusted).
phs_performance_metrics <- function(data, time = "age", event = "status", phs = "phs",
                                     hr_pairs = list(
                                       list(numerator = c(0.80, 1.00), denominator = c(0.00, 0.20)),
                                       list(numerator = c(0.80, 1.00), denominator = c(0.30, 0.70)),
                                       list(numerator = c(0.95, 1.00), denominator = c(0.30, 0.70))
                                     ),
                                     n_boot = 999, seed = 42) {
  form <- stats::as.formula(paste0("survival::Surv(", time, ", ", event, ") ~ ", phs))

  point_stat <- function(d) {
    cxph <- survival::coxph(form, data = d)
    hr_tbl <- purrr::map_dfr(hr_pairs, ~ calc_hr_band(cxph, d, phs, .x$numerator, .x$denominator))
    stats::setNames(
      c(hr_tbl$estimate, cxph$concordance[["concordance"]], calc_hr_sd(cxph, d, phs)),
      c(hr_tbl$metric, "C_index", "HR_SD")
    )
  }

  boot_est <- boot_ci(data, point_stat, n_boot = n_boot, seed = seed)

  cxph <- survival::coxph(form, data = data)
  n_tbl <- purrr::map_dfr(hr_pairs, ~ calc_hr_band(cxph, data, phs, .x$numerator, .x$denominator)) %>%
    dplyr::select(metric, n_numerator, n_denominator)

  boot_est %>%
    dplyr::left_join(n_tbl, by = "metric") %>%
    dplyr::mutate(
      method = dplyr::case_when(
        stringr::str_starts(metric, "HR\\[") ~ "continuous_group",
        metric == "C_index"                  ~ "harrell",
        TRUE                                  ~ NA_character_
      ),
      adjusted = FALSE
    )
}


