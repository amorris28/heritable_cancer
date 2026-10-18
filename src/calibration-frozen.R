library(tidyverse)
library(readxl)
library(survival)

# Usage: calibration-frozen.R <phenotypes.tsv> <output_prefix>
# External-validation calibration of absolute risk. Nothing is fitted in the cohort:
#   HR_i  = exp(phs_raw_i - c)   published PHS601 weights; c is chosen so the cohort's
#                                mean HR is 1 (registry rates describe the average man)
#   H0(a) = integral of the Norwegian registry rate (Huynh-Le 2020) from age 45 to a
#   risk_i(a) = 1 - exp(-H0(a) * HR_i)
# compared with the Kaplan-Meier observed risk by age a, by PHS group, both treating
# everyone as at risk from age 40 and (if entry age exists) with delayed entry.
# events_by_age / n are those of the cohort; n_at_risk is the KM risk set at age a.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: calibration-frozen.R <phenotypes.tsv> <output_prefix>")
metadata_file <- args[1]
output_prefix <- args[2]

screen_ages <- c(50, 55, 60, 65)
n_groups <- 5
start_age <- 45 # registry's youngest bin is 45-49

dat <- read_tsv(metadata_file, show_col_types = FALSE)
if (!"phs_raw" %in% names(dat)) {
  stop("phs_raw not found in ", metadata_file, "; rerun scripts/02_preprocess_metadata.sh")
}
dat <- dat %>% filter(age >= 40, !is.na(phs_raw))

# Entry age (delayed entry) is available when first-follow-up dates are in the file
if (!"entry_age" %in% names(dat) && all(c("fdato", "date_first_follow_up") %in% names(dat))) {
  dat <- dat %>%
    mutate(entry_age = lubridate::interval(as.Date(fdato), as.Date(date_first_follow_up)) /
             lubridate::years(1))
}
has_entry <- "entry_age" %in% names(dat) && any(!is.na(dat$entry_age))
if (has_entry) {
  n_before <- nrow(dat)
  # Men without a valid entry age, or diagnosed/censored before entry, cannot be in the risk set
  dat <- dat %>% filter(!is.na(entry_age), entry_age < age)
  message(sprintf("Delayed entry: kept %d of %d men (median entry age %.1f)",
                  nrow(dat), n_before, median(dat$entry_age)))
} else {
  message("No entry age available; only the from-age-40 analysis is run")
}

# Registry baseline (same construction as fit-models-plot-figures-incidence-rate.R)
registry <- read_excel("data/external/ASIRs_per100k_Huynh-Le_2020.xlsx", skip = 2,
  col_names = c("age", "low_risk", "fav_interm", "unfav_interm", "high_risk",
                "regional", "metastatic")) %>%
  rowwise() %>%
  mutate(
    total = sum(across(low_risk:metastatic)),
    age_nice = if_else(age == "<50", "45-49", if_else(age == "≥85", "85-89", age)),
    median_age = median(seq(as.numeric(substr(age_nice, 1, 2)),
                            as.numeric(substr(age_nice, 4, 5)))),
    rate_per_year = total / 1e5
  ) %>%
  ungroup()
reg_fit <- lm(log(rate_per_year) ~ poly(median_age, 2, raw = TRUE), data = registry)
reg_hazard <- function(a) exp(predict(reg_fit, newdata = data.frame(median_age = a)))
H0 <- function(a) {
  vapply(a, function(x) integrate(reg_hazard, start_age, x)$value, numeric(1))
}

# Centre so mean(HR) = 1 in the cohort
centre <- log(mean(exp(dat$phs_raw)))
dat <- dat %>%
  mutate(hr = exp(phs_raw - centre), phs_group = ntile(phs_raw, n_groups))
message(sprintf("Centring constant c = %.3f; HR range %.2f-%.2f", centre,
                min(dat$hr), max(dat$hr)))

km_risk <- function(d, a, entry) {
  surv_obj <- if (entry) Surv(d$entry_age, d$age, d$status) else Surv(d$age, d$status)
  s <- summary(survfit(surv_obj ~ 1), times = a, extend = TRUE)
  c(observed = 1 - s$surv, obs_se = s$std.err, n_at_risk = s$n.risk)
}

analyses <- c("From age 40" = FALSE, "Delayed entry" = TRUE)
if (!has_entry) analyses <- analyses[1]

calib <- map_dfr(names(analyses), function(an) map_dfr(screen_ages, function(a) {
  entry <- analyses[[an]]
  d <- dat %>% mutate(pred = 1 - exp(-H0(a) * hr))
  by_group <- d %>%
    group_by(phs_group) %>%
    group_modify(function(g, key) {
      k <- km_risk(g, a, entry)
      tibble(n = nrow(g), events_by_age = sum(g$status == 1 & g$age <= a),
             predicted = mean(g$pred), observed = k[["observed"]], obs_se = k[["obs_se"]],
             n_at_risk = k[["n_at_risk"]])
    }) %>%
    ungroup() %>%
    mutate(group = as.character(phs_group))
  k_all <- km_risk(d, a, entry)
  overall <- tibble(group = "Overall", n = nrow(d),
                    events_by_age = sum(d$status == 1 & d$age <= a),
                    predicted = mean(d$pred), observed = k_all[["observed"]],
                    obs_se = k_all[["obs_se"]], n_at_risk = k_all[["n_at_risk"]])
  bind_rows(by_group, overall) %>%
    mutate(analysis = an, screen_age = a,
           oe = observed / predicted,
           oe_lo = pmax(observed - 1.96 * obs_se, 0) / predicted,
           oe_hi = (observed + 1.96 * obs_se) / predicted)
}))
write_tsv(calib %>% select(-phs_group), paste0(output_prefix, "_calibration_frozen.tsv"))
print(as.data.frame(calib %>% select(analysis, screen_age, group, n, events_by_age, n_at_risk, predicted, observed, oe, oe_lo, oe_hi)))

grp <- calib %>% filter(group != "Overall")
p <- ggplot(grp, aes(predicted, observed)) +
  geom_abline(linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(ymin = pmax(observed - 1.96 * obs_se, 0),
                    ymax = observed + 1.96 * obs_se), width = 0, colour = "grey60") +
  geom_point(aes(colour = analysis)) +
  facet_wrap(~screen_age, scales = "free", labeller = label_both) +
  labs(x = "Predicted cumulative incidence (published PHS601 weights + registry baseline)",
       y = "Observed cumulative incidence (Kaplan-Meier)") +
  theme_minimal()
ggsave(paste0(output_prefix, "_calibration_frozen.png"), p, width = 8, height = 7)
