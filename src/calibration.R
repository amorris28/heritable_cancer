library(tidyverse)
library(survival)

# Usage: calibration.R <phenotypes.tsv> <output_prefix>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: calibration.R <phenotypes.tsv> <output_prefix> [n_boot]")
metadata_file <- args[1]
output_prefix <- args[2]
n_boot <- if (length(args) >= 3) as.integer(args[3]) else 200L

screen_ages <- c(50, 55, 60, 65)
n_groups <- 10

dat <- read_tsv(metadata_file, show_col_types = FALSE) %>%
  filter(age >= 40, !is.na(phs))

run_calibration <- function(dat) {
  # Constant-HR model
  fit <- coxph(Surv(age, status) ~ phs, data = dat)
  base <- basehaz(fit, centered = FALSE)
  H0 <- function(t, b = base) approx(b$time, b$hazard, xout = t, method = "constant",
                                     yleft = 0, rule = 2)$y

  # Age-varying HR: split follow-up into periods with a separate PHS effect each
  cuts <- c(40, 60, 70, Inf)
  n_per <- length(cuts) - 1
  dat$id <- seq_len(nrow(dat))
  split_dat <- survSplit(Surv(age, status) ~ ., data = dat, cut = cuts[2:3],
                         start = "tstart") %>%
    mutate(period = findInterval(tstart, c(0, cuts[2:3])))
  for (k in seq_len(n_per)) split_dat[[paste0("phs", k)]] <- split_dat$phs * (split_dat$period == k)
  fit_av <- coxph(Surv(tstart, age, status) ~ phs1 + phs2 + phs3, data = split_dat)
  base_av <- basehaz(fit_av, centered = FALSE)
  beta_av <- coef(fit_av)
  beta_av[is.na(beta_av)] <- 0 # period without events: no PHS effect estimable

  cum_haz_av <- function(phs, a) {
    out <- 0
    for (k in seq_len(n_per)) {
      lo <- cuts[k]; hi <- min(cuts[k + 1], a)
      if (a > lo) out <- out + (H0(hi, base_av) - H0(lo, base_av)) * exp(beta_av[k] * phs)
    }
    out
  }

  dat <- dat %>% mutate(lp = predict(fit, type = "lp"),
                        phs_group = ntile(phs, n_groups))

  calib <- map_dfr(screen_ages, function(a) {
    dat %>%
      mutate(pred_const = 1 - exp(-H0(a) * exp(lp)),
             pred_av = 1 - exp(-cum_haz_av(phs, a))) %>%
      pivot_longer(c(pred_const, pred_av), names_to = "model", values_to = "pred",
                   names_prefix = "pred_") %>%
      mutate(model = recode(model, const = "Constant HR", av = "Age-varying HR")) %>%
      group_by(model, phs_group) %>%
      group_modify(function(d, key) {
        km <- survfit(Surv(age, status) ~ 1, data = d)
        s <- summary(km, times = a, extend = TRUE)
        tibble(n = nrow(d), predicted = mean(d$pred),
               observed = 1 - s$surv, obs_se = s$std.err)
      }) %>%
      ungroup() %>%
      mutate(screen_age = a)
  })
    list(calib = calib, fit_av = fit_av)
}

res <- run_calibration(dat)
calib <- res$calib
print(summary(res$fit_av)$conf.int)

write_tsv(calib, paste0(output_prefix, "_calibration_by_age.tsv"))

# Calibration slope/intercept: regress observed on predicted across PHS groups
summarise_calib <- function(calib) calib %>%
  filter(!is.na(observed)) %>%
  group_by(model, screen_age) %>%
  summarise(
    intercept = coef(lm(observed ~ predicted, weights = n))[1],
    slope = coef(lm(observed ~ predicted, weights = n))[2],
    .groups = "drop"
  )

summ <- summarise_calib(calib)

# Bootstrap (resample men) to get percentile CIs for slope and intercept
set.seed(42)
boot <- map_dfr(seq_len(n_boot), function(i) {
  d <- dat[sample.int(nrow(dat), replace = TRUE), ]
  summarise_calib(run_calibration(d)$calib) %>% mutate(iter = i)
})
ci <- boot %>%
  group_by(model, screen_age) %>%
  summarise(across(c(intercept, slope),
                   list(lo = ~quantile(.x, 0.025, na.rm = TRUE),
                        hi = ~quantile(.x, 0.975, na.rm = TRUE))),
            .groups = "drop")
summ <- left_join(summ, ci, by = c("model", "screen_age"))
write_tsv(summ, paste0(output_prefix, "_calibration_slope_intercept.tsv"))
print(summ)

labels <- summ %>%
  mutate(label = sprintf("intercept = %.3f\nslope = %.2f (%.2f-%.2f)", intercept, slope,
                         slope_lo, slope_hi),
         y = ifelse(model == "Constant HR", Inf, Inf),
         vj = ifelse(model == "Constant HR", 1.2, 3.4))

p <- ggplot(calib, aes(predicted, observed)) +
  geom_abline(linetype = "dashed", colour = "grey50") +
  geom_smooth(aes(weight = n, colour = model), method = "lm", se = FALSE,
              na.rm = TRUE,
              linewidth = 0.7) +
  geom_errorbar(aes(ymin = pmax(observed - 1.96 * obs_se, 0),
                    ymax = observed + 1.96 * obs_se), width = 0, colour = "grey60") +
  geom_point(aes(colour = model)) +
  geom_text(data = labels, aes(x = -Inf, y = Inf, label = label, colour = model),
            hjust = -0.1, vjust = labels$vj, size = 3, inherit.aes = FALSE) +
  facet_wrap(~screen_age, scales = "free", labeller = label_both) +
  labs(x = "Predicted cumulative incidence",
       y = "Observed cumulative incidence (Kaplan-Meier)") +
  theme_minimal()

ggsave(paste0(output_prefix, "_calibration.png"), p, width = 8, height = 7)
