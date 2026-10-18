library(tidyverse)
library(survival)

# Usage: calibration.R <phenotypes.tsv> <output_prefix>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: calibration.R <phenotypes.tsv> <output_prefix>")
metadata_file <- args[1]
output_prefix <- args[2]

screen_ages <- c(50, 55, 60, 65)
n_groups <- 10

dat <- read_tsv(metadata_file, show_col_types = FALSE) %>%
  filter(age >= 40, !is.na(phs))

fit <- coxph(Surv(age, status) ~ phs, data = dat)
base <- basehaz(fit, centered = FALSE)
H0 <- function(t) approx(base$time, base$hazard, xout = t, method = "constant",
                         yleft = 0, rule = 2)$y

# Per-person predicted risk of diagnosis by each screening age (Cox, apparent)
dat <- dat %>% mutate(lp = predict(fit, type = "lp"),
                      phs_group = ntile(phs, n_groups))

calib <- map_dfr(screen_ages, function(a) {
  dat %>%
    mutate(pred = 1 - exp(-H0(a) * exp(lp))) %>%
    group_by(phs_group) %>%
    group_modify(function(d, key) {
      km <- survfit(Surv(age, status) ~ 1, data = d)
      s <- summary(km, times = a, extend = TRUE)
      tibble(n = nrow(d), predicted = mean(d$pred),
             observed = 1 - s$surv, obs_se = s$std.err)
    }) %>%
    ungroup() %>%
    mutate(screen_age = a)
})

write_tsv(calib, paste0(output_prefix, "_calibration_by_age.tsv"))

# Calibration slope/intercept: regress observed on predicted across PHS groups
summ <- calib %>%
  group_by(screen_age) %>%
  summarise(
    intercept = coef(lm(observed ~ predicted, weights = n))[1],
    slope = coef(lm(observed ~ predicted, weights = n))[2],
    .groups = "drop"
  )
write_tsv(summ, paste0(output_prefix, "_calibration_slope_intercept.tsv"))
print(summ)

labels <- summ %>%
  mutate(label = sprintf("intercept = %.3f\nslope = %.2f", intercept, slope))

p <- ggplot(calib, aes(predicted, observed)) +
  geom_abline(linetype = "dashed", colour = "grey50") +
  geom_smooth(aes(weight = n), method = "lm", se = FALSE, colour = "red",
              linewidth = 0.7) +
  geom_errorbar(aes(ymin = pmax(observed - 1.96 * obs_se, 0),
                    ymax = observed + 1.96 * obs_se), width = 0, colour = "grey60") +
  geom_point() +
  geom_text(data = labels, aes(x = -Inf, y = Inf, label = label),
            hjust = -0.1, vjust = 1.2, size = 3, inherit.aes = FALSE) +
  facet_wrap(~screen_age, scales = "free", labeller = label_both) +
  labs(x = "Predicted cumulative incidence",
       y = "Observed cumulative incidence (Kaplan-Meier)") +
  theme_minimal()

ggsave(paste0(output_prefix, "_calibration.png"), p, width = 8, height = 7)
