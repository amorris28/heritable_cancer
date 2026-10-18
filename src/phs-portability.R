library(tidyverse)
library(survival)

# Usage: phs-portability.R <phenotypes.tsv> <output_prefix>
# Fits Cox models on the unstandardized PHS (phs_raw, plink weighted sum) to ask
# whether the published PHS601 weights transport: if the Betas are log-HRs, a
# coefficient of 1 means the score's effect size carries over unchanged.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: phs-portability.R <phenotypes.tsv> <output_prefix>")
metadata_file <- args[1]
output_prefix <- args[2]

dat <- read_tsv(metadata_file, show_col_types = FALSE)
if (!"phs_raw" %in% names(dat)) {
  stop("phs_raw not found in ", metadata_file, "; rerun scripts/02_preprocess_metadata.sh")
}
dat <- dat %>% filter(age >= 40, !is.na(phs_raw))

endpoints <- list(
  any_pca    = c(time = "age", event = "status"),
  metastatic = c(time = "age_of_metastasis", event = "metastasis"),
  pca_death  = c(time = "age_at_pca_death", event = "pca_death")
)

fit_endpoint <- function(spec, name) {
  d <- dat %>% filter(!is.na(.data[[spec[["time"]]]]), !is.na(.data[[spec[["event"]]]]))
  form <- as.formula(sprintf("Surv(%s, %s) ~ phs_raw", spec[["time"]], spec[["event"]]))
  fit <- coxph(form, data = d)
  est <- summary(fit)$coefficients
  beta <- est[1, "coef"]
  se <- est[1, "se(coef)"]
  tibble(
    model = name,
    n = nrow(d),
    events = sum(d[[spec[["event"]]]]),
    phs_raw_mean = mean(d$phs_raw),
    phs_raw_sd = sd(d$phs_raw),
    slope = beta,
    slope_low = beta - 1.96 * se,
    slope_high = beta + 1.96 * se,
    p_slope_ne_1 = 2 * pnorm(-abs((beta - 1) / se)),
    hr_per_sd = exp(beta * sd(d$phs_raw)),
    hr_per_sd_low = exp((beta - 1.96 * se) * sd(d$phs_raw)),
    hr_per_sd_high = exp((beta + 1.96 * se) * sd(d$phs_raw)),
    # Expected HR per SD if the published weights transported (slope = 1)
    hr_per_sd_if_slope_1 = exp(sd(d$phs_raw))
  )
}

results <- imap_dfr(endpoints, fit_endpoint) %>% mutate(cohort = output_prefix)
print(as.data.frame(results))
write_tsv(results, paste0(output_prefix, "_phs_portability.tsv"))
