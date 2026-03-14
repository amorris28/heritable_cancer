library(tidyverse)
library(stringr)
library(hazrd)
library(survival)
library(patchwork)
source("src/functions.R")

###########################################
# Output flags

#########################################
# First argument: PHS file
# tab delimited file, no header, first column ID, second column PHS score
# Second argument: metadata file
# column names. Missing is okay, different order is okay, must match string exactly:
# pid	death	wgs	age_lastobs	age_diagnosis	status	age	bfile_id	age_of_metastasis
# metastasis_code	metastasis	conversion_to_metastasis_years	doed_av_pca	age_at_death
# diagnosis_to_death_years	brca2_carrier	hoxb13_carrier	hoxb13_snp
# tab delimited file with column names bfile_id, age, status.
# Additional columns will be ignored.
# Third argument:
# boolean indicating whether to inverse (x * -1) the PHS scores
# to reverse the direction of effect.

args <- commandArgs(trailingOnly = TRUE)

metadata_file <- args[1] # Phenotype data
perform_bool <- args[2]
output_prefix <- args[3]
bootstrap_iters <- args[4]

#########################################
# Import PHS data and metadata
combined_data <- read_tsv(metadata_file)

carriers <-
  combined_data %>%
  select(
    status, age, phs, age_of_metastasis, metastasis, wgs, ends_with("_carrier")
  ) %>%
  mutate(across(ends_with("_carrier"), ~ if_else(is.na(.x), 0, .x)))

wgs <-
  carriers %>%
  filter(wgs == 1) %>%
  # any_carrier is 1 for people who have any of the listed RPVs
  mutate(any_carrier = case_when(if_any(ends_with("_carrier"), ~ .x == 1) ~ 1, .default = 0))

wgs_summary <-
  wgs %>%
  select(status, wgs, ends_with("_carrier")) %>%
  summarize(across(everything(), list(sum)))


# Fisher test
variants <- c("BRCA1", "BRCA2", "ATM", "PALB2", "CHEK2", "HOXB13", "MLH1", "MSH2", "MSH6", "PMS2", "TP53", "EPCAM", "any") # , "all")
fisher_results <- list()
for (i in 1:length(variants)) {
  variant <- paste0(variants[i], "_carrier")
  print(variant)
  fisher_results[[variant]] <-
    wgs %>%
    select(status, all_of(variant)) %>%
    mutate(case_factor = factor(status, levels = c(1, 0), labels = c("Cases", "Controls"))) %>%
    mutate(carrier_factor = factor(.data[[variant]], levels = c(1, 0), labels = c("Carrier", "Not_Carrier"))) %>%
    select(case_factor, carrier_factor) %>%
    table() %>%
    fisher.test(alternative = "g")
}
fisher_tibble <- NULL
for (i in seq_along(fisher_results)) {
  temp <- broom::tidy(fisher_results[[i]])
  temp$gene <- names(fisher_results)[i]
  fisher_tibble <- rbind(fisher_tibble, temp)
}
fisher_tibble <-
  fisher_tibble %>%
  mutate(p.adjust = p.adjust(p.value, method = "fdr"))


# survival analysis
## PHS
phs_wgs_fit <- coxph(Surv(age, status) ~ phs, wgs)
summary(phs_wgs_fit)
phs_wgs_coxph <- broom::tidy(phs_wgs_fit, exponentiate = TRUE, conf.int = TRUE)
cox.zph(phs_wgs_fit)
phs_wgs_cindex <-
  phs_wgs_fit %>%
  broom::glance() %>%
  select(concordance, std.error.concordance)

## HOXB13
hoxb13_wgs_fit <- coxph(Surv(age, status) ~ HOXB13_carrier, wgs)
summary(hoxb13_wgs_fit)
hoxb13_wgs_coxph <- broom::tidy(hoxb13_wgs_fit,
  exponentiate = TRUE, conf.int = TRUE
)
cox.zph(hoxb13_wgs_fit)
hoxb13_wgs_cindex <-
  hoxb13_wgs_fit %>%
  broom::glance() %>%
  select(concordance, std.error.concordance)

## PHS + HOXB13
phs_hoxb13_wgs_fit <- coxph(Surv(age, status) ~ phs + HOXB13_carrier, wgs)
summary(phs_hoxb13_wgs_fit)
phs_hoxb13_wgs_coxph <- broom::tidy(phs_hoxb13_wgs_fit,
  exponentiate = TRUE, conf.int = TRUE
)
cox.zph(phs_hoxb13_wgs_fit)
phs_hoxb13_wgs_cindex <-
  phs_hoxb13_wgs_fit %>%
  broom::glance() %>%
  select(concordance, std.error.concordance)

phs_hoxb13_wgs_inter <- coxph(Surv(age, status) ~ phs * HOXB13_carrier, wgs)
summary(phs_hoxb13_wgs_inter)
phs_hoxb13_inter_coxph <- broom::tidy(phs_hoxb13_wgs_inter,
  exponentiate = TRUE, conf.int = TRUE
)
cox.zph(phs_hoxb13_wgs_inter)
phs_hoxb13_inter_cindex <-
  phs_hoxb13_wgs_inter %>%
  broom::glance() %>%
  select(concordance, std.error.concordance)

# HOXB13 frequency in each dataset
sum(wgs$HOXB13_carrier) / nrow(wgs) # 1.7%
sum(combined_data$hoxb13_snp) / nrow(combined_data) # 0.5%

# Performance metrics
perform <- phs_metrics(
  wgs,
  metrics  = c("HR", "C_index"),
  hr_pairs = list(
    list(numerator = c(100 - 0.017, 1.00), denominator = c(0.30, 0.70))
  ),
  bootstrap = TRUE,
  n_boot    = 999,
  seed      = 42
)
perform <- dplyr::mutate(perform, model = output_prefix)
write_tsv(perform, paste(output_prefix, "wgs_phs_performance_95.tsv",
  sep = "_"
))

# Cumulative incidence curve for PHS
wgs_inc <- phs_abs_risk(
  wgs,
  intervals = list(c(0.95, 1), c(0.80, 1), c(0.30, 0.70), c(0, 0.20))
)

# KM Curves for HOXB13
km_curves <- data.frame()
noncarrier_data <- wgs[wgs$HOXB13_carrier == 0, ]
carrier_data <- wgs[wgs$HOXB13_carrier == 1, ]
mod <- survfit(Surv(age, status) ~ 1, data = noncarrier_data)
mod_data <- data.frame(
  time = mod$time, n.risk = mod$n.risk,
  n.event = mod$n.event, n.censor = mod$n.censor, estimate = mod$surv,
  std.error = mod$std.err, conf.high = mod$upper, conf.low = mod$lower,
  cumhaz = mod$cumhaz
)
curven <- select(
  mod_data, time, estimate, conf.low,
  conf.high, cumhaz
)
curven$label <- "noncarrier"
km_curves <- rbind(km_curves, curven)
mod <- survfit(Surv(age, status) ~ 1, data = carrier_data)
mod_data <- data.frame(
  time = mod$time, n.risk = mod$n.risk,
  n.event = mod$n.event, n.censor = mod$n.censor, estimate = mod$surv,
  std.error = mod$std.err, conf.high = mod$upper, conf.low = mod$lower,
  cumhaz = mod$cumhaz
)
curven <- select(
  mod_data, time, estimate, conf.low,
  conf.high, cumhaz
)
curven$label <- "carrier"
km_curves <- rbind(km_curves, curven)
yaxis <- "Disease-free Survival"
ylims <- c(0, 1)

km_curves %>%
  mutate(
    estimate = (1 - estimate) * 100,
    conf.high = (1 - conf.high) * 100,
    conf.low = (1 - conf.low) * 100
  ) ->
km_curves
yaxis <- "Cumulative incidence, %\n(Absolute risk)"
ylims <- c(0, 100)

color_set <- RColorBrewer::brewer.pal(9, "Set1")
hoxb13_p <- ggplot(km_curves, aes(
  x = time,
  y = estimate,
  ymin = conf.low,
  ymax = conf.high,
  col = label,
  fill = label
)) +
  geom_ribbon(
    alpha = 0.1,
    color = 0
  ) +
  geom_step() +
  theme_classic() +
  lims(x = c(40, 100), y = ylims) +
  labs(x = "Age", y = yaxis) +
  scale_color_manual(
    values = c(color_set[5], color_set[2]),
    name = "HOXB13 Carrier"
  ) +
  scale_fill_manual(
    values = c(color_set[5], color_set[2]),
    name = "HOXB13 Carrier"
  ) +
  theme(aspect.ratio = 1) +
  theme(legend.position = "bottom")

wgs_p <-
  wgs_inc +
    hoxb13_p +
    plot_annotation(tag_levels = "A") &
    theme(legend.position = "bottom")

saveRDS(wgs_p, paste(output_prefix, "wgs_km_curve.rds", sep = "_"))


#################################################################################
# Does having a high PHS icnrease your risk of carrying a rare variant?
#################################################################################
binomial_phs_hoxb13 <- glm(HOXB13_carrier ~ phs, family = "binomial", data = wgs)
summary(binomial_phs_hoxb13)
broom::tidy(binomial_phs_hoxb13, conf.int = TRUE, exponentiate = TRUE) %>%
  write_tsv(paste(output_prefix, "binomial_phs_hoxb13.tsv", sep = "_"))

phs_wgs_coxph$model <- "phs_wgs"
hoxb13_wgs_coxph$model <- "hoxb13_wgs"
phs_hoxb13_wgs_coxph$model <- "phs_hoxb13"
phs_hoxb13_inter_coxph$model <- "phs_hoxb13_inter"

wgs_results <-
  rbind(
    phs_wgs_coxph,
    hoxb13_wgs_coxph,
    phs_hoxb13_wgs_coxph,
    phs_hoxb13_inter_coxph
  )
write_tsv(wgs_results, paste(output_prefix, "wgs.tsv", sep = "_"))

phs_wgs_cindex$model <- "phs_wgs"
hoxb13_wgs_cindex$model <- "hoxb13_wgs"
phs_hoxb13_wgs_cindex$model <- "phs_hoxb13"
phs_hoxb13_inter_cindex$model <- "phs_hoxb13_inter"

cindex_results <-
  rbind(
    phs_wgs_cindex,
    hoxb13_wgs_cindex,
    phs_hoxb13_wgs_cindex,
    phs_hoxb13_inter_cindex
  )

write_tsv(cindex_results, paste(output_prefix, "wgs_cindex.tsv", sep = "_"))
write_tsv(fisher_tibble, paste(output_prefix, "fisher.tsv", sep = "_"))
write_tsv(wgs_summary, paste(output_prefix, "wgs_summary.tsv", sep = "_"))
