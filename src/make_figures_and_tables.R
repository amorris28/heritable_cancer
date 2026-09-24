library(tidyverse)
library(patchwork)
library(hazrd)

source("src/functions.R")

# Import results
p2521_perform <- readr::read_tsv("data/public/p2521_any_pca_performance.tsv")
p212_perform <- readr::read_tsv("data/public/p212_any_pca_performance.tsv")
p2521_meta <- readr::read_tsv("data/public/p2521_metastatic_performance.tsv")
p212_meta <- readr::read_tsv("data/public/p212_metastatic_performance.tsv")

p212_zph <- readr::read_tsv("data/public/p212_phs_zph_tidy.tsv")
p2521_zph <- readr::read_tsv("data/public/p2521_phs_zph_tidy.tsv")
p2521_sample_size <- readr::read_tsv("data/public/p2521_sample_size.tsv")
p212_sample_size <- readr::read_tsv("data/public/p212_sample_size.tsv")

p212_coxph_results <- readr::read_tsv("data/public/p212_coxph.tsv")
p2521_coxph_results <- readr::read_tsv("data/public/p2521_coxph.tsv")
wgs_summary <- readr::read_tsv("data/public/p212_wgs_summary.tsv")
fisher <- readr::read_tsv("data/public/p212_fisher.tsv")

# wgs_phs <- readr::read_tsv("data/public/p212_wgs_phs_performance.tsv")
# wgs_phs_95 <- readr::read_tsv("data/public/p212_wgs_phs_performance_95.tsv")
binomial_phs_hoxb13 <- readr::read_tsv("data/public/p212_binomial_phs_hoxb13.tsv")
wgs_results <- readr::read_tsv("data/public/p212_wgs.tsv")

wgs_cindex <- readr::read_tsv("data/public/p212_wgs_cindex.tsv")

###############################################################
# PHS K-M Curve Figure for Any-PCa
###############################################################

p212_ar <- base::readRDS("data/public/p212_any_pca_abs_risk.rds")
p212_rt <- base::readRDS("data/public/p212_any_pca_risk_table.rds")
p2521_ar <- base::readRDS("data/public/p2521_any_pca_abs_risk.rds")
p2521_rt <- base::readRDS("data/public/p2521_any_pca_risk_table.rds")

p212_ar <- p212_ar +
  labs(title = "NFCS - Any-stage PCa",
      y = "Cumulative Incidence (%)",
      x = "Age",
      color = "PHS percentile",
      fill = "PHS percentile") +
      lims(x = c(45, 90)) +
      theme_classic()
p212_rt <- p212_rt + labs(x = NULL) +
      lims(x = c(45, 90))
p2521_ar <- p2521_ar +
  labs(title = "HUSK - Any-stage PCa",
      y = "Cumulative Incidence (%)",
      x = "Age",
      color = "PHS percentile",
      fill = "PHS percentile") +
      lims(x = c(45, 90)) +
      theme_classic()
p2521_rt <- p2521_rt + labs(x = NULL) +
      lims(x = c(45, 90))

p212 <- hazrd::attach_risk_table(p212_ar, p212_rt)
p2521 <- hazrd::attach_risk_table(p2521_ar, p2521_rt)

any_pca <- (p2521 | p212) +
  # patchwork::plot_annotation(tag_levels = "A") +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("figures_and_tables/km_curve.pdf", any_pca,
  height = 8 * 0.618, width = 8
)
ggsave("figures_and_tables/km_curve.png", any_pca,
  height = 8 * 0.618, width = 8
)
ggsave("figures_and_tables/km_curve.tiff", any_pca,
  height = 8 * 0.618, width = 8, dpi = 300
)

###############################################################
# PHS Any-PCa and Meta-PCa Table
###############################################################
p2521_sample_size$cohort <- "HUSK"
p212_sample_size$cohort <- "NFCS"
p2521_perform$cohort <- "HUSK"
p2521_perform$PCa <- "Any"
p212_perform$cohort <- "NFCS"
p212_perform$PCa <- "Any"
p2521_meta$cohort <- "HUSK"
p2521_meta$PCa <- "Aggr"
p212_meta$cohort <- "NFCS"
p212_meta$PCa <- "Aggr"
p212_coxph_results$cohort <- "NFCS"
p2521_coxph_results$cohort <- "HUSK"

coxph_results <-
  p2521_coxph_results %>%
  rbind(p212_coxph_results) %>%
  mutate(PCa = if_else(model == "any_pca", "Any",
    if_else(model == "metastatic", "Aggr", NA)
  )) %>%
  filter(model %in% c("any_pca", "metastatic")) %>%
  select(estimate, conf_low = conf.low, conf_high = conf.high, PCa, cohort)

coxph_results$metric <- "HR_cox_pval"

sample_size <-
  p2521_sample_size %>%
  rbind(p212_sample_size) %>%
  mutate(PCa = if_else(group == "any_pca", "Any",
    if_else(group == "metastatic", "Aggr", NA)
  ))

table_coxph <-
  rbind(p2521_perform, p212_perform, p2521_meta, p212_meta) %>%
  mutate(metric = dplyr::recode(metric,
    "HR[80-100]_[0-20]"  = "HR80_20",
    "HR[80-100]_[30-70]" = "HR80_50",
    "HR[95-100]_[30-70]" = "HR95_50",
    "C_index"            = "C_Index",
    "HR_SD"              = "HR_SD"
  )) %>%
  filter(metric %in% c("HR_SD", "HR80_20", "HR80_50", "HR95_50", "C_Index")) %>%
  select(metric, estimate, conf_low, conf_high, PCa, cohort) %>%
  mutate(value = sprintf("%.2f (%.2f, %.2f)", estimate, conf_low, conf_high)) %>%
  select(metric, value, PCa, cohort) %>%
  pivot_wider(names_from = "metric", values_from = "value") %>%
  left_join(sample_size, by = join_by(cohort, PCa)) %>%
  select(PCa,
    Cohort = cohort, n = total, Cases = cases,
    HR_SD, HR80_20, HR80_50, HR95_50, C_Index
  ) %>%
  rename(
    `HR_SD`   = HR_SD,
    `HR80/20` = HR80_20,
    `HR80/50` = HR80_50,
    `HR95/50` = HR95_50,
    `C-Index` = C_Index
  ) %>%
  select(-n)

write_tsv(table_coxph, "figures_and_tables/table_coxph.tsv")

###############################################################
# PHS K-M Curve Figure for Meta-PCa
###############################################################
p212_ar <- readRDS("data/public/p212_metastatic_abs_risk.rds")
p212_rt <- readRDS("data/public/p212_metastatic_risk_table.rds")
p2521_ar <- readRDS("data/public/p2521_metastatic_abs_risk.rds")
p2521_rt <- readRDS("data/public/p2521_metastatic_risk_table.rds")
p212_ar <- p212_ar +
  labs(title = "NFCS - Aggressive PCa",
  y = "Cumulative Incidence (%)",
  x = "Age",
      color = "PHS percentile",
      fill = "PHS percentile") +
  lims(x = c(45, 90)) +
  theme_classic()
p212_rt <- p212_rt + labs(x = NULL) +
  lims(x = c(45, 90))
p2521_ar <- p2521_ar +
  labs(title = "HUSK - Aggressive PCa",
  y = "Cumulative Incidence (%)",
  x = "Age",
      color = "PHS percentile",
      fill = "PHS percentile") +
  lims(x = c(45, 90)) +
  theme_classic()
p2521_rt <- p2521_rt + labs(x = NULL) +
  lims(x = c(45, 90))
p212 <- hazrd::attach_risk_table(p212_ar, p212_rt)
p2521 <- hazrd::attach_risk_table(p2521_ar, p2521_rt)
metastatic <- (p2521 | p212) +
  # patchwork::plot_annotation(tag_levels = "A") +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("figures_and_tables/metastatic_km_curve.pdf", metastatic,
  height = 8 * 0.618, width = 8
)
ggsave("figures_and_tables/metastatic_km_curve.png", metastatic,
  height = 8 * 0.618, width = 8
)
ggsave("figures_and_tables/metastatic_km_curve.tiff", metastatic,
  height = 8 * 0.618, width = 8, dpi = 300
)

###############################################################
# WGS HOXB13 K-M Curve figure
###############################################################
wgs_phs_ar <- readRDS("data/public/p212_wgs_phs_abs_risk.rds")
wgs_phs_rt <- readRDS("data/public/p212_wgs_phs_risk_table.rds")
wgs_hxb_ar <- readRDS("data/public/p212_wgs_hoxb13_abs_risk.rds")
wgs_hxb_rt <- readRDS("data/public/p212_wgs_hoxb13_risk_table.rds")

wgs_phs_ar <-
  wgs_phs_ar +
  labs(x = "Age",
      color = "PHS percentile",
      fill = "PHS percentile") +
  lims(x = c(45, 90)) +
  theme_classic()
wgs_phs_rt <- wgs_phs_rt + labs(x = NULL) +
  lims(x = c(45, 90))
wgs_hxb_ar <-
  wgs_hxb_ar +
  labs(y = "Cumulative Incidence (%)",
       x = "Age") +
  guides(color = guide_legend(title = "HOXB13"),
         fill = guide_legend(title = "HOXB13")) +
  lims(x = c(45, 90)) +
  labs(title = "NFCS - Any-stage PCa") +
  theme_classic()
wgs_hxb_rt <- wgs_hxb_rt + labs(x = NULL) +
  lims(x = c(45, 90))

p212 <- hazrd::attach_risk_table(p212_ar, p212_rt)
p2521 <- hazrd::attach_risk_table(p2521_ar, p2521_rt)

any_pca <- (p2521 | p212) +
  # patchwork::plot_annotation(tag_levels = "A") +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")


wgs_phs_ar <- wgs_phs_ar +
  labs(title = "NFCS - Any-stage PCa", y = "Cumulative Incidence (%)")
# wgs_hxb_ar <- wgs_hxb_ar +
#   labs(y = "Cumulative Incidence (%)")
wgs_km <-
  (hazrd::attach_risk_table(wgs_phs_ar, wgs_phs_rt) |
   hazrd::attach_risk_table(wgs_hxb_ar, wgs_hxb_rt)) +
  # plot_annotation(tag_levels = "A") +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave("figures_and_tables/p212_wgs_km_curve.pdf", wgs_km,
  height = 8 * 0.618, width = 8
)
ggsave("figures_and_tables/p212_wgs_km_curve.png", wgs_km,
  height = 8 * 0.618, width = 8
)
ggsave("figures_and_tables/p212_wgs_km_curve.tiff", wgs_km,
  height = 8 * 0.618, width = 8, dpi = 300
)

###############################################################
# Stratified incidence rate figure
###############################################################
anypca <- readRDS("data/public/p2521_stratified_incidence.rds")
metastatic <- readRDS("data/public/p2521_stratified_incidence_metastatic.rds")
anypca <- anypca +
  labs(title = "HUSK - Any-stage PCa", y = "Cumulative Incidence (%)") +
  # lims(x = c(45, 90)) +
  theme_classic()

metastatic <- metastatic +
  labs(title = "HUSK - Aggressive PCa", y = "Cumulative Incidence (%)") +
  # lims(x = c(45, 90)) +
  theme_classic()

stratified <-
  anypca +
    metastatic +
    # patchwork::plot_annotation(tag_levels = "A") +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

ggsave("figures_and_tables/stratified_incidence.pdf", stratified,
  height = 6 * 0.618, width = 8
)
ggsave("figures_and_tables/stratified_incidence.png", stratified,
  height = 6 * 0.618, width = 8
)
ggsave("figures_and_tables/stratified_incidence.tiff", stratified,
  height = 6 * 0.618, width = 8, dpi = 300
)

###############################################################
# Supplementary Figure 1
###############################################################
p212_zph_tt <- readRDS("data/public/p212_zph_tt.rds")
p2521_zph_tt <- readRDS("data/public/p2521_zph_tt.rds")
p212_zph_tt_aggr <- readRDS("data/public/p212_zph_tt_aggr.rds")
p2521_zph_tt_aggr <- readRDS("data/public/p2521_zph_tt_aggr.rds")

p2521_zph_tt <- p2521_zph_tt + labs(title = "HUSK - Any-stage PCa", x = "Age")
p212_zph_tt <- p212_zph_tt + labs(title = "NFCS - Any-stage PCa", x = "Age")
p2521_zph_tt_aggr <- p2521_zph_tt_aggr +
  labs(title = "HUSK - Aggressive PCa", x = "Age")
p212_zph_tt_aggr <- p212_zph_tt_aggr +
  labs(title = "NFCS - Aggressive PCa", x = "Age")

zph_tt <-
  p2521_zph_tt +
    p212_zph_tt +
    p2521_zph_tt_aggr +
    p212_zph_tt_aggr +
    plot_annotation(tag_levels = "A") &
    theme(aspect.ratio = 1)

ggsave("figures_and_tables/zph_tt.pdf", zph_tt,
  height = 8, width = 8
)
ggsave("figures_and_tables/zph_tt.png", zph_tt,
  height = 8, width = 8
)
ggsave("figures_and_tables/zph_tt.tiff", zph_tt,
  height = 8, width = 8, dpi = 300
)

###############################################################
# Supplementary Figure 2
###############################################################
p2521_anypca <- readRDS("data/public/p2521_cumulative_comparison.rds")
p2521_metastatic <-
  readRDS("data/public/p2521_cumulative_comparison_metastatic.rds")
p212_anypca <- readRDS("data/public/p212_cumulative_comparison.rds")
p212_metastatic <-
  readRDS("data/public/p212_cumulative_comparison_metastatic.rds")

p2521_anypca <- p2521_anypca +
  labs(y = "Cumulative Incidence (%)", title = "HUSK - Any-stage PCa") +
  lims(x = c(45, 90), y = c(0, 30))
p212_anypca <- p212_anypca +
  labs(y = "Cumulative Incidence (%)", title = "NFCS - Any-stage PCa") +
  lims(x = c(45, 90), y = c(0, 30))
p2521_metastatic <- p2521_metastatic +
  labs(y = "Cumulative Incidence (%)", title = "HUSK - Aggressive PCa") +
  lims(x = c(45, 90), y = c(0, 6))
p212_metastatic <- p212_metastatic +
  labs(y = "Cumulative Incidence (%)", title = "NFCS - Aggressive PCa") +
  lims(x = c(45, 90), y = c(0, 6))

comparison <-
  p2521_anypca +
  p212_anypca +
  p2521_metastatic +
  p212_metastatic +
  plot_annotation(tag_levels = "A")

ggsave("figures_and_tables/cumulative_comparison.pdf", comparison,
  height = 8, width = 8
)
ggsave("figures_and_tables/cumulative_comparison.png", comparison,
  height = 8, width = 8
)
ggsave("figures_and_tables/cumulative_comparison.tiff", comparison,
  height = 8, width = 8, dpi = 300
)

###############################################################
# Supplementary Table - cox zph test
###############################################################
table_zph <-
  p212_zph %>%
  bind_rows(p2521_zph) %>%
  mutate(cohort = if_else(cohort == "p212", "NFCS", "HUSK")) %>%
  mutate(PCa = str_to_title(PCa)) %>%
  mutate(chisq = round(chisq, 2)) %>%
  mutate(p = pvalr(p)) %>%
  arrange(desc(PCa), cohort) %>%
  select(PCa, Cohort = cohort, ChiSq = chisq, df, `P-value` = p)

write_tsv(table_zph, "figures_and_tables/table_zph.tsv")

###############################################################
# Supplementary Table - fishers tests
###############################################################
table_fisher <-
  fisher %>%
  mutate(OR = sprintf("%.2f (%.2f)", estimate, conf.low)) %>%
  select(-conf.low, -conf.high) %>%
  mutate(gene = str_remove(gene, "_carrier")) %>%
  mutate(gene = if_else(gene == "any", "Any", gene)) %>%
  mutate(p.value = pvalr(p.value), p.adjust = pvalr(p.adjust)) %>%
  select(RPV = gene, `OR (LCB)` = OR, `P-value` = p.value, `Adj. P-value` = p.adjust)

write_tsv(table_fisher, "figures_and_tables/table_fisher.tsv")

###############################################################
# Supplementary Table - time-transform
###############################################################
coxph_results <-
  p2521_coxph_results %>%
  rbind(p212_coxph_results) %>%
  filter(model %in% c("phs_tt")) %>%
  select(term, value = estimate, conf.low, conf.high, cohort)

coxph_results$index <- "HR"

sample_size <-
  p2521_sample_size %>%
  rbind(p212_sample_size) %>%
  filter(group == "any_pca") %>%
  select(cohort, cases)

table_tt <-
  coxph_results %>%
  pivot_wider(values_from = c(value, conf.low, conf.high), names_from = "term") %>%
  mutate(value_phs = sprintf("%.2f (%.2f, %.2f)", value_phs, conf.low_phs, conf.high_phs)) %>%
  mutate(`value_tt(phs)` = sprintf("%.3f (%.3f, %.3f)", `value_tt(phs)`, `conf.low_tt(phs)`, `conf.high_tt(phs)`)) %>%
  select(-conf.low_phs, -`conf.low_tt(phs)`, -conf.high_phs, -`conf.high_tt(phs)`) %>%
  pivot_longer(starts_with("value_")) %>%
  mutate(name = str_remove(name, "value_")) %>%
  mutate(name = if_else(name == "phs", "PHS", "tt(PHS)")) %>%
  left_join(sample_size, by = join_by(cohort)) %>%
  select(Cohort = cohort, Cases = cases, Term = name, `HR (95% CI)` = value)

write_tsv(table_tt, "figures_and_tables/table_tt.tsv")

###############################################################
# Supplementary Table - wgs subset
###############################################################
wgs_results <-
  wgs_results %>%
  left_join(wgs_cindex)

wgs_results$model <-
  c(
    "PHS", "HOXB13", "PHS + HOXB13", "PHS + HOXB13", "PHS x HOXB13",
    "PHS x HOXB13", "PHS x HOXB13"
  )

wgs_results$term <-
  c("PHS", "HOXB13", "PHS", "HOXB13", "PHS", "HOXB13", "PHS:HOXB13")

table_gws <-
  wgs_results %>%
  mutate(HR = sprintf("%.2f (%.2f, %.2f)", estimate, conf.low, conf.high)) %>%
  mutate(cindex = sprintf(
    "%.3f (%.3f)",
    concordance,
    std.error.concordance
  )) %>%
  mutate(p.value = pvalr(p.value)) %>%
  select(
    Model = model, Term = term, `HR (95% CI)` = HR, `P-value` = p.value,
    `C-index (SE)` = cindex
  ) %>%
  filter(Model %in% c("PHS", "HOXB13")) %>%
  select(-Term)

write_tsv(table_gws, "figures_and_tables/table_wgs.tsv")

###############################################################
# Supplementary Table - all pathogenic variants
###############################################################
all_variants <- read_tsv("data/public/all_variants.tsv")

table_variants <-
  all_variants %>%
  filter(FINAL_CLASSIFICATION %in% c("Pathogenic")) %>%
  rename(
    Symbol = SYMBOL,
    "Genomic Change" = GENOMIC_CHANGE,
    Alteration = ALTERATION,
    "dbSNP RSID" = DBSNP_RSID,
    "ClinVar Classification" = FINAL_CLASSIFICATION,
    N = n
  )

write_tsv(table_variants, "figures_and_tables/table_variants.tsv")

###############################################################
# Supplementary Table - all nonpathogenix hoxb13 variants
###############################################################

table_hoxb13 <-
  all_variants %>%
  filter(SYMBOL == "HOXB13") %>%
  rename(
    Symbol = SYMBOL,
    "Genomic Change" = GENOMIC_CHANGE,
    Alteration = ALTERATION,
    "dbSNP RSID" = DBSNP_RSID,
    "ClinVar Classification" = FINAL_CLASSIFICATION,
    N = n
  )

write_tsv(table_hoxb13, "figures_and_tables/table_hoxb13.tsv")

###############################################################
# Supplementary Table - age summary
###############################################################
p212_age_summary     <- read_tsv("data/public/p212_age_summary.tsv")
p212_age_summary_wgs <- read_tsv("data/public/p212_age_summary_wgs.tsv")
p2521_age_summary    <- read_tsv("data/public/p2521_age_summary.tsv")

p212_age_summary$cohort     <- "NFCS"
p212_age_summary_wgs$cohort <- "NFCS (WGS)"
p2521_age_summary$cohort    <- "HUSK"

table_age_summary <-
  p2521_age_summary %>%
  rbind(p212_age_summary) %>%
  rbind(p212_age_summary_wgs) %>%
  mutate(value = case_when(
    variable == "n_total"         ~ sprintf("%.0f", n),
    str_detect(variable, "^pct_") ~ sprintf("%.0f (%.1f%%)", n, pct),
    variable == "birth_year"      ~ sprintf("%.0f (%.0f, %.0f)", median, q25, q75),
    TRUE                          ~ sprintf("%.2f (%.2f, %.2f)", median, q25, q75)
  )) %>%
  mutate(variable = case_match(
    variable,
    "birth_year"           ~ "Birth Year",
    "age_first_follow_up"  ~ "Age at First Follow-up",
    "age_diagnosis"        ~ "Age at Diagnosis",
    "age_of_metastasis"    ~ "Age at Aggressive Diagnosis",
    "age_lastobs_censored" ~ "Age at Last Follow-up (controls)",
    "phs"                  ~ "PHS (z-score)",
    "n_total"              ~ "n",
    "pct_pca"              ~ "PCa cases, n (%)",
    "pct_aggressive"       ~ "Aggressive PCa cases, n (%)",
    "pct_censored"         ~ "Censored, n (%)",
    "pct_death_before_dx"  ~ "Died before diagnosis, n (%)"
  )) %>%
  select(cohort, Variable = variable, `Median (Q1, Q3)` = value) %>%
  pivot_wider(names_from = cohort, values_from = `Median (Q1, Q3)`) %>%
  mutate(across(where(is.character), ~ if_else(str_detect(.x, "^NA"), "\u2014", .x))) %>%
  add_row(
    Variable = "Recruitment Period",
    HUSK = "1992\u20131999", NFCS = "1980\u20131995", `NFCS (WGS)` = "\u2014",
    .after = 1
  )

write_tsv(table_age_summary, "figures_and_tables/table_age_summary.tsv")
