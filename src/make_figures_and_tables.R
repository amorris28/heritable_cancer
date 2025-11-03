library(ggplot2)
library(tidyverse)
library(patchwork)
#library(kableExtra)

source("src/config.R")
source("src/functions.R")
model = "PHS601"

# Import results
p2521_PHS601 = readr::read_tsv("data/public/p2521_PHS601_any_pca_performance.tsv")
p212_PHS601 = readr::read_tsv("data/public/p212_PHS601_any_pca_performance.tsv")
p2521_PHS601_meta = readr::read_tsv("data/public/p2521_PHS601_metastatic_performance.tsv")
p212_PHS601_meta = readr::read_tsv("data/public/p212_PHS601_metastatic_performance.tsv")

p212_zph = readr::read_tsv("data/public/p212_PHS601_phs_zph_tidy.tsv")
p2521_zph = readr::read_tsv("data/public/p2521_PHS601_phs_zph_tidy.tsv")
p2521_sample_size = readr::read_tsv("data/public/p2521_PHS601_sample_size.tsv")
p212_sample_size = readr::read_tsv("data/public/p212_PHS601_sample_size.tsv")
p212_coxph_results = readr::read_tsv("data/public/p212_PHS601_coxph.tsv")
p2521_coxph_results = readr::read_tsv("data/public/p2521_PHS601_coxph.tsv")
wgs_summary = readr::read_tsv("data/public/p212_PHS601_wgs_summary.tsv")
fisher = readr::read_tsv("data/public/p212_PHS601_fisher.tsv")
wgs_phs = readr::read_tsv("data/public/p212_PHS601_wgs_phs_performance.tsv")
wgs_phs_95 = readr::read_tsv("data/public/p212_PHS601_wgs_phs_performance_95.tsv")
wgs_hoxb13 = readr::read_tsv("data/public/p212_PHS601_wgs_hoxb13_performance.tsv")
binomial_phs_hoxb13 = readr::read_tsv("data/public/p212_PHS601_binomial_phs_hoxb13.tsv")
wgs_results = readr::read_tsv("data/public/p212_PHS601_wgs.tsv")

###############################################################
# PHS K-M Curve Figure for Any-PCa
###############################################################

p212 <- readRDS("data/public/p212_phs601_any_pca_km_curve.rds")
p2521 <- readRDS("data/public/p2521_phs601_any_pca_km_curve.rds")
# p212 = p212 + 
#     labs(y = "Proportion cancer-free, any stage (NFCS)") +
    #   lims(y = c(0.25, 1))
# p2521 = p2521 + 
#     labs(y = "Proportion cancer-free, any stage (HUSK)") +
    #   lims(y = c(0.25, 1))
any_pca = p2521 +
    p212 +
    plot_annotation(tag_levels = 'A') +
    plot_layout(guides='collect') &
    theme(legend.position = "bottom")

ggsave("figures_and_tables/km_curve.pdf", any_pca,
       height = 8 * 0.618, width = 8)
ggsave("figures_and_tables/km_curve.png", any_pca,
       height = 8 * 0.618, width = 8)
#saveRDS(any_pca, here(paste0(results, paste(model, "km_curve.rds", sep = "_"))))

###############################################################
# PHS Any-PCa and Meta-PCa Table
###############################################################
p2521_sample_size$cohort = "HUSK"
p212_sample_size$cohort = "NFCS"
p2521_PHS601$cohort = "HUSK"
p2521_PHS601$PCa = "Any"
p212_PHS601$cohort = "NFCS"
p212_PHS601$PCa = "Any"
p2521_PHS601_meta$cohort = "HUSK"
p2521_PHS601_meta$PCa = "Meta"
p212_PHS601_meta$cohort = "NFCS"
p212_PHS601_meta$PCa = "Meta"

p2521_sample_size %>%
    rbind(p212_sample_size) %>%
    mutate(PCa = if_else(group == "any_pca", "Any", if_else(group == "metastatic", "Meta", NA))) ->
    sample_size

rbind(p2521_PHS601, p212_PHS601, p2521_PHS601_meta, p212_PHS601_meta) %>%
    filter(index %in% c("HRSD", "HR80_20", "C_Index")) %>%
    mutate(value = sprintf("%.2f (%.2f,%.2f)", value, conf.low, conf.high)) %>%
    select(-conf.low, -conf.high) %>%
    pivot_wider(names_from = "index", values_from = "value") %>%
    left_join(sample_size, by = join_by(cohort, PCa)) %>%
    select(PCa, cohort, n = total, cases, HRSD, HR80_20, C_Index) ->
    table_1

write_tsv(table_1, "figures_and_tables/table_1.tsv")

###############################################################
# PHS K-M Curve Figure for Meta-PCa
###############################################################
p212 <- readRDS("data/public/p212_phs601_metastatic_inc_curve.rds")
p2521 <- readRDS("data/public/p2521_phs601_metastatic_inc_curve.rds")
# p212 = p212 + labs(y = "Proportion cancer-free, metastatic (NFCS)") +
#     lims(y = c(0.80, 1.0))
# p2521 = p2521 + labs(y = "Proportion cancer-free, metastatic (HUSK)") +
#     lims(y = c(0.80, 1.0))
metastatic = p2521 + 
    p212 +
    plot_annotation(tag_levels = 'A') +
    plot_layout(guides='collect') &
    theme(legend.position = "bottom")

ggsave("figures_and_tables/metastatic_km_curve.pdf", metastatic,
           height = 8 * 0.618, width = 8)
ggsave("figures_and_tables/metastatic_km_curve.png", metastatic,
           height = 8 * 0.618, width = 8)
#saveRDS(metastatic, here(paste0(results, paste(model, "metastatic_km_curve.rds", sep = "_"))))


###############################################################
# WGS HOXB13 K-M Curve figure 
###############################################################

# library(patchwork)
# margin = theme(plot.margin = unit(c(2,2,2,2), "cm"))
# wgs_hoxb13 <- readRDS("data/results/p212_PHS601_wgs_km_curve.rds")
# wgs_phs_km <- readRDS("data/results/p212_PHS601_phs_wgs_km_curve.rds")
# wgs_hoxb13 = wgs_hoxb13 + 
#     ggtitle(NULL) +
#     theme(legend.position = "bottom") 
# wgs_phs_km = wgs_phs_km + 
#     ggtitle(NULL) +
#     theme(legend.position = "bottom") 
# #wgs_hoxb13
# #wgs_phs_km
# wgs_hoxb13 + wgs_phs_km  +
#     plot_annotation(tag_levels = 'A')
# # gridExtra::grid.arrange(wgs_hoxb13, wgs_phs_km,
# #   widths = c(5, 1, 5),
# #   layout_matrix = rbind(c(1, NA, 2))
# # )

wgs_km <- readRDS("data/public/p212_phs601_wgs_km_curve.rds")
wgs_km = wgs_km & labs(y = "Proportion cancer-free, any stage (NFCS)")
ggsave("figures_and_tables/p212_phs601_wgs_km_curve.pdf", wgs_km,
           height = 8 * 0.618, width = 8)
ggsave("figures_and_tables/p212_phs601_wgs_km_curve.png", wgs_km,
           height = 8 * 0.618, width = 8)


###############################################################
# WGS HOXB13 Table
###############################################################

# rbind(p2521_PHS601, p212_PHS601, p2521_PHS601_meta, p212_PHS601_meta) %>%
wgs_hoxb13$model = "hoxb13"
wgs_phs_95$model = "phs"
# options(knitr.kable.NA = '-')
rbind(wgs_hoxb13, wgs_phs_95) %>%
    filter(index %in% c("HR", "C_Index")) %>%
    mutate(value = sprintf("%.2f (%.2f, %.2f)", value, conf.low, conf.high)) %>%
    select(-conf.low, -conf.high) %>%
    pivot_wider(names_from = "index", values_from = "value") %>%
    select(model, HR, C_Index) %>%
    arrange(desc(model)) ->
    table_2
#write_tsv(table_2, "table_2.tsv")


###############################################################
# Stratified incidence rate figure
###############################################################

anypca <- readRDS("data/public/p2521_phs601_stratified_incidence.rds")
metastatic <- readRDS("data/public/p2521_phs601_stratified_incidence_metastatic.rds")
anypca = anypca + 
    labs(y = "Predicted cumulative incidence,\nany stage, %", title = "") +
    lims(x = c(40, 90), y = c(0, 40))
metastatic = metastatic + 
    labs(y = "Predicted cumulative incidence,\nmetastatic, %", title = "") +
    lims(x = c(40, 90), y = c(0, 8))

stratified = 
    anypca + 
    metastatic +
    plot_annotation(tag_levels = 'A') +
    plot_layout(guides='collect') &
    theme(legend.position = "bottom") #&
    # theme(aspect.ratio=1)

ggsave("figures_and_tables/stratified_incidence.pdf", stratified,
       height = 8 * 0.618, width = 8)
ggsave("figures_and_tables/stratified_incidence.png", stratified,
       height = 8 * 0.618, width = 8)
#saveRDS(stratified, here(paste0(results, paste(model, "stratified_incidence.rds", sep = "_"))))


###############################################################
#Supplementary Figure 1
###############################################################
# change these to p2521:
p2521_anypca <- readRDS("data/public/p2521_phs601_cumulative_comparison.rds")
p2521_metastatic <- readRDS("data/public/p2521_phs601_cumulative_comparison_metastatic.rds")
p212_anypca <- readRDS("data/public/p212_phs601_cumulative_comparison.rds")
p212_metastatic <- readRDS("data/public/p212_phs601_cumulative_comparison_metastatic.rds")

p2521_anypca = p2521_anypca + 
    labs(y = "Cumulative incidence,\nany stage, % (HUSK)", title = "") +
    lims(x = c(40, 90), y = c(0, 30))
p212_anypca = p212_anypca + 
    labs(y = "Cumulative incidence,\nany stage, % (NFCS)", title = "") +
    lims(x = c(40, 90), y = c(0, 30))
p2521_metastatic = p2521_metastatic + 
    labs(y = "Cumulative incidence,\nmetastatic, % (HUSK)", title = "") +
    lims(x = c(40, 90), y = c(0, 6))
p212_metastatic = p212_metastatic + 
    labs(y = "Cumulative incidence,\nmetastatic, % (NFCS)", title = "") +
    lims(x = c(40, 90), y = c(0, 6))
# p2521_anypca = p2521_anypca + lims(y = c(0, 30000)) + labs(y = "Cumulative PCa incidence rate,\nany stage, per 100,000 (HUSK)")
# p212_anypca = p212_anypca + lims(y = c(0, 30000)) + labs(y = "Cumulative PCa incidence rate,\nany stage, per 100,000 (NFCS)")
# p2521_metastatic = p2521_metastatic + lims(y = c(0, 6500)) + labs(y = "Cumulative PCa incidence rate,\nmetastatic, per 100,000 (HUSK)")
# p212_metastatic = p212_metastatic + lims(y = c(0, 6500)) + labs(y = "Cumulative PCa incidence rate,\nmetastatic, per 100,000 (NFCS)")


comparison = 
    p2521_anypca + 
    p212_anypca +
    p2521_metastatic +
    p212_metastatic +
    plot_annotation(tag_levels = 'A')

ggsave("figures_and_tables/cumulative_comparison.pdf", comparison,
       height = 8, width = 8)
ggsave("figures_and_tables/cumulative_comparison.png", comparison,
       height = 8, width = 8)
#saveRDS(comparison, here(paste0(results, paste(model, "cumulative_comparison.rds", sep = "_"))))
#knitr::include_graphics(here("data/results/PHS601_cumulative_comparison.png"))


#file.copy("data/public/p212_phs601_wgs_km_curve.pdf", "figures_and_tables/p212_phs601_wgs_km_curve.pdf")
#file.copy("data/public/p212_phs601_wgs_km_curve.png", "figures_and_tables/p212_phs601_wgs_km_curve.png")


