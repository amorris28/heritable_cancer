library(tidyverse)
library(stringr)
library(hazrd)
library(survival)
library(patchwork)
source("src/00_config/config.R")
source("src/00_config/functions.R")

###########################################
# Output flags
bootstrap_iters = 999
perform_bool = TRUE

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

args = commandArgs(trailingOnly=TRUE)

model_file = args[1] # PHS file
metadata_file = args[2] # Phenotype data
inverse = args[3]
tsd_proj = args[4]

if (Sys.info()["nodename"] == "p2521-hpc-01.tsd.usit.no") {
    metadata_file = "data/processed/husk_phenotypes.tsv"
    tsd_proj = "p2521"
} else if (Sys.info()["nodename"] == "p212-hpc-01.tsd.usit.no") {
    metadata_file = "data/processed/heritable_cancer_phenotypes.tsv" # Phenotype data
    tsd_proj = "p212"
}

model_file = "data/processed/PHS601_phs" # PHS file
model = gsub("(.+?)(\\_.*)", "\\1", basename(model_file))
inverse = FALSE
#########################################
# Import PHS data and metadata
phs_data = import_phs_data(model_file, inverse)

if (Sys.info()["nodename"] == "p2521-hpc-01.tsd.usit.no") {
phs_data2 = read_table("src/xx_merge_snps/data_prs_husk.profile")
} else if (Sys.info()["nodename"] == "p212-hpc-01.tsd.usit.no") {
phs_data2 = read_table("src/xx_merge_snps/data_prs_heritable_cancer.profile")
}
phs_data2 %>%
  select(bfile_id = IID, phs = SCORESUM) ->
  phs_data2

metadata = read_tsv(metadata_file) # Phenotype data

combined_data = process_combined_data(metadata, phs_data2)

#########################################
# Analyze PHS data

# Coxph
combined_data %>%
  filter(status == 1) %>%
  summarize(min_age = min(age)) # 40ish

combined_data %>%
  filter(age >= 40) %>%
  mutate(age_minus_40 = age - 40) ->
  combined_data

fit = coxph(Surv(age, status) ~ phs, data = combined_data)
zph_fit = cox.zph(fit)

plot(zph_fit)
summary(fit)

zph_fit %>%
  print() %>%
  as_tibble ->
phs_zph_tidy

write_tsv(phs_zph_tidy, paste0(results, paste(tsd_proj, model, "phs_zph_tidy.tsv", sep = "_")))

any_pca_coxph = broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)

# Time transformation
phs_fit_tt = coxph(Surv(age_minus_40, status) ~ phs + tt(phs), combined_data,
                  tt = function(x, t, ...) x * t)
phs_tt_tidy = broom::tidy(phs_fit_tt, exponentiate = TRUE, conf.int = TRUE)

summary(phs_fit_tt)
# not significant
coef(phs_fit_tt)[[1]]
zph_tbl = tibble(time = zph_fit$time,
      beta_t = zph_fit$y[, "phs"])

tt_comp_f = function(t) {
    coef(phs_fit_tt)[[1]] + coef(phs_fit_tt)[[2]] * t
    # 0.894841 - 0.011519 * t
}
zph_tbl$tt_comp = tt_comp_f(zph_tbl$time - 40)

ggplot(zph_tbl, aes(x = time, y = beta_t)) +
    geom_point() +
    stat_smooth() +
    geom_line(aes(y = tt_comp), color = "red") +
    theme_classic()

# Performance metrics

if (perform_bool) {
  perform = get_performance_metrics(combined_data, bootstrap_iters, model)
  write_tsv(perform, paste0(results, paste(tsd_proj, model, "any_pca_performance.tsv", sep = "_")))
}

# Plot KM Curves
km_curve = plot_km_curves(combined_data)

ggsave(paste0(results, paste(tsd_proj, model, "km_curve.pdf", sep = "_")), km_curve,
       height = 6 * 0.618, width = 6)
ggsave(paste0(results, paste(tsd_proj, model, "km_curve.png", sep = "_")), km_curve,
       height = 6 * 0.618, width = 6)
saveRDS(km_curve, paste0(results, paste(tsd_proj, model, "any_pca_km_curve.rds", sep = "_")))

# Time to metastatic diseaes
metastasis_fit = coxph(Surv(age_of_metastasis, metastasis) ~ phs, combined_data)
metastatic_coxph = broom::tidy(metastasis_fit, exponentiate = TRUE, conf.int = TRUE)

# Performance metrics

if (perform_bool) {
  perform = get_performance_metrics(combined_data, bootstrap_iters, model, age = "age_of_metastasis", status = "metastasis")
  write_tsv(perform, paste0(results, paste(tsd_proj, model, "metastatic_performance.tsv", sep = "_")))
}

# KM Curves for metastatic disease
km_curve = plot_km_curves(combined_data, age = "age_of_metastasis", status = "metastasis")
km_curve = km_curve + lims(y = c(0.80, 1.0))
ggsave("kmcurve_metastasis_strict.png", km_curve, height = 6 * 0.618, width = 6)

ggsave(paste0(results, paste(tsd_proj, model, "km_curve_metastatic.pdf", sep = "_")), km_curve,
       height = 6 * 0.618, width = 6)
ggsave(paste0(results, paste(tsd_proj, model, "km_curve_metastatic.png", sep = "_")), km_curve,
       height = 6 * 0.618, width = 6)
saveRDS(km_curve, paste0(results, paste(tsd_proj, model, "metastatic_km_curve.rds", sep = "_")))


# Gather sample sizes
combined_data %>%
  select(phs, age, status) ->
  any_pca_sample_size
any_pca_sample_size_cases = sum(any_pca_sample_size$status)
any_pca_sample_size_total = nrow(any_pca_sample_size)
any_pca_sample_size_controls = any_pca_sample_size_total - any_pca_sample_size_cases


combined_data %>%
  select(phs, age_of_metastasis, metastasis) -> #%>%
  # drop_na() ->
  metastatic_sample_size
metastatic_sample_size_cases = sum(metastatic_sample_size$metastasis)
metastatic_sample_size_total = nrow(metastatic_sample_size)
metastatic_sample_size_controls = 
  metastatic_sample_size_total - 
  metastatic_sample_size_cases

sample_size = tibble(group = c("any_pca", "metastatic"),
       total = c(any_pca_sample_size_total, metastatic_sample_size_total),
       cases = c(any_pca_sample_size_cases, metastatic_sample_size_cases),
       controls = c(any_pca_sample_size_controls, metastatic_sample_size_controls))

write_tsv(sample_size, paste0(results, paste(tsd_proj, model, "sample_size.tsv", sep = "_")))

# Combine Coxph results
any_pca_coxph$model = "any_pca"
metastatic_coxph$model = "metastatic"
phs_tt_tidy$model = "phs_tt"

rbind(any_pca_coxph, 
      metastatic_coxph, 
      phs_tt_tidy) %>%
write_tsv(paste0(results, paste(tsd_proj, model, "coxph.tsv", sep = "_")))




