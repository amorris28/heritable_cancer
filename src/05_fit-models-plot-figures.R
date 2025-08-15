library(tidyverse)
library(stringr)
library(hazrd)
library(survival)
library(patchwork)
source("src/00_config.R")

###########################################
# Output flags
bootstrap_iters = 1
perform_bool = FALSE


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
model_file = "data/processed/PHS601_phs"
# model = gsub("^.*?_", "", basename(model_file))
# model = gsub("\\_.*", "", model) # filename for output
model = gsub("(.+?)(\\_.*)", "\\1", basename(model_file))

metadata_file = args[2] # Phenotype data
metadata_file = "data/processed/heritable_cancer_phenotypes.tsv"
# metadata_file = "data/processed/husk_phenotypes.tsv"
inverse = args[3]
inverse = FALSE
tsd_proj = args[4]
tsd_proj = "p212"
# tsd_proj = "p2521"



#########################################
# Import PHS data

phs_data = read_tsv(model_file, col_names = c("bfile_id", "phs"))

metadata = read_tsv(metadata_file) # Phenotype data

# colnames(metadata) = c("bfile_id", "age", "status")
phs_data %>%
  mutate(phs = as.vector(scale(phs, center = TRUE, scale = TRUE))) ->
  phs_data

if (inverse) {mutate(phs_data, phs = phs * -1) -> phs_data}

metadata %>%
  left_join(phs_data, by = "bfile_id") %>%
  filter(!(is.na(phs))) ->
  combined_data

# Remove samples with missing genotype data or phenotype data
combined_data %>%
  drop_na(phs, age, status) ->
  combined_data

#########################################
# Analyze PHS data


# PHS Histogram
phs_hist = 
	phs_hist(combined_data, normalize = TRUE) +
	labs(title = model) + theme(aspect.ratio=1)

# ggsave(paste0(results, paste(tsd_proj, model, "hist.png", sep = "_")), phs_hist, 
#        units = "px", height = 1300, width = 1300)

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

print(zph_fit)
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
# zph not defined for tt() fit
# zph_fit_tt = cox.zph(phs_fit_tt)
# print(zph_fit_tt)
# plot(zph_fit_tt)
summary(phs_fit_tt)
# not significant
coef(phs_fit_tt)[[1]]
zph_tbl = tibble(time = zph_fit$time,
      beta_t = zph_fit$y[,"phs"])

tt_comp_f = function(t) {
    coef(phs_fit_tt)[[1]] + coef(phs_fit_tt)[[2]] * t
    # 0.894841 - 0.011519 * t
}
zph_tbl$tt_comp = tt_comp_f(zph_tbl$time-40)

ggplot(zph_tbl, aes(x = time, y = beta_t)) +
    geom_point() + 
    stat_smooth() +
    geom_line(aes(y = tt_comp), color = "red") +
    theme_classic()


# Performance metrics
HR80_20 = get_hr(combined_data, CI = TRUE, boot = bootstrap_iters)
HR80_60 = get_hr(combined_data, 
                 lower_interval = c(0.2, 0.80),
                 upper_interval = c(0.8, 1.0),
                 CI = TRUE, 
                 boot = bootstrap_iters)

ages = c(60, 70, 80)
OR80_20 = lapply(ages, function(x) {get_or(combined_data, or_age = x, CI = TRUE, boot = bootstrap_iters)})
c_index = get_cindex(combined_data, CI = TRUE, boot = bootstrap_iters)

HR80_20$age = NA
HR80_60$age = NA
HR80_60$index = "HR80_60"
c_index$age = NA
indices = list(HR80_20, HR80_60, OR80_20[[1]], OR80_20[[2]], OR80_20[[3]], c_index)

perform = do.call(rbind.data.frame, lapply(indices, function(x){data.frame(index = x$index, 
                                                                 value = x$value, 
                                                                 conf.low = x$conf.low, 
                                                                 conf.high = x$conf.high, 
                                                                 age = x$age)}))
perform = cbind(model, perform)
if (perform_bool) {
  write_tsv(perform, paste0(results, paste(tsd_proj, model, "any_pca_performance.tsv", sep = "_")))
}

# KM Curves
curves = data.frame(curve0_20 = c(0,   0.2),
                    curve20_80 = c(0.2, 0.8),
                    curve80_100 = c(0.8, 1.0))

label_generator = function(x, y) {
  x = x * 100
  y = y * 100
  # out = paste0("PHS ", x, "-", y, "th centile")
  out = paste0(x, "-", y, "th")
  return(out)
}

km_curves = data.frame()
for (i in seq_len(length(curves))) {
  curven <- km_curve(data = combined_data,
                     interval = curves[[i]],
                     age_range = 40:100, 
                     scale = FALSE, 
                     inverse = FALSE)
  curven$label = label_generator(curves[1, i], curves[2, i])
  km_curves = rbind(km_curves, curven)
}


km_curve = ggplot(km_curves, aes(x = time, 
                      y = estimate,
                      ymin = conf.low,
                      ymax = conf.high,
                      col = label,
                      fill = label)) +
  geom_ribbon(alpha = 0.1,
              color = 0) +
  geom_step() +
  theme_classic() +
  xlim(40, 100) + 
  ylim(0, 1) +
  labs(x = "Age", y = "Disease-free Survival") +
  scale_color_brewer(palette = "Set1", direction = -1,
                     name = "PHS Centile") +
  scale_fill_brewer(palette = "Set1", direction = -1,
                    name = "PHS Centile")  + 
                    theme(aspect.ratio=1)

ggsave(paste0(results, paste(tsd_proj, model, "km_curve.pdf", sep = "_")), km_curve,
       height = 6 * 0.618, width = 6)
ggsave(paste0(results, paste(tsd_proj, model, "km_curve.png", sep = "_")), km_curve,
       height = 6 * 0.618, width = 6)

saveRDS(km_curve, paste0(results, paste(tsd_proj, model, "any_pca_km_curve.rds", sep = "_")))

# Time to metastatic diseaes

combined_data %>%
  mutate(metastasis = if_else(is.na(metastasis), 0, metastasis)) %>%
  mutate(age_of_metastasis = if_else(is.na(age_of_metastasis), age_lastobs, age_of_metastasis)) ->
  combined_data
# Coxph

metastasis_fit = coxph(Surv(age_of_metastasis, metastasis) ~ phs, combined_data)
metastatic_coxph = broom::tidy(metastasis_fit, exponentiate = TRUE, conf.int = TRUE)

# combined_data %>%
#   select(age_lastobs, age_diagnosis, age, age_of_metastasis, age_at_death) %>%
#   rowwise() %>%
#   mutate(max = max(c_across(age_lastobs:age_at_death), na.rm = T)) ->
#   ages

# Performance metrics
HR80_20 = get_hr(combined_data,
                 age = "age_of_metastasis",
                 status = "metastasis",
                 CI = TRUE, 
                 boot = bootstrap_iters)
HR80_60 = get_hr(combined_data,
                 age = "age_of_metastasis",
                 status = "metastasis",
                 lower_interval = c(0.2, 0.80),
                 upper_interval = c(0.8, 1.0),
                 CI = TRUE, 
                 boot = bootstrap_iters)

ages = c(60, 70, 80)
OR80_20 = lapply(ages, function(x) {get_or(combined_data,
                                           age = "age_of_metastasis",
                                           status = "metastasis", 
                                           or_age = x, 
                                           CI = TRUE, 
                                           boot = bootstrap_iters)})
c_index = get_cindex(combined_data,
                     age = "age_of_metastasis",
                     status = "metastasis", 
                     CI = TRUE, 
                     boot = bootstrap_iters)

HR80_20$age = NA
HR80_60$age = NA
HR80_60$index = "HR80_60"
c_index$age = NA
indices = list(HR80_20, HR80_60, OR80_20[[1]], OR80_20[[2]], OR80_20[[3]], c_index)

perform = do.call(rbind.data.frame, lapply(indices, function(x){data.frame(index = x$index, 
                                                                 value = x$value, 
                                                                 conf.low = x$conf.low, 
                                                                 conf.high = x$conf.high, 
                                                                 age = x$age)}))
perform = cbind(model, perform)
if (perform_bool) {
  write_tsv(perform, paste0(results, paste(tsd_proj, model, "metastatic_performance.tsv", sep = "_")))
}
# KM Curves for metastatic disease
curves = data.frame(curve0_20 = c(0,   0.2),
                    curve20_80 = c(0.2, 0.8),
                    curve80_100 = c(0.8, 1.0))


km_curves = data.frame()
for (i in seq_len(length(curves))) {
  curven <- km_curve(data = combined_data,
                     age = "age_of_metastasis",
                     status = "metastasis",
                     interval = curves[[i]],
                     age_range = 40:100, 
                     scale = FALSE, 
                     inverse = FALSE)
  curven$label = label_generator(curves[1, i], curves[2, i])
  km_curves = rbind(km_curves, curven)
}


km_curve = ggplot(km_curves, aes(x = time, 
                      y = estimate,
                      ymin = conf.low,
                      ymax = conf.high,
                      col = label,
                      fill = label)) +
  geom_ribbon(alpha = 0.1,
              color = 0) +
  geom_step() +
  theme_classic() +
  xlim(40, 100) + 
  ylim(0.8, 1) +
  labs(x = "Age", y = "Disease-free Survival") +
  scale_color_brewer(palette = "Set1", direction = -1,
                     name = "PHS Centile") +
  scale_fill_brewer(palette = "Set1", direction = -1,
                    name = "PHS Centile")  + theme(aspect.ratio=1)

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
  select(phs, age_of_metastasis, metastasis) %>%
  drop_na() ->
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
#####################################################################################
# RPVs
#####################################################################################
# BRCA2, HOXB13
if (tsd_proj == "p212") {

  combined_data %>%
    select(status, age, phs, age_of_metastasis, metastasis, wgs, brca2_carrier:hoxb13_snp) %>%
    mutate(across(brca2_carrier:atm_carrier, ~ if_else(is.na(.x), 0, .x))) ->
    carriers

  carriers %>%
    filter(wgs == 1) %>%
    # any_carrier is 1 for people who have any of the listed RPVs
    mutate(any_carrier = case_when(if_any(brca2_carrier:atm_carrier, ~ .x == 1) ~ 1, .default = 0)) ->#%>%
    # all_carrier is the sum of all RPVs that each individual has (in this case, nobody has
    # more than one)
    # mutate(all_carrier = rowSums(across(brca2_carrier:atm_carrier))) ->
    wgs

  wgs %>%
    select(status, wgs:atm_carrier, any_carrier) %>%
    summarize(across(everything(), list(sum))) ->
    wgs_summary

  write_tsv(wgs_summary, paste0(results, paste(tsd_proj, "wgs_summary.tsv", sep = "_")))


  wgs %>%
    select(-msh2_carrier, -msh6_carrier) ->
    wgs
    # 0 subjects with these mutations

  #Fisher test
  variants = c("brca2", "hoxb13", "brca1", "mlh1", "atm", "any")#, "all")
  fisher_results = list()
  for (i in 1:length(variants)) {
    variant = paste0(variants[i], "_carrier")
    print(variant)
    wgs %>%
      select(status, all_of(variant)) %>%
      mutate(case_factor = factor(status, levels = c(1, 0), labels = c("Cases", "Controls"))) %>%
      mutate(carrier_factor = factor(.data[[variant]], levels = c(1, 0), labels = c("Carrier", "Not_Carrier"))) %>%
      select(case_factor, carrier_factor) %>%
      table() %>%
      fisher.test(alternative = "g") ->
      # fisher.test() ->
      fisher_results[[variant]]
  }
  fisher_tibble = NULL
  for (i in seq_along(fisher_results)) {
    temp = broom::tidy(fisher_results[[i]])
    temp$gene = names(fisher_results)[i]
    fisher_tibble = rbind(fisher_tibble, temp)
  }
  fisher_tibble %>%
    mutate(p.adjust = p.adjust(p.value, method = "fdr")) ->
    fisher_tibble
  fisher_tibble
  write_tsv(fisher_tibble, paste0(results, paste(tsd_proj, "fisher.tsv", sep = "_")))


  # survival analysis
  phs_wgs_fit = coxph(Surv(age, status) ~ phs, wgs)
  summary(phs_wgs_fit)
  phs_wgs_coxph = broom::tidy(phs_wgs_fit, exponentiate = TRUE, conf.int = TRUE)
  cox.zph(phs_wgs_fit)

  hoxb13_wgs_fit = coxph(Surv(age, status) ~ hoxb13_carrier, wgs)
  summary(hoxb13_wgs_fit)
  hoxb13_wgs_coxph = broom::tidy(hoxb13_wgs_fit, exponentiate = TRUE, conf.int = TRUE)
  cox.zph(hoxb13_wgs_fit)
  # hoxb13_wgs_fit = coxph(Surv(age_of_metastasis, metastasis) ~ hoxb13_carrier, wgs)
  # summary(hoxb13_wgs_fit)

  phs_hoxb13_wgs_fit = coxph(Surv(age, status) ~ phs + hoxb13_carrier, wgs)
  summary(phs_hoxb13_wgs_fit)
  phs_hoxb13_wgs_coxph = broom::tidy(phs_hoxb13_wgs_fit, exponentiate = TRUE, conf.int = TRUE)

  cox.zph_test <- cox.zph(phs_hoxb13_wgs_fit)
  plot(cox.zph_test)

  phs_hoxb13_wgs_inter = coxph(Surv(age, status) ~ phs * hoxb13_carrier, wgs)
  summary(phs_hoxb13_wgs_inter)
  phs_hoxb13_inter_coxph = broom::tidy(phs_hoxb13_wgs_inter, exponentiate = TRUE, conf.int = TRUE)


  # phs_hoxb13_snp_fit = coxph(Surv(age, status) ~ phs + hoxb13_snp, carriers)
  # summary(phs_hoxb13_snp_fit)
  # hoxb13_snp_fit = coxph(Surv(age, status) ~ hoxb13_snp, carriers)
  # summary(hoxb13_snp_fit)

  ## New performance indices for PHS 

  # Performance metrics
  HR80_20 = get_hr(wgs,
                  CI = TRUE, 
                  boot = bootstrap_iters)
  HR80_60 = get_hr(wgs,
                  lower_interval = c(0.2, 0.80),
                  upper_interval = c(0.8, 1.0),
                  CI = TRUE, 
                  boot = bootstrap_iters)

  ages = c(60, 70, 80)
  OR80_20 = lapply(ages, function(x) {get_or(wgs,
                                            or_age = x, 
                                            CI = TRUE, 
                                            boot = bootstrap_iters)})
  c_index = get_cindex(wgs,
                      CI = TRUE, 
                      boot = bootstrap_iters)

  HR80_20$age = NA
  HR80_60$age = NA
  HR80_60$index = "HR80_60"
  c_index$age = NA
  indices = list(HR80_20, HR80_60, OR80_20[[1]], OR80_20[[2]], OR80_20[[3]], c_index)

  perform = do.call(rbind.data.frame, lapply(indices, function(x){data.frame(index = x$index, 
                                                                  value = x$value, 
                                                                  conf.low = x$conf.low, 
                                                                  conf.high = x$conf.high, 
                                                                  age = x$age)}))
  perform = cbind(model, perform)
  if (perform_bool) {
    write_tsv(perform, paste0(results, paste(tsd_proj, model, "wgs_phs_performance.tsv", sep = "_")))
  }

  ## New plot of PHS km curves but on the subsetted WGS dataset

  # KM Curves
  curves = data.frame(curve0_20 = c(0,   0.2),
                      curve20_80 = c(0.2, 0.8),
                      curve80_100 = c(0.8, 1))


  km_curves = data.frame()
  for (i in seq_len(length(curves))) {
    curven <- km_curve(data = wgs,
                      interval = curves[[i]],
                      age_range = 40:100,
                      scale = FALSE, 
                      inverse = FALSE)
    curven$label = label_generator(curves[1, i], curves[2, i])
    km_curves = rbind(km_curves, curven)
  }


  wgs_km = ggplot(km_curves, aes(x = time,
                        y = estimate,
                        ymin = conf.low,
                        ymax = conf.high,
                        col = label,
                        fill = label)) +
    geom_ribbon(alpha = 0.1,
                color = 0) +
    geom_step() +
    theme_classic() +
    xlim(40, 100) + 
    ylim(0, 1) +
    labs(x = "Age", y = "Disease-free Survival") +
    scale_color_brewer(palette = "Set1", direction = -1,
                      name = "PHS Centile") +
    scale_fill_brewer(palette = "Set1", direction = -1,
                      name = "PHS Centile") + 
    theme(aspect.ratio = 1) +
        theme(legend.position = "bottom") 

  # ggsave(paste0(results, paste(tsd_proj, model, "phs_wgs_km_curve.png", sep = "_")), wgs_km,
  #        height = 6 * 0.618, width = 6)

  saveRDS(km_curve, paste0(results, paste(tsd_proj, model, "phs_wgs_km_curve.rds", sep = "_")))

  ## Plot of HOXB13 km curve on WGS data
  # KM Curves
    km_curves = data.frame()
    noncarrier_data <- wgs[wgs$hoxb13_carrier == 0, ]
    carrier_data <- wgs[wgs$hoxb13_carrier == 1, ]
    mod <- survfit(Surv(age, status) ~ 1, data = noncarrier_data)
    mod_data <- data.frame(time = mod$time, n.risk = mod$n.risk, 
        n.event = mod$n.event, n.censor = mod$n.censor, estimate = mod$surv, 
        std.error = mod$std.err, conf.high = mod$upper, conf.low = mod$lower, 
        cumhaz = mod$cumhaz)
    curven = select(mod_data, time, estimate, conf.low, 
        conf.high, cumhaz)
    curven$label = "noncarrier"
    km_curves = rbind(km_curves, curven)
    mod <- survfit(Surv(age, status) ~ 1, data = carrier_data)
    mod_data <- data.frame(time = mod$time, n.risk = mod$n.risk, 
        n.event = mod$n.event, n.censor = mod$n.censor, estimate = mod$surv, 
        std.error = mod$std.err, conf.high = mod$upper, conf.low = mod$lower, 
        cumhaz = mod$cumhaz)
    curven = select(mod_data, time, estimate, conf.low, 
        conf.high, cumhaz)
    curven$label = "carrier"
    km_curves = rbind(km_curves, curven)
    
  hoxb13_p = ggplot(km_curves, aes(x = time,
                        y = estimate,
                        ymin = conf.low,
                        ymax = conf.high,
                        col = label,
                        fill = label)) +
    geom_ribbon(alpha = 0.1,
                color = 0) +
    geom_step() +
    theme_classic() +
    xlim(40, 100) + 
    ylim(0, 1) +
    labs(x = "Age", y = "Disease-free Survival") +
    scale_color_brewer(palette = "Set1", direction = -1,
                      name = "PHS Centile") +
    scale_fill_brewer(palette = "Set1", direction = -1,
                      name = "PHS Centile") + 
    theme(aspect.ratio = 1) +
        theme(legend.position = "bottom") 
    
  # fit = coxph(Surv(age, status) ~ hoxb13_carrier, data = wgs)
  # #fit = coxph(Surv(age, status) ~ hoxb13_snp, data = carriers)
  # summary(fit)
  # new_data0 = data.frame(hoxb13_carrier = c(0))
  # new_data1 = data.frame(hoxb13_carrier = c(1))
  # fit_surv0 = survfit(fit, newdata = new_data0)
  # df0 = data.frame(hoxb13_carrier = 0, 
  #                 surv = fit_surv0$surv, 
  #                 time = fit_surv0$time,
  #                 lower = fit_surv0$lower, 
  #                 upper = fit_surv0$upper)
  # fit_surv1 = survfit(fit, newdata = new_data1)
  # df1 = data.frame(hoxb13_carrier = 1, 
  #                 surv = fit_surv1$surv, 
  #                 time = fit_surv1$time,
  #                 lower = fit_surv1$lower, 
  #                 upper = fit_surv1$upper)
  # prediction = rbind(df0, df1) %>%
  #     mutate(hoxb13_carrier = as.factor(hoxb13_carrier))

  # hoxb13_p = ggplot(prediction, aes(x = time, 
  #                       y = surv, 
  #                       ymin = lower,
  #                       ymax = upper,
  #                       color = hoxb13_carrier,
  #                       fill = hoxb13_carrier)) +
  #     geom_ribbon(alpha = 0.1,
  #                 color = 0) +
  #     geom_step() +
  #     theme_classic() +
  #     xlim(40, 100) + 
  #     ylim(0, 1) +
  #     labs(x = "Age", y = "Disease-free Survival") +
  #     scale_color_brewer(palette = "Set1",
  #                       name = "HOXB13 Carrier") +
  #     scale_fill_brewer(palette = "Set1",
  #                       name = "HOXB13 Carrier") + 
  #   theme(aspect.ratio = 1) +
  #       theme(legend.position = "bottom") 

  wgs_p = wgs_km +
      hoxb13_p + 
      plot_annotation(tag_levels = 'A')
  ggsave(paste0(results, paste(tsd_proj, model, "wgs_km_curve.pdf", sep = "_")), wgs_p,
            height = 8 * 0.618, width = 8)
  ggsave(paste0(results, paste(tsd_proj, model, "wgs_km_curve.png", sep = "_")), wgs_p,
            height = 8 * 0.618, width = 8)
  saveRDS(wgs_p, paste0(results, paste(tsd_proj, model, "hoxb13_km_curve.rds", sep = "_")))

  #####################################################################################
  # HOXB13 SNPs
  #####################################################################################


  # survival analysis
  phs_fit = coxph(Surv(age, status) ~ phs, combined_data)
  summary(phs_fit)


  hoxb13_snp_fit = coxph(Surv(age, status) ~ hoxb13_snp, combined_data)
  summary(hoxb13_snp_fit)
  hoxb13_snp_coxph = broom::tidy(hoxb13_snp_fit, exponentiate = TRUE, conf.int = TRUE)
  cox.zph(hoxb13_snp_fit)
  # hoxb13_snp_fit = coxph(Surv(age_of_metastasis, metastasis) ~ hoxb13_snp, combined_data)
  # summary(hoxb13_snp_fit)

  phs_hoxb13_snp_fit = coxph(Surv(age, status) ~ phs + hoxb13_snp, combined_data)
  summary(phs_hoxb13_snp_fit)
  phs_hoxb13_snp_coxph = broom::tidy(phs_hoxb13_snp_fit, exponentiate = TRUE, conf.int = TRUE)

  cox.zph_test <- cox.zph(phs_hoxb13_snp_fit)
  print(cox.zph_test)
  plot(cox.zph_test)

  phs_hoxb13_snp_inter = coxph(Surv(age, status) ~ phs * hoxb13_snp, combined_data)
  summary(phs_hoxb13_wgs_inter)
  phs_hoxb13_snp_inter_coxph = broom::tidy(phs_hoxb13_snp_inter, exponentiate = TRUE, conf.int = TRUE)

  ## Plot of HOXB13 SNPs km curve on entire data

  fit = coxph(Surv(age, status) ~ hoxb13_snp, data = combined_data)
  summary(fit)
  new_data0 = data.frame(hoxb13_snp = c(0))
  new_data1 = data.frame(hoxb13_snp = c(1))
  fit_surv0 = survfit(fit, newdata = new_data0)
  df0 = data.frame(hoxb13_snp = 0, 
                  surv = fit_surv0$surv, 
                  time = fit_surv0$time,
                  lower = fit_surv0$lower, 
                  upper = fit_surv0$upper)
  fit_surv1 = survfit(fit, newdata = new_data1)
  df1 = data.frame(hoxb13_snp = 1, 
                  surv = fit_surv1$surv, 
                  time = fit_surv1$time,
                  lower = fit_surv1$lower, 
                  upper = fit_surv1$upper)
  prediction = rbind(df0, df1) %>%
      mutate(hoxb13_snp = as.factor(hoxb13_snp))

  hoxb13_snp_p = ggplot(prediction, aes(x = time, 
                        y = surv, 
                        ymin = lower,
                        ymax = upper,
                        color = hoxb13_snp,
                        fill = hoxb13_snp)) +
      geom_ribbon(alpha = 0.1,
                  color = 0) +
      geom_step() +
      theme_classic() +
      xlim(40, 100) + 
      ylim(0, 1) +
      labs(x = "Age", y = "Disease-free Survival") +
      scale_color_brewer(palette = "Set1",
                        name = "HOXB13 Carrier") +
      scale_fill_brewer(palette = "Set1",
                        name = "HOXB13 Carrier") + 
    theme(aspect.ratio = 1) +
        theme(legend.position = "bottom") 

  fig_width = 6
  ggsave(paste0(results, paste(tsd_proj, model, "hoxb13_snp_km_curve.pdf", sep = "_")), hoxb13_snp_p,
            height = fig_width * 0.618, width = fig_width)
  ggsave(paste0(results, paste(tsd_proj, model, "hoxb13_snp_km_curve.png", sep = "_")), hoxb13_snp_p,
            height = fig_width * 0.618, width = fig_width)
  saveRDS(hoxb13_snp_p, paste0(results, paste(tsd_proj, model, "hoxb13_snp_km_curve.rds", sep = "_")))



  #####################################################################################
  # # Conversion to metastasis
  #####################################################################################

  # combined_data %>%
  #   select(conversion_to_metastasis_years, metastasis, age_lastobs, age_diagnosis) %>%
  #   print(n=100, width = 1000)

  # combined_data %>%
  #   filter(status == 1) %>%
  #   mutate(conversion_to_metastasis_years = if_else(metastasis == 1, 
  #       conversion_to_metastasis_years, age_lastobs - age_diagnosis)) %>%
  #   select(phs, conversion_to_metastasis_years, metastasis) ->
  #   conversion

  # fit = coxph(Surv(conversion_to_metastasis_years, metastasis) ~ phs, data = conversion)
  # summary(fit)


  # Does having a high PHS icnrease your risk of carrying a rare variant?
  binomial_phs_hoxb13 <- glm(hoxb13_carrier ~ phs, family="binomial", data=wgs)
  summary(binomial_phs_hoxb13)
  broom::tidy(binomial_phs_hoxb13, conf.int = TRUE, exponentiate = TRUE) %>%
    write_tsv(paste0(results, paste(tsd_proj, model, "binomial_phs_hoxb13.tsv", sep = "_")))

  phs_wgs_coxph$model = "phs_wgs"
  hoxb13_wgs_coxph$model = "hoxb13_wgs"
  phs_hoxb13_wgs_coxph$model = "phs_hoxb13"
  phs_hoxb13_inter_coxph$model = "phs_hoxb13_inter"
  hoxb13_snp_coxph$model = "hoxb13_snp"

  rbind(phs_wgs_coxph, 
  hoxb13_wgs_coxph, 
  phs_hoxb13_wgs_coxph, 
  phs_hoxb13_inter_coxph,
  hoxb13_snp_coxph) ->
 wgs_results
  write_tsv(wgs_results, paste0(results, paste(tsd_proj, model, "wgs.tsv", sep = "_")))

}

# Combine Coxph results
any_pca_coxph$model = "any_pca"
metastatic_coxph$model = "metastatic"
phs_tt_tidy$model = "phs_tt"

rbind(any_pca_coxph, 
metastatic_coxph, 
phs_tt_tidy) ->
coxph_results
write_tsv(coxph_results, paste0(results, paste(tsd_proj, model, "coxph.tsv", sep = "_")))



