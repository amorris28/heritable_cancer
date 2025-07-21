library(tidyverse)
library(stringr)
library(hazrd)
library(survival)
source("src/00_config.R")

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
inverse = args[3]
inverse = FALSE
tsd_proj = args[4]
tsd_proj = "p212"
output = "reports/"

figs = "reports/figures/"

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

#########################################
# Analyze PHS data

# PHS Histogram
phs_hist = 
	phs_hist(combined_data, normalize = TRUE) +
	labs(title = model) + theme(aspect.ratio=1)

ggsave(paste0(figs, paste(tsd_proj, model, "hist.png", sep = "_")), phs_hist, 
       units = "px", height = 1300, width = 1300)

# Performance metrics


HR80_20 = get_hr(combined_data, CI = TRUE, boot = 300)
HR80_60 = get_hr(combined_data, 
                 lower_interval = c(0.2, 0.80),
                 upper_interval = c(0.8, 1.0),
                 CI = TRUE, 
                 boot = 300)

ages = c(60, 70, 80)
OR80_20 = lapply(ages, function(x) {get_or(combined_data, or_age = x, CI = TRUE, boot = 300)})
c_index = get_cindex(combined_data, CI = TRUE, boot = 300)

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
write_tsv(perform, paste0(output, paste(tsd_proj, model, "performance.tsv", sep = "_")))


# KM Curves
curves = data.frame(curve0_20 = c(0,   0.2),
                    curve20_80 = c(0.2, 0.8),
                    curve80_100 = c(0.8, 1.0))

label_generator = function(x, y) {
  x = x * 100
  y = y * 100
  out = paste0("PHS ", x, "-", y, "th centile")
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
  labs(x = "Age", y = "Disease-free Survival", title = model) +
  scale_color_brewer(palette = "Set1", direction = -1,
                     name = "Centile") +
  scale_fill_brewer(palette = "Set1", direction = -1,
                    name = "Centile")  + theme(aspect.ratio=1)

ggsave(paste0(figs, paste(tsd_proj, model, "km_curve.png", sep = "_")), km_curve,
       units = "px", height = 1300, width = 1600)


# BRCA2, HOXB13
combined_data %>%
  select(status, age, phs, wgs, brca2_carrier:hoxb13_snp) %>%
  mutate(across(brca2_carrier:atm_carrier, ~ if_else(is.na(.x), 0, .x))) ->
  carriers

carriers %>%
  filter(wgs == 1) %>%
  mutate(any_carrier = case_when(if_any(brca2_carrier:atm_carrier, ~ .x == 1) ~ 1, .default = 0)) %>%
  mutate(all_carrier = rowSums(across(brca2_carrier:atm_carrier))) ->
  wgs
wgs %>%
  select(-msh2_carrier, -msh6_carrier) ->
  wgs
  # 0 subjects with these mutations

 #Fisher test
variants = c("brca2", "hoxb13", "brca1", "mlh1", "atm", "any", "all")
fisher_results = list()
for (i in 1:length(variants)) {
  variant = paste0(variants[i], "_carrier")
  print(variant)
  wgs %>%
    select(status, all_of(variant)) %>%
    table() %>%
    fisher.test() ->
    fisher_results[[variant]]
}

# survival analysis

phs_hoxb13_wgs_fit = coxph(Surv(age, status) ~ phs + hoxb13_carrier, wgs)
summary(phs_hoxb13_wgs_fit)
hoxb13_wgs_fit = coxph(Surv(age, status) ~ hoxb13_carrier, wgs)
summary(hoxb13_wgs_fit)

phs_hoxb13_snp_fit = coxph(Surv(age, status) ~ phs + hoxb13_snp, carriers)
summary(phs_hoxb13_snp_fit)
hoxb13_snp_fit = coxph(Surv(age, status) ~ hoxb13_snp, carriers)
summary(hoxb13_snp_fit)


## Write a function that returns a pretty ggplot with curves separated by covariate


fit = coxph(Surv(age, status) ~ hoxb13_carrier, data = wgs)
#fit = coxph(Surv(age, status) ~ hoxb13_snp, data = carriers)
summary(fit)
new_data0 = data.frame(hoxb13_carrier = c(0))
new_data1 = data.frame(hoxb13_carrier = c(1))
fit_surv0 = survfit(fit, newdata = new_data0)
df0 = data.frame(hoxb13_carrier = 0, 
                 surv = fit_surv0$surv, 
                 time = fit_surv0$time,
                 lower = fit_surv0$lower, 
                 upper = fit_surv0$upper)
fit_surv1 = survfit(fit, newdata = new_data1)
df1 = data.frame(hoxb13_carrier = 1, 
                 surv = fit_surv1$surv, 
                 time = fit_surv1$time,
                 lower = fit_surv1$lower, 
                 upper = fit_surv1$upper)
prediction = rbind(df0, df1) %>%
    mutate(hoxb13_carrier = as.factor(hoxb13_carrier))

p = ggplot(prediction, aes(x = time, 
                       y = surv, 
                       ymin = lower,
                       ymax = upper,
                       color = hoxb13_carrier,
                       fill = hoxb13_carrier)) +
    geom_ribbon(alpha = 0.1,
                color = 0) +
    geom_step() +
    theme_classic() +
    xlim(40, 106) + 
    ylim(0, 1) +
    labs(x = "Age", y = "Disease-free Survival") +
    scale_color_brewer(palette = "Set1",
                       name = "Hoxb13 Carrier") +
    scale_fill_brewer(palette = "Set1",
                      name = "Hoxb13 Carrier")

ggsave(paste0(figs, paste(tsd_proj, model, "hoxb13_km_curve.png", sep = "_")), p,
       units = "px", height = 1300, width = 1600)


# Time to metastatic diseaes

combined_data %>%
  mutate(metastasis = if_else(is.na(metastasis), 0, metastasis)) %>%
  mutate(age_of_metastasis = if_else(is.na(age_of_metastasis), age_lastobs, age_of_metastasis)) ->
  combined_data
metastasis_fit = coxph(Surv(age_of_metastasis, metastasis) ~ phs, combined_data)
summary(metastasis_fit)

combined_data %>%
  select(age_lastobs, age_diagnosis, age, age_of_metastasis, age_at_death) %>%
  rowwise() %>%
  mutate(max = max(c_across(age_lastobs:age_at_death), na.rm = T)) ->
  ages



# KM Curves for metastatic disease
curves = data.frame(curve0_20 = c(0,   0.2),
                    curve20_80 = c(0.2, 0.8),
                    curve80_100 = c(0.8, 1.0))

label_generator = function(x, y) {
  x = x * 100
  y = y * 100
  out = paste0("PHS ", x, "-", y, "th centile")
  return(out)
}

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
  ylim(0, 1) +
  labs(x = "Age", y = "Disease-free Survival", title = model) +
  scale_color_brewer(palette = "Set1", direction = -1,
                     name = "Centile") +
  scale_fill_brewer(palette = "Set1", direction = -1,
                    name = "Centile")  + theme(aspect.ratio=1)

ggsave(paste0(figs, paste(tsd_proj, model, "km_curve_metastatic.png", sep = "_")), km_curve,
       units = "px", height = 1300, width = 1600)
