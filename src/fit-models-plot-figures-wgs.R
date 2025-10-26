library(tidyverse)
library(stringr)
library(hazrd)
library(survival)
library(patchwork)
source("src/config.R")
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

args = commandArgs(trailingOnly=TRUE)

model_file = args[1] # PHS file
metadata_file = args[2] # Phenotype data
inverse = args[3]
perform_bool = args[4]


#if (Sys.info()["nodename"] == "p2521-hpc-01.tsd.usit.no") {
#    metadata_file = "data/processed/husk_phenotypes.tsv"
#    tsd_proj = "p2521"
#} else if (Sys.info()["nodename"] == "p212-hpc-01.tsd.usit.no") {
#    metadata_file = "data/processed/heritable_cancer_phenotypes.tsv" # Phenotype data
#    tsd_proj = "p212"
#}
model = tools::file_path_sans_ext(basename(model_file))

#model_file = "data/processed/PHS601_phs" # PHS file
#model = gsub("(.+?)(\\_.*)", "\\1", basename(model_file))
#bootstrap_iters = 999
#########################################
# Import PHS data and metadata
phs_data = import_phs_data(model_file, inverse)

metadata = read_tsv(metadata_file) # Phenotype data

combined_data = process_combined_data(metadata, phs_data)

combined_data %>%
    select(status, age, phs, age_of_metastasis, metastasis, wgs, brca2_carrier:atm_carrier) %>%
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

write_tsv(wgs_summary, paste(model, "wgs_summary.tsv", sep = "_"))

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
write_tsv(fisher_tibble, paste(model, "fisher.tsv", sep = "_"))

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

cox.zph(phs_hoxb13_wgs_fit)

phs_hoxb13_wgs_inter = coxph(Surv(age, status) ~ phs * hoxb13_carrier, wgs)
summary(phs_hoxb13_wgs_inter)
phs_hoxb13_inter_coxph = broom::tidy(phs_hoxb13_wgs_inter, exponentiate = TRUE, conf.int = TRUE)

cox.zph(phs_hoxb13_wgs_inter)

# phs_hoxb13_snp_fit = coxph(Surv(age, status) ~ phs + hoxb13_snp, carriers)
# summary(phs_hoxb13_snp_fit)
# hoxb13_snp_fit = coxph(Surv(age, status) ~ hoxb13_snp, carriers)
# summary(hoxb13_snp_fit)

#################################################################################
# Performance metrics, PHS

if (perform_bool) {
    perform = get_performance_metrics(wgs, bootstrap_iters, model)
    write_tsv(perform, paste(model, "wgs_phs_performance.tsv", sep = "_"))
}


sum(wgs$hoxb13_carrier) / nrow(wgs) # 1.7%
HR95_50 = get_hr(wgs, 
                lower_interval = c(0.2, 0.80), # median
                upper_interval = c(1.00-.017, 1),  # top 1.7%
                CI = TRUE, 
                bootstrap_iterations = bootstrap_iters)
c_index = get_cindex(wgs,
                  CI = TRUE, 
                  bootstrap_iterations = bootstrap_iters)
HR95_50$index = "HR"

indices = list(HR95_50, c_index)

perform = do.call(rbind.data.frame, lapply(indices, function(x){data.frame(index = x$index, 
                                                                value = x$value, 
                                                                conf.low = x$conf.low, 
                                                                conf.high = x$conf.high)}))
perform = cbind(model, perform)
write_tsv(perform, paste(model, "wgs_phs_performance_95.tsv", sep = "_"))

#################################################################################
# Performance metrics, HOXB13
# HR, C Index
if (perform_bool) {
    get_hr_hoxb13 = function(wgs, bootstrap_iterations = 999) {
        return_hr = function(wgs) {summary(coxph(Surv(age, status) ~ hoxb13_carrier, data = wgs))$coefficients[2]}
        HR = return_hr(wgs)
        boot_out = boot_conf(wgs, bootstrap_iters, return_hr)
        return(list(index = "HR", 
            value = HR, 
            conf.low = boot_out$quantiles[[1]], 
            conf.high = boot_out$quantiles[[2]], 
            iters = boot_out$iters))
    }
    return_cindex = function(wgs) {as.numeric(coxph(Surv(age, status) ~ hoxb13_carrier, data = wgs)$concordance['concordance'])}
    get_cindex_hoxb13 = function(wgs, bootstrap_iterations = 999) {
        c_index = return_cindex(wgs)
        boot_out = boot_conf(wgs, bootstrap_iters, return_cindex)
        return(list(index = "C_Index", 
            value = c_index, 
            conf.low = boot_out$quantiles[[1]], 
            conf.high = boot_out$quantiles[[2]], 
            iters = boot_out$iters))
    }
    get_performance_metrics_hoxb13 = function(wgs, bootstrap_iters, model) {

    HR = get_hr_hoxb13(wgs, bootstrap_iters)

    c_index = get_cindex_hoxb13(wgs,
                                bootstrap_iterations = bootstrap_iters)

    indices = list(HR, c_index)

    perform = do.call(rbind.data.frame, lapply(indices, function(x){data.frame(index = x$index, 
                                                                    value = x$value, 
                                                                    conf.low = x$conf.low, 
                                                                    conf.high = x$conf.high)}))
    perform = cbind(model, perform)

    return(perform)
    }
    perform = get_performance_metrics_hoxb13(wgs, boostrap_iters, model)
    write_tsv(perform, paste(model, "wgs_hoxb13_performance.tsv", sep = "_"))
}
## New plot of PHS km curves but on the subsetted WGS dataset

# KM Curves
wgs_km = plot_km_curves(wgs)

# ggsave(paste0(results, paste(tsd_proj, model, "phs_wgs_km_curve.png", sep = "_")), wgs_km,
#        height = 6 * 0.618, width = 6)

#saveRDS(km_curve, paste0(results, paste(tsd_proj, model, "phs_wgs_km_curve.rds", sep = "_")))

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

color_set = RColorBrewer::brewer.pal(9, "Set1")
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
    scale_color_manual(values = c(color_set[5], color_set[2]),
                        name = "HOXB13 Carrier") +
    scale_fill_manual(values = c(color_set[5], color_set[2]),
                        name = "HOXB13 Carrier") + 
    # scale_color_brewer(palette = "Set1", direction = -1,
    #                     name = "HOXB13 Carrier") +
    # scale_fill_brewer(palette = "Set1", direction = -1,
    #                     name = "HOXB13 Carrier") + 
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

wgs_p = 
    wgs_km +
    hoxb13_p + 
    plot_annotation(tag_levels = 'A') &
    theme(legend.position = "bottom") 
wgs_p
ggsave(paste(model, "wgs_km_curve.pdf", sep = "_"), wgs_p,
        height = 8 * 0.618, width = 8)
ggsave(paste(model, "wgs_km_curve.png", sep = "_"), wgs_p,
        height = 8 * 0.618, width = 8)
saveRDS(wgs_p, paste(model, "wgs_km_curve.rds", sep = "_"))

#####################################################################################
# HOXB13 SNPs
#####################################################################################


## survival analysis
#phs_fit = coxph(Surv(age, status) ~ phs, combined_data)
#summary(phs_fit)
#
#
#hoxb13_snp_fit = coxph(Surv(age, status) ~ hoxb13_snp, combined_data)
#summary(hoxb13_snp_fit)
#hoxb13_snp_coxph = broom::tidy(hoxb13_snp_fit, exponentiate = TRUE, conf.int = TRUE)
#cox.zph(hoxb13_snp_fit)
## hoxb13_snp_fit = coxph(Surv(age_of_metastasis, metastasis) ~ hoxb13_snp, combined_data)
## summary(hoxb13_snp_fit)
#
#phs_hoxb13_snp_fit = coxph(Surv(age, status) ~ phs + hoxb13_snp, combined_data)
#summary(phs_hoxb13_snp_fit)
#phs_hoxb13_snp_coxph = broom::tidy(phs_hoxb13_snp_fit, exponentiate = TRUE, conf.int = TRUE)
#
#cox.zph_test <- cox.zph(phs_hoxb13_snp_fit)
#print(cox.zph_test)
#plot(cox.zph_test)
#
#phs_hoxb13_snp_inter = coxph(Surv(age, status) ~ phs * hoxb13_snp, combined_data)
#summary(phs_hoxb13_wgs_inter)
#phs_hoxb13_snp_inter_coxph = broom::tidy(phs_hoxb13_snp_inter, exponentiate = TRUE, conf.int = TRUE)
#
### Plot of HOXB13 SNPs km curve on entire data
#
#fit = coxph(Surv(age, status) ~ hoxb13_snp, data = combined_data)
#summary(fit)
#new_data0 = data.frame(hoxb13_snp = c(0))
#new_data1 = data.frame(hoxb13_snp = c(1))
#fit_surv0 = survfit(fit, newdata = new_data0)
#df0 = data.frame(hoxb13_snp = 0, 
#                surv = fit_surv0$surv, 
#                time = fit_surv0$time,
#                lower = fit_surv0$lower, 
#                upper = fit_surv0$upper)
#fit_surv1 = survfit(fit, newdata = new_data1)
#df1 = data.frame(hoxb13_snp = 1, 
#                surv = fit_surv1$surv, 
#                time = fit_surv1$time,
#                lower = fit_surv1$lower, 
#                upper = fit_surv1$upper)
#prediction = rbind(df0, df1) %>%
#    mutate(hoxb13_snp = as.factor(hoxb13_snp))
#
#hoxb13_snp_p = ggplot(prediction, aes(x = time, 
#                    y = surv, 
#                    ymin = lower,
#                    ymax = upper,
#                    color = hoxb13_snp,
#                    fill = hoxb13_snp)) +
#    geom_ribbon(alpha = 0.1,
#                color = 0) +
#    geom_step() +
#    theme_classic() +
#    xlim(40, 100) + 
#    ylim(0, 1) +
#    labs(x = "Age", y = "Disease-free Survival") +
#    scale_color_brewer(palette = "Set1",
#                    name = "HOXB13 Carrier") +
#    scale_fill_brewer(palette = "Set1",
#                    name = "HOXB13 Carrier") + 
#theme(aspect.ratio = 1) +
#    theme(legend.position = "bottom") 
#
#fig_width = 6
#ggsave(paste0(results, paste(tsd_proj, model, "hoxb13_snp_km_curve.pdf", sep = "_")), hoxb13_snp_p,
#        height = fig_width * 0.618, width = fig_width)
#ggsave(paste0(results, paste(tsd_proj, model, "hoxb13_snp_km_curve.png", sep = "_")), hoxb13_snp_p,
#        height = fig_width * 0.618, width = fig_width)
#saveRDS(hoxb13_snp_p, paste0(results, paste(tsd_proj, model, "hoxb13_snp_km_curve.rds", sep = "_")))



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

#####################################################################################
# Does having a high PHS icnrease your risk of carrying a rare variant?
#####################################################################################
binomial_phs_hoxb13 <- glm(hoxb13_carrier ~ phs, family="binomial", data=wgs)
summary(binomial_phs_hoxb13)
broom::tidy(binomial_phs_hoxb13, conf.int = TRUE, exponentiate = TRUE) %>%
write_tsv(paste(model, "binomial_phs_hoxb13.tsv", sep = "_"))

phs_wgs_coxph$model = "phs_wgs"
hoxb13_wgs_coxph$model = "hoxb13_wgs"
phs_hoxb13_wgs_coxph$model = "phs_hoxb13"
phs_hoxb13_inter_coxph$model = "phs_hoxb13_inter"
#hoxb13_snp_coxph$model = "hoxb13_snp"

rbind(phs_wgs_coxph, 
hoxb13_wgs_coxph, 
phs_hoxb13_wgs_coxph, 
phs_hoxb13_inter_coxph) ->
#hoxb13_snp_coxph) ->
wgs_results
write_tsv(wgs_results, paste(model, "wgs.tsv", sep = "_"))
