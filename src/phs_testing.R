library(tidyverse)
library(survival)
#library(survminer)
#library(broom)

source("src/00_config.R")


#########################################
# First argument: PHS file
# tab delimited file, no header, first column ID, second column PHS score
# Second argument: metadata file
# tab delimited file with column names bfile_id, age, status. Additional
# columns will be ignored.
# Third argument: 
# boolean indicating whether to inverse (x * -1) the PHS scores
# to reverse the direction of effect.

args = commandArgs(trailingOnly=TRUE)

model_file = args[1] # PHS file

model = basename(args[1]) # filename for output

metadata = read_tsv(args[2]) # Phenotype data

inverse = args[3]

output = "reports/"

figs = "reports/figures/"

#########################################
# Import PHS data

phs = read_tsv(model_file, col_names = c("bfile_id", "phs"))

phs$phs = scale(phs$phs, center = TRUE, scale = TRUE)

if (inverse) {phs$phs = phs$phs * -1}

metadata %>%
  left_join(phs, by = "bfile_id") %>%
  filter(!(is.na(phs))) ->
  combined_data

phs_hist = ggplot(combined_data, aes(x = phs, group = status, fill = as.character(status))) +
  geom_histogram() + 
  scale_fill_manual(values = c("#132B43", "#56B1F7")) +
  theme_minimal()

ggsave(paste0(figs, model, "_hist.png"), phs_hist)

# Fit model

agevals = 40:100
reference_PHS = combined_data$phs
prctile_list = c(0, 0.2, 0.8, 1)
prctile_labels = c('0-20th', '20-80th', '80-100th')
combined_data$percentile = 
  with(combined_data,
       cut(phs, breaks = quantile(phs, probs = prctile_list),
                                      labels = prctile_labels,
                                      include.lowest = TRUE))

cox_model = coxph(Surv(age, status) ~ phs, data = combined_data)

bt = cox_model$coefficients[["phs"]]

combined_data %>% 
  filter(percentile %in% c('0-20th', '20-80th', '80-100th')) ->
  combined_data

combined_data %>% 
  group_by(percentile) %>% 
  summarize(median_phs = median(phs)) ->
  median_tibble


median_phs_individual = quantile(reference_PHS, 0.50)
critvalvec = median_tibble$median_phs

#########################
# Calculate idealized curves
#
# Basically copied this from:
# https://github.com/cmig-research-group/phs/blob/master/phs_model_example/phs_generate_new_demo.m
a = 0.0700 / 100;
b = 0.0753
pcainc = a * exp(b * (agevals - 40)) # this is the population incidence rate for PCa
H0 = a / b * exp(b * (agevals - 40)) # compute baseline hazard
H0 = H0 - H0[1] # assume cumulative hazard == 0 before age 40

H = matrix(nrow = length(H0), ncol = length(critvalvec))
S = matrix(nrow = length(H0), ncol = length(critvalvec))
for (cnt in 1:length(critvalvec)) {
  H[, cnt] = H0 * exp(bt * critvalvec[cnt] - bt * median_phs_individual[[1]])
  S[, cnt] = exp(-H[, cnt])
}

S_tibble = as_tibble(S, .name_repair = ~ as.character(median_tibble$percentile))
S_tibble$agevals = agevals
S_tibble %>%
  pivot_longer(-agevals, names_to = "percentile") ->
  S_tibble

##########################

# Plot K-M curves for centiles
mod = survfit(Surv(age, status) ~ percentile, data = combined_data)

tidy_strata = function(strata) {
  temp = c()
  for (i in seq(strata)) {
    temp = c(temp, rep_len(names(strata)[i], strata[i]))
  }
  return(temp)
}

mod_data = tibble(time = mod$time,
           n.risk = mod$n.risk,
           n.event = mod$n.event,
           n.censor = mod$n.censor,
           estimate = mod$surv,
           std.error = mod$std.err,
           conf.high = mod$upper,
           conf.low = mod$lower,
           strata = tidy_strata(mod$strata))
#mod_data = tidy(mod)

tibble(strata = unique(mod_data$strata),
       percentile = unique(S_tibble$percentile)) ->
  key

mod_data %>% 
  left_join(key, by = 'strata') ->
  mod_data

kmcurve = ggplot(mod_data, aes(x = time, y = estimate, color = percentile)) +
  geom_step() +
  theme_classic() +
  xlim(40, 100) + 
#  geom_line(aes(x = agevals, y = value, color = percentile), S_tibble) + 
  scale_color_brewer(palette = "Set1") +
  labs(x = "Age", y = "PCa-free Survival")

ggsave(paste0(figs, model, "_km_curve.png"),
       kmcurve)

combined_data %>% 
  group_by(status) %>% 
  summarize(n = n()) ->
  cases_controls

## Hazard ratios
source("src/RK_get_perf.R")
perform = with(combined_data, 
               RK_get_perf(phs, 
                           age, 
                           status, 
                           ref = phs, 
                           swc.popnumcases = cases_controls[[2, 2]], 
                           swc.popnumcontrols = cases_controls[[1, 2]]))
perform = stack(perform)[, c(2, 1)]
write_tsv(perform, paste0(output, model, "_performance.txt"))
