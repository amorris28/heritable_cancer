library(tidyverse)
library(readxl)
library(survival)

source("src/00_config.R")

allData = read_excel("/ess/p212/home/p212-andrewmo/FW_Poly_tables/poly_final_220425_estimated_lastobs.xls", sheet = 1)
caseData = read_excel("/ess/p212/home/p212-andrewmo/FW_Poly_tables/poly_final_220425_estimated_lastobs.xls", sheet = 2)
rareVariants = read_excel("/ess/p212/home/p212-andrewmo/FW_Poly_tables/poly_final_220425_estimated_lastobs.xls", sheet = 3)

phenotypes = read_tsv("data/processed/heritable_cancer_phenotypes.tsv")

phs = read_tsv(paste0(PROCESSED_DATA, "p212_PHS601_phs"), col_names = c("bfile_id", "phs"))

# Merge phenotype+phs data
phenotypes %>%
    full_join(phs, by = "bfile_id") ->
    full_data

# Remove the 129 people missing genotype data
full_data %>%
    filter(!is.na(phs)) ->
    full_data
# Filter individuals from rare variants in phenotype data
rareVariants %>%
    mutate(bfile_id = paste0("Eivind_", PID)) %>%
    filter(bfile_id %in% phenotypes$bfile_id) ->
    rare_variants_filtered

# It looks like the CGEN variants are classified as class 3 to 5
rare_variants_filtered %>%
    select(bfile_id, SOURCE:MANDY) %>%
    print(n = 1000)

# It appears the most frequent mutation among the PCa cases and controls
# are BRCA1 and 2 within the CGEN testing. 
# For DECODE, unfortunately, there are very few variants
rare_variants_filtered %>%
    select(bfile_id, SOURCE, GENE) %>%
    group_by(SOURCE, GENE) %>%
    summarize(n = n()) %>%
    arrange(desc(n)) %>%
    print(n = 100)

# Summarize rare variant data as 0/1 with any variant
rare_variants_filtered %>%
    select(bfile_id, SOURCE) %>%
    mutate(carrier = 1) ->
    any_variant

# It appears as though 6 subjects ahve 2 variants, all from "panel testing"
any_variant %>%
    group_by(bfile_id, SOURCE) %>%
    summarize(n = n()) %>%
    arrange(desc(n))
any_variant %>%
    group_by(bfile_id) %>%
    summarize(n = n()) %>%
    arrange(desc(n))

# remove source info
any_variant %>%
    select(bfile_id, carrier) %>%
    unique() ->
    any_variant

full_data %>%
    full_join(any_variant, by = "bfile_id") %>%
    mutate(carrier = if_else(is.na(carrier), 0, 1)) %>%
    filter(!is.na(phs)) -> #remove more missing genotype people
    combined_any_variant
sum(is.na(combined_any_variant$phs))
# Test if PHS predicts any variant
# n.s.
anyvariant_fit <- glm(carrier~phs, family="binomial", data=combined_any_variant)
summary(anyvariant_fit)
binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)
}
p = ggplot(combined_any_variant, aes(x = phs, y = carrier)) +
    geom_point() +
    binomial_smooth()
ggsave(paste0(FIGS, "binomial_rare_variant.png"), p)

# Test if PHS predicts any DECODE variant
rare_variants_filtered %>%
    select(bfile_id, SOURCE) %>%
    filter(SOURCE == "DECODE") %>%
    mutate(carrier = 1) %>%
    select(bfile_id, carrier) %>%
    unique() ->
    any_decode_variant
full_data %>%
    full_join(any_decode_variant, by = "bfile_id") %>%
    mutate(carrier = if_else(is.na(carrier), 0, 1)) %>%
    filter(!is.na(phs)) ->
    combined_any_decode_variant
any_decode_variant_fit <- glm(carrier~phs, family="binomial", data=combined_any_decode_variant)
summary(any_decode_variant_fit)

# How does coxph performance compare with PHS versus PHS and rare variants?
cph_fit_phs = coxph(Surv(age, status) ~ phs, combined_any_variant)
summary(cph_fit_phs)
cph_fit_phs_rv = coxph(Surv(age, status) ~ phs + carrier, combined_any_variant)
summary(cph_fit_phs_rv)

anova(cph_fit_phs, cph_fit_phs_rv, test = "LRT")
cph_fit_rv = coxph(Surv(age, status) ~ carrier, combined_any_variant)
summary(cph_fit_rv)
cph_fit_rv_phs = coxph(Surv(age, status) ~ carrier + phs, combined_any_variant)
summary(cph_fit_rv_phs)

anova(cph_fit_rv, cph_fit_rv_phs, test = "LRT")


# How does coxph performance compare with and without subjects with a rare variant?
combined_any_variant %>%
    filter(carrier == 0) ->
    no_rvs
cph_fit_phs_norvs = coxph(Surv(age, status) ~ phs, no_rvs)
summary(cph_fit_phs_norvs)
summary(cph_fit_phs)

# No difference in Concordance index if you remove those with a rare variant

# What about if we only use the decode data
cph_fit_phs = coxph(Surv(age, status) ~ phs, combined_any_decode_variant)
summary(cph_fit_phs)
cph_fit_phs_rv = coxph(Surv(age, status) ~ phs + carrier, combined_any_decode_variant)
summary(cph_fit_phs_rv)

anova(cph_fit_phs, cph_fit_phs_rv, test = "LRT")

cph_fit_rv = coxph(Surv(age, status) ~ carrier, combined_any_decode_variant)
summary(cph_fit_rv)
cph_fit_rv_phs = coxph(Surv(age, status) ~ carrier + phs, combined_any_decode_variant)
summary(cph_fit_rv_phs)

anova(cph_fit_rv, cph_fit_rv_phs, test = "LRT")

# How does coxph performance compare with and without subjects with a DECODE rare variant?
combined_any_decode_variant %>%
    filter(carrier == 0) ->
    no_rvs
cph_fit_phs_norvs = coxph(Surv(age, status) ~ phs, no_rvs)
summary(cph_fit_phs_norvs)
summary(cph_fit_phs)

# Again, no difference, pretty much exactly the same

# What about with only known PCa pathogenic variants
# just BRCA2
rare_variants_filtered %>%
    filter(GENE == "BRCA2") %>%
    mutate(carrier = 1) %>%
    select(bfile_id, carrier) ->
    just_BRCA2

full_data %>%
    full_join(just_BRCA2, by = "bfile_id") %>%
    mutate(carrier = if_else(is.na(carrier), 0, 1)) %>%
    filter(!is.na(phs)) ->
    combined_just_BRCA2

cph_fit_phs = coxph(Surv(age, status) ~ phs, combined_just_BRCA2)
summary(cph_fit_phs)
cph_fit_phs_rv = coxph(Surv(age, status) ~ phs + carrier, combined_just_BRCA2)
summary(cph_fit_phs_rv)

anova(cph_fit_phs, cph_fit_phs_rv, test = "LRT")

cph_fit_rv = coxph(Surv(age, status) ~ carrier, combined_just_BRCA2)
summary(cph_fit_rv)
cph_fit_rv_phs = coxph(Surv(age, status) ~ carrier + phs, combined_just_BRCA2)
summary(cph_fit_rv_phs)

anova(cph_fit_rv, cph_fit_rv_phs, test = "LRT")



# all pathogenic variants
pca_pathogenic_variants = c("BRCA2", "BRCA1", "HOXB13", "MLH1", "MSH2", "MSH6")
pca_pathogenic_variants = c("BRCA2", "BRCA1", "HOXB13", "MLH1", "MSH2", "MSH6", "SAMHD1", "ATM", "CHEK2", "FAM111A", "BIK")
rare_variants_filtered %>%
    filter(GENE %in% pca_pathogenic_variants) %>%
    mutate(carrier = 1) %>%
    select(bfile_id, carrier) ->
    pathogenic_variants

full_data %>%
    full_join(pathogenic_variants, by = "bfile_id") %>%
    mutate(carrier = if_else(is.na(carrier), 0, 1)) %>%
    filter(!is.na(phs)) ->
    combined_pathogenic_variants

cph_fit_phs = coxph(Surv(age, status) ~ phs, combined_pathogenic_variants)
summary(cph_fit_phs)
cph_fit_phs_rv = coxph(Surv(age, status) ~ phs + carrier, combined_pathogenic_variants)
summary(cph_fit_phs_rv)

anova(cph_fit_phs, cph_fit_phs_rv, test = "LRT")

cph_fit_rv = coxph(Surv(age, status) ~ carrier, combined_pathogenic_variants)
summary(cph_fit_rv)
cph_fit_rv_phs = coxph(Surv(age, status) ~ carrier + phs, combined_pathogenic_variants)
summary(cph_fit_rv_phs)

anova(cph_fit_rv, cph_fit_rv_phs, test = "LRT")