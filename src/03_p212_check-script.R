library(tidyverse)
library(readxl)
library(genio)

source("src/00_config.R")

allData = read_excel("/ess/p212/home/p212-andrewmo/FW_Poly_tables/poly_final_220425_estimated_lastobs.xls", sheet = 1)
caseData = read_excel("/ess/p212/home/p212-andrewmo/FW_Poly_tables/poly_final_220425_estimated_lastobs.xls", sheet = 2)
rareVariants = read_excel("/ess/p212/home/p212-andrewmo/FW_Poly_tables/poly_final_220425_estimated_lastobs.xls", sheet = 3)
fam_file = read_table("data/raw/moba_imputed_all.fam", col_names = FALSE) %>% as_tibble()

phenotypes = read_tsv("data/processed/heritable_cancer_phenotypes.tsv")


# Check which subjects are MISSING genotype data
allData %>%
	pull("PID") ->
	all_subject_ids

phenotypes %>%
	mutate(ids_in_genotype = bfile_id %in% fam_file$X2) %>%
	filter(!ids_in_genotype) %>%
	pull(bfile_id) ->
	subjects_notin_genotype

# Check if WGS subjects match the WGS rare variant table

rareVariants %>%
	filter(SOURCE == "DECODE") %>%
	pull(PID) %>%
	unique() ->
	wgs_subjects_with_variants

allData %>%
	select(PID, WGS) %>%
	filter(WGS == 1) %>%
	pull(PID) %>%
	unique() ->
	wgs_subjects

# 247 of the WGS subjects are listed with a pathogenic variant in the rare variant table
wgs_subjects[wgs_subjects %in% wgs_subjects_with_variants]

# All of the rare variant subjects are marked as WGS in the phenotype data
wgs_subjects_with_variants %in% wgs_subjects

# There are 88 unique genes with variants identified in the WGS subjects
rareVariants %>%
	filter(SOURCE == "DECODE") %>%
	pull(GENE) %>%
	unique() ->
	unique_genes_in_wgs

# The table doesn't appear to list what the specific mutation is so I can't
# answer how many mutations or how many mutations per gene
rareVariants %>%
	filter(SOURCE == "DECODE") %>%
	print(n=1000)

# There are up to 16 mutations per gene
rareVariants %>%
	filter(SOURCE == "DECODE") %>%
	group_by(GENE) %>%
	summarize(n = n()) %>%
	arrange(desc(n)) ->
	gene_counts
gene_counts %>%
	print(n = 1000)

rareVariants %>%
	filter(SOURCE == "CGEN") %>%
	group_by(GENE) %>%
	summarize(n = n()) %>%
	arrange(desc(n)) ->
	gene_counts_cgen
gene_counts_cgen %>%
	print(n = 1000)
# One subject has 4 mutations listed, one has 3, and many have 2
rareVariants %>%
	filter(SOURCE == "DECODE") %>%
	group_by(PID) %>%
	summarize(n = n()) %>%
	arrange(desc(n)) ->
	subject_counts

# Three subjects have 2 mutations in the same gene
rareVariants %>%
	filter(SOURCE == "DECODE") %>%
	group_by(PID, GENE) %>%
	summarize(n = n()) %>%
	arrange(desc(n)) ->
	subject_gene_counts

# Some subjects are repeated multiple times within an individual source (decode, panel testing, or cgen)
rareVariants %>%
	select(PID, SOURCE) %>%
	group_by(PID, SOURCE) %>%
	summarize(n = n()) %>%
	arrange(desc(n))

# There were 2847 unique subjects with CGEN testing, possibly 68 with panel testing, and 247 with WGS from decode
rareVariants %>%
	select(PID, SOURCE) %>%
	unique() %>%
	group_by(SOURCE) %>%
	summarize(n = n()) ->
	n_subjects_per_source

rareVariants %>%
	filter(SOURCE == "CGEN") %>%
	select(GENE) %>%
	unique()