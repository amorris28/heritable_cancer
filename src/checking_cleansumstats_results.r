
library(tidyverse)
library(readxl)

genotype_moba = read_tsv("data/processed/moba_imputed_all_cs.bim", col_names = c("CHR", "RSID", "DUMMY", "POS", "EFFECT", "OTHER"))
genotype_norgene = read_tsv("data/processed/norgene_merged_cs_filtered.bim", col_names = c("CHR", "RSID", "DUMMY", "POS", "EFFECT", "OTHER"))
pre_cs = read_excel("data/raw/01_geneticScore.finalModel.snplist.xlsx")
post_cs = read_tsv("data/processed/PHS601.GRCh37.gz")

pre_cs %>%
    rename(POS = `Position (GRCh37)`) %>%
    select(-`Position (GRCh38)`) ->
    pre_cs
pre_cs %>% 
    filter(SNP %in% post_cs$RSID) ->
    pre_cs_filtered
post_cs %>%
    filter(RSID %in% pre_cs_filtered$SNP) ->
    post_cs_filtered
# Check that pre-post cs matches in sumstats

pre_cs_filtered %>%
    select(snp = SNP, pre_b = Beta, pre_allele = Effect) ->
    pre

post_cs_filtered %>%
    select(snp = RSID, post_b = B, post_allele = EffectAllele) ->
    post
pre %>%
    full_join(post, by = join_by(snp)) ->
    pre_post

pre_post %>%
    filter(snp == "rs138213197")

pre_post %>%
    filter(pre_allele == post_allele) %>%
    mutate(pre_b = round(pre_b, 4)) %>%
    mutate(post_b = round(post_b, 4)) %>%
    filter(pre_b != post_b)
    # Ok, this shows that the matching alleles have the same effect direction

pre_post %>%
    filter(pre_allele != post_allele) %>%
    mutate(pre_b = round(pre_b, 6)) %>%
    mutate(post_b = round(post_b, 6)) %>%
    filter(pre_b == post_b)
    # Ok, this shows that the reversed alleles also have reversed beta

# Next I should check if effect allele matches in genotype and sumstats
genotype_moba %>%
    filter(RSID %in% post_cs_filtered$RSID) ->
    genotype_moba_filtered

genotype_moba_filtered %>%
    select(snp = RSID, genotype_allele = EFFECT) %>%
    full_join(pre_post, by = join_by(snp)) %>%
    filter(genotype_allele == pre_allele) %>%
    mutate(pre_b = round(pre_b, 6)) %>%
    mutate(post_b = round(post_b, 6)) %>%
    filter(pre_b != post_b)
    # This means that the alleles that were not flipped have the same effect direction before and after

genotype_moba_filtered %>%
    select(snp = RSID, genotype_allele = EFFECT) %>%
    full_join(pre_post, by = join_by(snp)) %>%
    filter(genotype_allele != pre_allele) %>%
    mutate(pre_b = round(pre_b, 6)) %>%
    mutate(post_b = round(post_b, 6)) %>%
    filter(pre_b == post_b)
    # This means the alleles that were switched also switched their effect direction. good

genotype_norgene %>%
    filter(RSID %in% post_cs_filtered$RSID) ->
    genotype_norgene_filtered


genotype_norgene_filtered %>%
    select(snp = RSID, genotype_allele = EFFECT) %>%
    full_join(pre_post, by = join_by(snp)) %>%
    filter(genotype_allele == pre_allele) %>%
    mutate(pre_b = round(pre_b, 6)) %>%
    mutate(post_b = round(post_b, 6)) %>%
    filter(pre_b != post_b)
    # This means that the alleles that were not flipped have the same effect direction before and after

genotype_norgene_filtered %>%
    select(snp = RSID, genotype_allele = EFFECT) %>%
    full_join(pre_post, by = join_by(snp)) %>%
    filter(genotype_allele != pre_allele) %>%
    mutate(pre_b = round(pre_b, 6)) %>%
    mutate(post_b = round(post_b, 6)) %>%
    filter(pre_b == post_b)
    # This means the alleles that were switched also switched their effect direction. good

genotype_norgene_filtered %>%
    select(snp = RSID, norgene_effect = EFFECT, norgene_other = OTHER) ->
    norgene
genotype_moba_filtered %>%
    select(snp = RSID, moba_effect = EFFECT, moba_other = OTHER) ->
    moba
norgene %>%
    full_join(moba, by = join_by(snp)) %>%
    filter(norgene_effect != moba_effect)


genotype_norgene %>%
    filter(RSID == "rs138213197")

pre_cs_filtered %>%
    arrange(desc(Beta)) ->
    pre_cs_filtered_sorted

post_cs_filtered %>%
    arrange(B) ->
    post_cs_filtered_sorted


post_cs %>%
    mutate(B_negative = B * -1) %>%
    arrange(desc(B_negative)) ->
    post_cs_sorted
pre_cs_filtered_sorted$Beta == post_cs_sorted$B_negative
