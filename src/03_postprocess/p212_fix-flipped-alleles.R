
library(tidyverse)
library(readxl)
hoxb13 = "rs138213197"
norgene = read_tsv("data/processed/norgene_merged_cs_filtered.bim", col_names = c("CHR", "RSID", "DUNNO", "POS", "EFFECT", "REF"))
moba = read_tsv("data/processed/moba_imputed_all_cs.bim", col_names = c("CHR", "RSID", "DUNNO", "POS", "EFFECT", "REF"))
original_model = read_excel("data/external/01_geneticScore.finalModel.snplist.xlsx")
genotype = genio::read_plink("data/processed/norgene_merged_cs_filtered")$X
sum(genotype[hoxb13, ])
sumstats = read_tsv("data/processed/sumstats_basenames", col_names = "model")

fix_genotype_ref = NULL
for (i in sumstats$model) {
    print(i)
    model = read_tsv(paste0("data/processed/", i, ".GRCh37.gz"))
    norgene %>% 
        inner_join(model) %>%
    mutate(new_effect = if_else(EFFECT == EffectAllele, EFFECT,
                            if_else(REF == EffectAllele, REF, NA))) %>%
    mutate(new_ref = if_else(REF == OtherAllele, REF,
                            if_else(EFFECT == OtherAllele, EFFECT, NA))) %>%
    filter(complete.cases(.)) %>%
    select(CHR, POS, RSID, new_effect, new_ref) %>%
    mutate(model = i) ->
    new_alleles
    fix_genotype_ref = rbind(fix_genotype_ref, new_alleles)
}

fix_genotype_ref %>%
    pivot_wider(names_from = model, values_from = new_effect) %>%
    unite("new_effect", PGS000067:PHS601, na.rm = T) %>%
    rowwise() %>%
    mutate(effect = unique(strsplit(new_effect, split = "_"))) %>%
    unnest(effect) %>%
    select(CHR, POS, RSID, new_ref, new_effect = effect) ->
    new_genotype_ref

norgene %>%
    full_join(new_genotype_ref) %>%
    mutate(EFFECT = new_effect) %>%
    mutate(REF = new_ref)

PHS601 = read_tsv("data/processed/PHS601.GRCh37.gz")

original_model %>%
    select(RSID = SNP, CHR = Chromosome, POS = `Position (GRCh37)`, orig_effect = Effect, orig_ref = Ref) ->
    original_model_clean

norgene %>%
    filter(RSID %in% PHS601$RSID) %>%
    select(CHR, POS, RSID, cs_norgene_effect = EFFECT, cs_norgene_ref = REF) ->
    norgene_filt

PHS601 %>%
    filter(RSID %in% norgene$RSID) %>%
    select(CHR, POS, RSID, cs_model_effect = EffectAllele, cs_model_other = OtherAllele, B) ->
    PHS601_filt

PHS601_filt %>%
    left_join(norgene_filt, by = join_by(CHR, POS, RSID)) %>%
    left_join(original_model_clean) ->
    combined
combined %>%
    select(orig_effect, cs_model_effect, cs_norgene_effect, RSID, B) %>%
    print(n = 20)

combined %>%
    mutate(new_effect = if_else(cs_norgene_effect == cs_model_effect, cs_norgene_effect,
                            if_else(cs_norgene_ref == cs_model_effect, cs_norgene_ref, NA))) %>%
    mutate(new_ref = if_else(cs_norgene_ref == cs_model_other, cs_norgene_ref,
                            if_else(cs_norgene_effect == cs_model_other, cs_norgene_effect, NA))) %>%
    filter(complete.cases(.)) %>%
    select(CHR, POS, RSID, new_effect, new_ref) ->
    new_alleles

norgene %>%
    left_join(new_alleles)
