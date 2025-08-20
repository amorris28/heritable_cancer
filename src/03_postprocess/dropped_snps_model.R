library(tidyverse)
library(readxl)
library(genio)

pre = read_excel("data/external/01_geneticScore.finalModel.snplist.xlsx")
post = read_tsv("data/processed/PHS601.GRCh37.gz")
genotype = genio::read_plink("data/processed/moba_imputed_all_cs")$X
#########################################################################
# The model SNPs went through cleansumstats first to standardize the RSIDs
# and remove SNPs. The genotype data went through cleansumstats to do the same
# then the processed genotype data was filtered for SNPs in the models. So,
# a few numbers that I need to report: 1. # SNPs in the original genotype data
# 2. # SNPs in genotype after cleansumstats
# 3. # of SNPs in the original model PHS601 (601)
# 4. # of SNPs in the model after css
# 5. # of SNPs from model that are missing in genotype data


geno_post = rownames(genotype)

geno_post[geno_post %in% pre$SNP]

pre_snps = nrow(pre)

post_snps = nrow(post)

removed_css = pre_snps - post_snps

missing_geno = nrow(post) - length(geno_post[geno_post %in% post$RSID])

end_snps = nrow(post) - missing_geno

css_summary = tibble(pre_snps, post_snps, removed_css, missing_geno, end_snps)

write_tsv(css_summary, "data/results/p212_PHS601_css_summary.tsv")
