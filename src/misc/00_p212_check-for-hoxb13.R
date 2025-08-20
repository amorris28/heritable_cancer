library(tidyverse)
#library(genio)

#source("src/00_config.R")

# Input files
# bfiles
# phenotype data

###################################################
# commandArgs block
#
#
#genotype_files = paste0(PROCESSED_DATA, "moba_imputed_all_cs")
#genotype_files = "/ess/p212/cluster/moba_qc_imputation/BA/release6/moba_imputed_all"
#metadata_file = paste0(PROCESSED_DATA, "heritable_cancer_phenotypes.tsv")
#model_file = paste0(PROCESSED_DATA, "PHS601.GRCh37.gz")
###################################################
#phs_model = read_tsv(model_file)
hoxb13 = "rs138213197"

#phenotype = read_tsv(metadata_file) # Phenotype data

#genotype = genio::read_plink(genotype_files)$X
bimfile = read_tsv("/ess/p212/cluster/moba_qc_imputation/BA/release6/moba_imputed_all.bim", col_names = FALSE)
print(paste0("Does hoxb13 exist in the moba data?: ", hoxb13 %in% bimfile$X2))

bimfile = read_tsv("/ess/p212/cluster/antenor_pipeline/antegenes_imputed_allchr-temporary.bim", col_names = FALSE)
print(paste0("Does hoxb13 exist in the antenor data?: ", hoxb13 %in% bimfile$X2))


#genotype[rownames(genotype) == hoxb13, ]
