#!/bin/bash
#export COMORMENT=/ess/p2521/cluster/users/espehage/github/comorment/containers
#export PLINK="singularity exec --bind /ess/p2521:/ess/p2521 --home $PWD $COMORMENT/containers/1.9.0/gwas.sif plink"
export PLINK="singularity exec --bind /ess/p212:/ess/p212 --home $PWD $COMORMENT/containers/singularity/gwas.sif plink"


#$PLINK \
#  --bfile /ess/p212/cluster/moba_qc_imputation/BA/release6/moba_imputed_all \
#  --score src/xx_merge_snps/phs601_heritable_cancer.csv 6 7 5 header sum \
#  --out src/xx_merge_snps/data_prs_heritable_cancer
#$PLINK \
#  --bfile /ess/p2521/data/durable/HuskMent/genotypes/imputed/v240619/husk_all.autosomes.v240619.plink1 \
#  --score src/xx_merge_snps/phs601_husk.csv 6 7 5 header sum \
#  --out src/xx_merge_snps/data_prs_husk
$PLINK \
  --bfile data/processed/chrALL.imputed.phs601 \
  --score data/processed/phs601_heritable_cancer.csv 6 7 5 header sum \
  --out data/processed/PHSes/data_prs_p212_chrALL
