#!/bin/bash
export COMORMENT=/ess/p2521/cluster/users/espehage/github/comorment/containers
export PLINK="singularity exec --bind /ess/p2521:/ess/p2521 --home $PWD $COMORMENT/containers/1.9.0/gwas.sif plink"

$PLINK \
  --bfile /ess/p2521/data/durable/HuskMent/genotypes/imputed/v240619/husk_all.autosomes.v240619.plink1 \
  --score src/xx_merge_snps/phs601_husk.csv 6 7 5 header sum \
  --out src/xx_merge_snps/data_prs_husk

