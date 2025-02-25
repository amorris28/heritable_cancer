#!/bin/bash
export INTERIM_DIR=${PWD}/data/interim

# Pull out model SNP RSIDs
while read VARIABLE;
do
  echo $VARIABLE
  zcat data/interim/sumstats/${VARIABLE}/${VARIABLE}.GRCh37.gz | awk '{ print $3 }' | tail -n +2 >> data/processed/phs_rsids
  cp data/interim/sumstats/${VARIABLE}/${VARIABLE}.GRCh37.gz data/processed/${VARIABLE}.txt.gz
done < data/processed/sumstats_basenames
