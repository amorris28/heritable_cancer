#!/bin/bash
export INTERIM_DIR=${PWD}/data/interim

# Pull out model SNP RSID
while read VARIABLE;
do
  echo $VARIABLE
  zcat data/interim/sumstats/${VARIABLE}/${VARIABLE}.GRCh37.gz | awk '{ print $3 }' | tail -n +2 >> data/processed/phs_rsids_GRCh37
  cp data/interim/sumstats/${VARIABLE}/${VARIABLE}.GRCh37.gz data/processed/${VARIABLE}.GRCh37.gz
done < data/processed/sumstats_basenames

# Pull out model SNP RSID
while read VARIABLE;
do
  echo $VARIABLE
  zcat data/interim/sumstats/${VARIABLE}/${VARIABLE}.GRCh38.gz | awk '{ print $3 }' | tail -n +2 >> data/processed/phs_rsids_GRCh38
  cp data/interim/sumstats/${VARIABLE}/${VARIABLE}.GRCh38.gz data/processed/${VARIABLE}.GRCh38.gz
done < data/processed/sumstats_basenames
