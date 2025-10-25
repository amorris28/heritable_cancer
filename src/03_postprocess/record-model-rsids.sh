#!/bin/bash
export INTERIM_DIR=${PWD}/data/interim

## Pull out model SNP RSID
#while read VARIABLE;
#do
#  echo $VARIABLE
#  zcat ${INTERIM_DIR}/sumstats/${VARIABLE}/${VARIABLE}.GRCh37.gz | awk '{ print $3 }' | tail -n +2 >> ${INTERIM_DIR}/phs_rsids
#  cp ${INTERIM_DIR}/sumstats/${VARIABLE}/${VARIABLE}.GRCh37.gz data/processed/${VARIABLE}.GRCh37.gz
#done < data/processed/sumstats_basenames
# Pull out model SNP RSID
while read VARIABLE;
do
  echo $VARIABLE
  gzcat ${INTERIM_DIR}/${VARIABLE}.GRCh37.gz | awk '{ print $3 }' | tail -n +2 >> ${INTERIM_DIR}/phs_rsids
#  cp ${INTERIM_DIR}/${VARIABLE}.GRCh37.gz data/processed/${VARIABLE}.GRCh37.gz
done < data/processed/sumstats_basenames
