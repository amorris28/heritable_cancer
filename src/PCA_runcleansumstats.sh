#!/bin/bash

test $SCRATCH && module load singularity/3.7.3 

# choose location (both p33 and p697 project have the pipeline deployed)
export CLEANSUMSTATS=/ess/p697/cluster/projects/cleansumstats_v1.7.0
export CLEANSUMSTATS=/ess/p33/cluster/groups/biostat/cleansumstats_v1.7.0

export METADATA_YAML=<your .yaml file>
export OUT_FOLDER=<your output folder>

#no need to define for SLURM jobs - i.e. this can be job's $SCRATCH
export SCRATCH=<your scratch area>

cd ${CLEANSUMSTATS} && ./cleansumstats.sh \
-i ${METADATA_YAML} \
-d ${CLEANSUMSTATS}/out_dbsnp \
-k ${CLEANSUMSTATS}/out_1kgp \
-o ${OUT_FOLDER} \
-w $SCRATCH/work