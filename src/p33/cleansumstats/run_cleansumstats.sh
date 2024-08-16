#!/bin/bash
# args: $1: meta, $2: OUT_DIR, $3: WORK_DIR, $4: p33/p697, $5: SUMSTAT_DIR

###### user input ######
export META_YAML=${1}
export OUT_DIR=${2}
export WORK_DIR=${3}
export TSD_PROJ=${4}
export SUMSTAT_DIR=${5} # optional

###### Initialisation ######
set -o errexit

# set vars dependendent on the TSD project
if [[ "${TSD_PROJ}" == "p697" ]]; then
    export COMORMENT=/ess/p697/data/durable/s3-api/github/comorment
    export CLEANSUMSTATS=/ess/p697/cluster/projects/cleansumstats_v1.7.0
    export PYTHON="singularity exec --bind /ess/p697:/ess/p697 --home $PWD:/home $COMORMENT/containers/singularity/python3.sif python"
    export RSCRIPT="singularity exec --bind /ess/p697:/ess/p697 --home $PWD:/home $COMORMENT/containers/singularity/r.sif Rscript"
else 
    export COMORMENT=/ess/p33/cluster/github/comorment
    export CLEANSUMSTATS=/ess/p33/cluster/groups/biostat/cleansumstats_v1.7.0
    export PYTHON="singularity exec --bind /ess/p33:/ess/p33 --home $PWD:/home $COMORMENT/containers/singularity/python3.sif python"
    export RSCRIPT="singularity exec --bind /ess/p33:/ess/p33 --home $PWD:/home $COMORMENT/containers/singularity/r.sif Rscript"
fi
export REF_FILE="$COMORMENT/matlabruntime/reference/pleiofdr/9545380_ref/9545380.ref"

# Path and variable definition
export TRAIT=$(basename "${META_YAML%.*}") && echo Trait name: "${TRAIT}"
export WORK_DIR=${WORK_DIR}/${TRAIT} && mkdir -p ${WORK_DIR} # individual work directory for parallel execution on app nodes
export META_YAML_TMP=${WORK_DIR}/${TRAIT}.yaml

# get script location
SOURCE=${BASH_SOURCE[0]}
while [ -L "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  SCRIPT_DIR=$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )
  SOURCE=$(readlink "$SOURCE")
  [[ $SOURCE != /* ]] && SOURCE=$SCRIPT_DIR/$SOURCE # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
SCRIPT_DIR=$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )

###### Pre-processing ######
if [[ -z "${SUMSTAT_DIR}" ]]; then
   $PYTHON ${SCRIPT_DIR}/prep_sum.py --meta ${META_YAML} --out ${WORK_DIR} 
   else
   $PYTHON ${SCRIPT_DIR}/prep_sum.py --meta ${META_YAML} --out ${WORK_DIR}  --sumstat ${SUMSTAT_DIR}
fi

###### Run cleansumstats ######
cd ${CLEANSUMSTATS} && ./cleansumstats.sh -i ${META_YAML_TMP} -d ${CLEANSUMSTATS}/out_dbsnp -k ${CLEANSUMSTATS}/out_1kgp -o ${WORK_DIR} -w ${WORK_DIR}/work -b ${WORK_DIR}/tmp

# Concatenate cleaned_GRCh37.gz and cleaned_GRCh38.gz
cd ${WORK_DIR} && paste <(zcat cleaned_GRCh37.gz) <(cut -f 5- <(zcat cleaned_GRCh38.gz) ) | awk '$1 != "NA"' | gzip >  cleaned_GRCh37.sumstats.gz

###### Post-processing ######
mkdir -p "${OUT_DIR}"/"${TRAIT}" && echo Writing files to output directory: "${OUT_DIR}"/"${TRAIT}" # create subfolder in output folder using yaml prefix

# the output file from cleansumstats
FILE=cleaned_GRCh37.sumstats.gz

# copy cleaned_GRCh37.sumstats.gz to out directory
cp ${FILE} ${OUT_DIR}/${TRAIT}/${TRAIT}.GRCh37.gz

# copy cleaned_GRCh38.gz to out directory
cp cleaned_GRCh38.gz ${OUT_DIR}/${TRAIT}/${TRAIT}.GRCh38.gz

# no MHC
$PYTHON ${SCRIPT_DIR}/remove_region.py --input ${WORK_DIR}/$FILE --output ${OUT_DIR}/${TRAIT}/${TRAIT}.GRCh37.nomhc.gz # --chr 6 --start 26000000 --end 34000000

# mat file export
#$PYTHON ${SCRIPT_DIR}/sumstats2mat.py --sumstats ${WORK_DIR}/$FILE --ref ${REF_FILE} --out ${OUT_DIR}/${TRAIT}/${TRAIT}.ref9545380.mat --force    

# QQ and Manhatten plot
#$RSCRIPT ${SCRIPT_DIR}/qq_manhattan.R ${WORK_DIR}/$FILE ${OUT_DIR}/${TRAIT}/${TRAIT}

# copy details of cleansumstats pipeline to the output folder
mkdir ${OUT_DIR}/${TRAIT}/log_cleansumstats && echo Writing cleansumstats log files to ${OUT_DIR}/${TRAIT}/log_cleansumstats
zcat cleaned_GRCh38.gz | cut -f 3 | tail -n+2 | gzip > ${OUT_DIR}/${TRAIT}/log_cleansumstats/removed_lines_map_to_raw.gz
cp details/{stepwise_overview.txt,selected_source_stats.txt,removed_lines_per_type_table.txt,removed_lines.gz} ${OUT_DIR}/${TRAIT}/log_cleansumstats

echo Done!
