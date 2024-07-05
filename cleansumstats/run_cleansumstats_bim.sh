#!/bin/bash
# args: $1: BFILE, $2: OUT_DIR, $3: WORK_DIR, $4: p33/p697

###### user input ######
export BFILE=${1}
export OUT_DIR=${2}
export WORK_DIR=${3}
export TSD_PROJ=${4}
export REF=GRCh37 # choose reference build (GRCh37 or GRCh38) according to the bim build

###### Initialisation ######
set -o errexit

# set vars dependendent on the TSD project
if [[ "${TSD_PROJ}" == "p697" ]]; then
    export COMORMENT=/ess/p697/data/durable/s3-api/github/comorment
    export CLEANSUMSTATS=/ess/p697/cluster/projects/cleansumstats_v1.7.0
    export PYTHON="singularity exec --bind /ess/p697:/ess/p697 --home $PWD:/home $COMORMENT/containers/singularity/python3.sif python"
    export PLINK2="singularity exec --bind /ess/p697:/ess/p697 --home $PWD:/home $COMORMENT/containers/singularity/gwas.sif plink2"
else 
    export COMORMENT=/ess/p33/cluster/github/comorment
    export CLEANSUMSTATS=/ess/p33/cluster/groups/biostat/cleansumstats_v1.7.0
    export PYTHON="singularity exec --bind /ess/p33:/ess/p33,/usr/bin:/usr/bin --home $PWD:/home $COMORMENT/containers/singularity/python3.sif python"
    export PLINK2="singularity exec --bind /ess/p33:/ess/p33 --home $PWD:/home $COMORMENT/containers/singularity/gwas.sif plink2"
fi

# Path and variable definition
export TRAIT=$(basename "${BFILE}") && echo Trait name: "${TRAIT}"
export WORK_DIR=${WORK_DIR}/${TRAIT} && mkdir -p ${WORK_DIR} # individual work directory for parallel execution on app nodes
export META_YAML="${WORK_DIR}"/"$(basename ${BFILE})".yaml # created in bim2sum.py
export BFILE_OUT="$(basename ${BFILE})"_cs


# Get script location
SOURCE=${BASH_SOURCE[0]}
while [ -L "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  SCRIPT_DIR=$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )
  SOURCE=$(readlink "$SOURCE")
  [[ $SOURCE != /* ]] && SOURCE=$SCRIPT_DIR/$SOURCE # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
SCRIPT_DIR=$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )

###### Pipeline ######
# Create meta data file and bim sumstat version
${PYTHON} ${SCRIPT_DIR}/bim2sum.py --bfile ${BFILE} --out ${WORK_DIR}

# Run the pipeline
cd ${CLEANSUMSTATS} && ./cleansumstats.sh -i ${META_YAML} -d ${CLEANSUMSTATS}/out_dbsnp -k ${CLEANSUMSTATS}/out_1kgp -o ${WORK_DIR} -w ${WORK_DIR}/work -b ${WORK_DIR}/tmp

# Concatenate cleaned_GRCh37.gz and cleaned_GRCh38.gz
cd ${WORK_DIR} && paste <(zcat cleaned_GRCh37.gz) <(cut -f 3,5- <(zcat cleaned_GRCh38.gz) ) | awk '$1 != "NA"' | gzip > ${TRAIT}.GRCh37.gz
cd ${WORK_DIR} && zcat cleaned_GRCh38.gz | awk '$1 != "NA"' | gzip > ${TRAIT}.GRCh38.gz

# Reshape output bim to original PLINK dimension
${PYTHON} ${SCRIPT_DIR}/sum2bim.py --bfile ${BFILE} --bim ${WORK_DIR}/${TRAIT}.${REF}.gz --out ${WORK_DIR}

# copy files from ${WORK_DIR} to ${OUT_DIR} for external PLINK2 call (optional)
mkdir -p ${OUT_DIR} && cp ${WORK_DIR}/{${BFILE_OUT}_tmp.bim,forceALT_betan_${BFILE_OUT}.txt} ${OUT_DIR}

# Rotate alleles in bed file which were flipped by cleansumstats (files needed for external PLINK call: ${BFILE_OUT}_tmp.bim and forceALT_betan_${BFILE_OUT}.txt)
#$PLINK2 --bed ${BFILE}.bed \
#    --bim ${WORK_DIR}/${BFILE_OUT}_tmp.bim \
#    --fam ${BFILE}.fam \
#    --alt1-allele force ${WORK_DIR}/forceALT_betan_${BFILE_OUT}.txt \
#    --make-bed \
#    --out ${OUT_DIR}/${BFILE_OUT}

# Remove line suffix from RSIDS for multi-allelic SNPs (optional)
#awk '{split($2,a,"_"); print $1 "\t" a[1] "\t" $3 "\t" $4 "\t" $5 "\t" $6}' ${OUT_DIR}/${BFILE_OUT}.bim > ${OUT_DIR}/${BFILE_OUT}.bim_tmp
#mv ${OUT_DIR}/${BFILE_OUT}.bim_tmp ${OUT_DIR}/${BFILE_OUT}.bim

# Test for bim_clsum==bim_plink
#echo -e '
#import pandas as pd
#bim_org=pd.read_csv("'${BFILE}.bim'", sep=None, header=None, names=["CHR","RSID","POSITION","BP","ALT","REF"], engine="python")
#bim_plink=pd.read_csv("'${OUT_DIR}/${BFILE_OUT}.bim'", sep=None, header=None, names=["CHR","RSID","POSITION","BP","ALT","REF"], engine="python")
#bim_clsum=pd.read_csv("'${WORK_DIR}/${TRAIT}.${REF}.gz'", sep=None, engine="python")
#bim_out=pd.read_csv("'${WORK_DIR}/${BFILE_OUT}'_tmp.bim", sep=None, names=["CHR","RSID","POSITION","BP","ALT","REF"], engine="python", index_col=0, nrows=0)
#bim_clsum=bim_clsum.rename(columns={"Direction": "POSITION", "POS": "BP", "EffectAllele": "ALT", "OtherAllele": "REF"})
#bim_clsum=bim_clsum.loc[bim_clsum["CHR"]==bim_clsum["CHR"].value_counts().index[0]] if bim_clsum["CHR"].value_counts().index.shape[0]>1 and bim_org["CHR"].value_counts().index.shape[0]<1 else bim_clsum
#suf=("_clsum", "_plink")
#ref=pd.merge(bim_clsum, bim_plink, on="RSID", suffixes=suf)
#check_clsum_plink=[ref[cols[0]].equals(ref[cols[1]]) for cols in zip(bim_out.columns.drop("RSID")+suf[0], bim_out.columns.drop("RSID")+suf[1])]
#assert ref.shape[0]==bim_clsum.shape[0] and ref.shape[0]==bim_plink.shape[0] and check_clsum_plink
#' | ${PYTHON} 
#if [ $? -eq 1 ]; then echo "Error: cleansumstats bim != plink bim." >&2;fi

echo Done!
