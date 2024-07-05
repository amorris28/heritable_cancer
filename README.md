# Project directory for the Prostate Cancer Polygenic/Multimodal Hazard Score manuscript

## To-do

- Currently, all scripts require modification based on specific filenames to run. Make this a command line argument instead so they are generalizable.
- Add output filename option as well

## Overview

Polygenic hazard score (PHS) model parameters downloaded from the [PGS Catalog](https://www.pgscatalog.org/) or from the supplementary data of individual publications. 

Uses the [cleansumstats wrapper](https://github.com/precimed/GWAS_SUMSTAT/tree/master/JF) from Julian (and dependent on the original [cleansumstats pipeline](https://github.com/BioPsyk/cleansumstats)) available using command `run_cleansumstats.sh` on TSD p33.

## Steps to reproduce analysis

### 1. Remove leading comment lines

The cleansumstats pipeline may not handle leading comment lines beginning with "#"
so run this command on each sumstat file to remove them.

```bash
gzcat raw_data/sumstats/PGS000067_hmPOS_GRCh37.txt.gz \
  | sed '/^#/d' \
  | gzip \
  > derived_data/sumstats/PRACTICAL_PCa_Seibert2018.txt.gz
```

### 2. Run the cleansumstats pipeline on p33

You can use the cleansumstats pipeline driectly, but I have modified this
wrapper script to facilitate running the pipeline on TSD p33. Place all scripts from `cleansumstats/` in a directory on p33 and modify the `.job` scripts to run in slurm.

```bash
sbatch run_cleansumstats_bulk.job
sbatch run_cleansumstats_bim_bulk.job
```

Export `${BFILE_OUT}_tmp.bim`, `forceALT_betan_${BFILE_OUT}.txt`, and `sumstats.GRCh37.gz` from p33.
Import to the local p project for your B files (e.g., p212).

### Prepare data before analysis on p212

Identify sample IDs of prostate cancer cases and cancer-free controls.

```bash
module load R/4.2.1-foss-2022a
Rscript identify_cases_controls.R
````

For the antegenes/1KG imputation Bfiles: Update the Sample IDs since the fid and iid's have been created incorrectly.

```bash
Rscript update_ids_antegenes.R
````

Identify RSIDs of SNPs in models. sub `gzcat` if on macos.

```bash
SUMSTATS=derived_data/sumstats/
zcat ${SUMSTATS}PRACTICAL_PCa_Seibert2018.GRCh37.gz | awk 'NR > 1 { print $3 }' > derived_data/rsids_seibert
zcat ${SUMSTATS}PRACTICAL_PCa_Huynh-Le2021.GRCh37.gz | awk 'NR > 1 { print $3 }' > derived_data/rsids_huynle2021
zcat ${SUMSTATS}PRACTICAL_PCa_Huynh-Le2022.GRCh37.gz | awk 'NR > 1 { print $3 }' > derived_data/rsids_huynle2022
zcat ${SUMSTATS}PRACTICAL_PCa_Karunamuni2021.GRCh37.gz | awk 'NR > 1 { print $3 }' > derived_data/rsids_karunamuni

awk '{ print $2 }' out/moba_imputed_all_cs.bim > derived_data/RSIDS_MOBA
awk '{ print $2 }' out/antegenes_imputed_allchr-temporary_cs.bim > derived_data/RSIDS_ANTE
```

Compare RSIDs in sumstats to Bfiles and output missing SNPs.

```bash
Rscript R/compare_rsids.R
```

Use [LDproxy](https://ldlink.nih.gov/?tab=ldproxy) to find SNPs in LD with the missing SNPs from the model. Some SNPs may be monoallelic in this population or missing from the reference panel.



Use the output from `run_cleansumstats_bim.sh` to rotate alleles in BED and FAM files. Filter BIM file for cases/controls and SNPs in sumstats. 

```bash
sbatch rotate_alleles.batch
```

### Calculate PHS scores, fit model, and test performance on p212

This batch script calls two R scripts, one to calculate PHS scores for one set of sumstats + genotype and one to fit a K-M curve and calculate performance statistics. Multiple such analysis for multiple PHS models can be run in parallel using the `--array` argument to sbatch. The custome script `multijoin.sh` can merge the resulting performance tables into one table.

```bash
sbatch bash/phs.batch
./bash/multijoin.sh derived_data/*performance.txt > performance.txt
```

## Directories

`bash/` - contains `.sh` scripts and batch scripts
`R/` - contains `.R` scripts
`raw_data/` - Raw data in the form it was received. Do not modify.
`derived_data/` - Data that has been processed via scripts. Can be deleted/recreated.
`env/` - Python `conda` environment information
`renv/` - R `renv` environment information
`manuscript/`- Contains files needed to render the manuscript in Docx format.