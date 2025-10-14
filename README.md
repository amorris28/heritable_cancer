# Prostate Cancer Polygenic Hazard Score in Norwegian Cohorts


## Data Availability

Polygenic hazard score (PHS) model parameters downloaded from the
[PGS Catalog](https://www.pgscatalog.org/) or through personal communication with the authors.

## Required software

Recreating the analysis requires a recent version of [Nextflow](https://github.com/nextflow-io/nextflow).
Analyses were developed using Nextflow [version 24.04.2](https://github.com/nextflow-io/nextflow/releases/tag/v24.04.2).

All R, python3, and plink calls were performed using
[the COMORMENT containers](https://github.com/comorment/containers) with either 
singularity or docker. Analyses were developed using the COMORMENT containers 
[version 1.9.0](https://github.com/comorment/containers/releases/tag/v1.9.0).

## Commands to reproduce analysis

To run the analysis, run the 
```
main
```




Scripts within `src/` are named with fields separated by underscores `_` and
words within a field separated by dashes `-`. The leading number is the TSD 
project number, i.e., the remote cluster where the data are saved and where
analyses were conducted:

`p212` - The heritable cancer cohort
`p2521` - The HUSK cohort

The following sequential numbers (e.g., `01`) is the
order scripts run. The third field is a slug to describe the analysis
step.

### Preprocessing steps specific to the PCA project

### 1. Download and preprocess sumstats

Update `traitlist.txt` with the URLs for your PHS models on PHS Catalog.

Run `make sumstats` to download the sumstats files from PGS, preprocess them
for input into cleansumstats, and automatically generate yaml files. 

### 2. Download the bim files manually from their respective TSD project

Download bim files from each TSD project and place in `data/raw`.

Update `bimlist.txt` with the absolute paths for your bim files (without file
extensions).

Run `make zip` and upload the zip file to p33.

### 3. Run the cleansumstats pipeline on p33

Unzip the zip file in your working directory.

Run `make cleansumstats`

### 4. Export files from p33 and reupload them to their respective TSD project

Export `${BFILE_OUT}_tmp.bim`, `forceALT_betan_${BFILE_OUT}.txt`, and
`sumstats.GRCh37.gz` from p33:

Run `make zip`. Download the zip file to your local PC. Then, upload the zip
file to each TSD project.

### 5. Prepare data before analysis on p212

We need to harmonize the old bed and fam files with the updated bim file. We
also need to identify subjects that are cases and controls, their censoring age
(diagnosis/last observed), and filtered the bfiles for only the subjects and
RSIDs relevant to us.

Do all of this by running `make metadata`

### 6. Calculate PHS scores

Run `make phs` to calculate PHS scores for each combination of model/bfile
using matrix multiplication.

This concludes the project-specific prep-processing. The remainded of the
analysis should be generalizable to any dataset.

## PHS Testing Pipeline

This batch script fits a K-M curve and calculate performance statistics.
Multiple such analysis for multiple PHS models can be run in parallel using the
`--array` argument to sbatch. The custome script `multijoin.sh` can merge the
resulting performance tables into one table.

```bash
sbatch bash/phs.batch
./bash/multijoin.sh derived_data/*performance.txt > performance.txt
```
### Test performance on p212

## Project Organization

```
├── Makefile           <- Makefile with convenience commands like `make 
│                         sumstats` or `make metadata`
├── README.md          <- The top-level README for developers using this
│                         project.
├── data
│   ├── external       <- Data from third party sources.
│   ├── interim        <- Intermediate data that has been transformed.
│   ├── processed      <- The final, canonical data sets for modeling.
│   └── raw            <- The original, immutable data dump.
│
├── notebooks          <- Quarto notebooks.
│
├── manuscript         <- Manuscript documents.
│
└── src                <- Source code for use in this project. scripts are
    │	                  numbered in the order they run.
    └── cleansumstats  <- Custom cleansumstats wrapper scripts
```

## Someday/Maybe

### Identify missing RSIDs and then fine alleles in LD with them

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

Use [LDproxy](https://ldlink.nih.gov/?tab=ldproxy) to find SNPs in LD with the
missing SNPs from the model. Some SNPs may be monoallelic in this population or
missing from the reference panel.


