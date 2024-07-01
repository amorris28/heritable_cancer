# Project directory for the Prostate Cancer Polygenic/Multimodal Hazard Score manuscript

## Overview

Polygenic hazard score (PHS) model parameters downloaded from the [PGS Catalog](https://www.pgscatalog.org/) or from the supplementary data of individual publications. 

Uses the [cleansumstats wrapper](https://github.com/precimed/GWAS_SUMSTAT/tree/master/JF) from Julian (and dependent on the original [cleansumstats pipeline](https://github.com/BioPsyk/cleansumstats)) available using command `run_cleansumstats.sh` on TSD p33.

## Directories

`bash/` - contains `.sh` scripts and batch scripts
`R/` - contains `.R` scripts
`raw_data/` - Raw data in the form it was received. Do not modify.
`derived_data/` - Data that has been processed via scripts. Can be deleted/recreated.
`env/` - Python `conda` environment information
`renv/` - R `renv` environment information
`manuscript/`- Contains files needed to render the manuscript in Docx format.