# Log

## Workflow (updated 27.09.2024)

1. LOCAL: Run `PCA_sumstats.r` to download sumstats, preprocess them for cleansumstats, and create YAML files. Copy the directory to p33
2. p33: Run `run_cleansumstats_bulk.sh` to harmonize all sumstats in `traitlist.txt`

## Log


### 27.09.2024

After cleansumstats ran successfully, copied the `heritable_cancer` directory back to local and copied the `data/sumstats` directory into local. Also moved the updated `run_cleansumstats_bulk.sh` script into the local `src/cleansumstats`.

### 25.09.2024

Copied `data/interim` and `src/cleansumstats` to p33 and ran `run_cleansumstats_bulk.sh` to clean sumstats files. This script required some modifications to run correctly. Also, updated the preprocessing script, `PCA_sumstats.R` because the extra columns in the sumstats files was breaking cleansumstats.