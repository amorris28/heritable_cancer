# Log

## Workflow (updated 30.09.2024)

1. LOCAL: Run `PCA_sumstats.r` to download sumstats, preprocess them for cleansumstats, and create YAML files.
2. p####: Copy `.bim` file from each TSD project to LOCAL.
    - p212: ssh into `p212-hpc-01` and find the .bim file at `/cluster/projects/p212/antenor_pipeline/antegenes_imputed_allchr-temporary.bim`
3. Copy the local repo including scripts, bim files, and sumstats files to p33.
3. p33: Run `run_cleansumstats_bulk.job` to harmonize all sumstats in `traitlist.txt`. Copy sumstats back to LOCAL.
4. p33: Run `run_cleansumstats_bim_bulk.job` to harmonize all bim files in `bimlist.txt`
5. Export these files to local and then import to their respective TSD project
    `${BFILE_OUT}_tmp.bim,forceALT_betan_${BFILE_OUT}.txt`

## Log

### 30.09.2024

`cleansumstats_bim.sh` failed because I was missing a script, which I added, and now I am rerunning it. I also updated my local repo to reflect any changes in the scripts as they are on p33. Downloaded the resulting files.

Next step: Upload to p212 and run the code to rotate alleles in the B files

### 27.09.2024

After cleansumstats ran successfully, copied the `heritable_cancer` directory back to local and copied the `data/sumstats` directory into local. Also moved the updated `run_cleansumstats_bulk.sh` script into the local `src/cleansumstats`. Coped the familial cancer bim file from p212 to local and then to p33 plus the cleansumstats scripts. Ran cleansumstats on the bim file.

### 25.09.2024

Copied `data/interim` and `src/cleansumstats` to p33 and ran `run_cleansumstats_bulk.sh` to clean sumstats files. This script required some modifications to run correctly. Also, updated the preprocessing script, `PCA_sumstats.R` because the extra columns in the sumstats files was breaking cleansumstats.