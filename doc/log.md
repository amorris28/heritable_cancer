# Log

## Workflow (updated 30.09.2024)

1. LOCAL: Run `PCA_sumstats.r` to download sumstats, preprocess them for cleansumstats, and create YAML files.
2. TSD: Copy `.bim` file from each TSD project to LOCAL.
    - p212: ssh into `p212-hpc-01` and find the .bim file at `/cluster/projects/p212/antenor_pipeline/antegenes_imputed_allchr-temporary.bim`
3. Copy the local repo including scripts, bim files, and sumstats files to p33.
3. p33: Run `run_cleansumstats_bulk.job` to harmonize all sumstats in `traitlist.txt`. Copy sumstats back to LOCAL.
4. p33: Run `run_cleansumstats_bim_bulk.job` to harmonize all bim files in `bimlist.txt`
5. Export these files to local and then import to their respective TSD project:
    `${BFILE_OUT}_tmp.bim,forceALT_betan_${BFILE_OUT}.txt`
6. p212: Run `rotate_alleles.batch` 

## To Do:

1. Check which RSIDs are missing from the bfiles (not just the bim file but the bed file too)
2. Use LDlink to find alternative SNPs
3. Write a script to substitute in these IDs

Scripts to update in local:

- `src/PCA_rotate_alleles.batch`
- `src/PCA_update_ids_antegenes.R`
- `src/PCA_identify_cases_and_controls.R`
- `src/PCA_sumstats.R`
- basically everything in `src`

Other files:
- Some lists in `data/interim`

## Log

### 04.10.2024

Today I got the script running for testing the PHS models. Next I should probably do a sync up between p212 and local just to backup and version control my changes. Then, I should modify the PHS_testing file to align with Alex's decisions about how to go about testing.

### 03.10.2024

I have rewritten the code to create the combined phenotype data frame. No idea where the original code is/was/went. Now, I have bfiles filtered by RSIDs and sample IDs, case/controls with age of diagnosis or age of last obs. I think now I can run the matrix multiplication script.

Okay, I did all of this but the sumstats are still in INTERIM. I feel like they shoul be in the root of PROCESSED. So, at the end of the cleansumstats pipeline I should move just the sumstat files to PROCESSED.

I updated the bfiles batch script to move those to processed, which feels janky but whatever. Now everything is in processed so I need to figure out how to run the PHS calculation script again.

Okay, I did that. Now I have genotypes, phenotypes, and PHS scores. This is everything to go into the PHS "package." Now what I need to do is actually write the scripts for testing these models in new cohorts.

It appears several of the R packages I used in this analysis previously are not on TSD (broom, survminer). I will either need to install them from source myself or rewrite the code to not rely on them. But that is a task for tomorrow.

### 02.10.2024

I actually don't think I can just plug in this piece of code to pull out the rsids and filter on that. There are actually multiple steps:

1. Check which RSIDs in the sumstats files are missing from the BED file and save these missing RSIDs somewhere
2. Manually (or maybe with a CLI?) using LDLink to find SNPs in LD with these SNPs
3. Manually edit the sumstats file to replace the RSIDs (and the alleles??) with the ones from the LD SNPs
4. Filter the bfiles on this new list

Alternatively, since that is a bit sketchy, I could just plug in the RSID code as a first pass...

I will do that. What I need to do next is fix my "find cases and controls" scripts to output the appropriate "antegenes_id_samples.txt" with 0s in the FID column and IDs in the IID column and stop outputting all of the intermediate files. For now, this is what I have done. I don't know why I was generating all of those intermeidate files, but I have removed those steps.

Next, I have to prepare the inputs to the script that calculate PHS scores. The last step before the PHS validation package. I do not know how I created the `combined_data_simplified.tsv` file which contains IDs, age, and status. Possibly manually from the XLSX files, but that is bad, I should rewrite it as a script.

### 01.10.2024

Keep track of changes to src files and data to track those changes in local. (Maybe only edit in local?) Added `rotate_alleles.batch` to `src` from the old pipeline. Also added a `logs` directory. Need to update batch script to use this. Makefile. I cannot remember why I made pgen files intermediate to the rotate alleles step in the 1KG-imputed data. I see now, it is because this dataset has a split chromosome so I had to make-pgen and sort-vars first.

I was able to re-implement the `rotate_alleles.batch` script. Next, I need to setup the script to filter for relevant samples and filter the bfiles by this sample ID list.

I have filtered the sample IDs.

### 30.09.2024

`cleansumstats_bim.sh` failed because I was missing a script, which I added, and now I am rerunning it. I also updated my local repo to reflect any changes in the scripts as they are on p33. Downloaded the resulting files.

Next step: Upload to p212 and run the code to rotate alleles in the B files

### 27.09.2024

After cleansumstats ran successfully, copied the `heritable_cancer` directory back to local and copied the `data/sumstats` directory into local. Also moved the updated `run_cleansumstats_bulk.sh` script into the local `src/cleansumstats`. Coped the familial cancer bim file from p212 to local and then to p33 plus the cleansumstats scripts. Ran cleansumstats on the bim file.

### 25.09.2024

Copied `data/interim` and `src/cleansumstats` to p33 and ran `run_cleansumstats_bulk.sh` to clean sumstats files. This script required some modifications to run correctly. Also, updated the preprocessing script, `PCA_sumstats.R` because the extra columns in the sumstats files was breaking cleansumstats.