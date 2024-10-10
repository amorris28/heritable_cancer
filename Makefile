#################################################################################
# GLOBALS                                                                       #
#################################################################################

PROJECT_NAME = heritable_cancer

#################################################################################
# COMMANDS                                                                      #
#################################################################################

## Delete all external and interim data files
.PHONY: clean_data
clean_data:
	find data/external/. -type f -name "*.txt.gz" -delete
	find data/interim/. -type f -name "*.yaml" -delete
	find data/interim/. -type f -name "*.txt.gz" -delete

## Delete all slurm log files in the logs/ directory
.PHONY: clean_logs
clean_logs:
	find logs/. -type f -name "*.err" -delete
	find logs/. -type f -name "*.out" -delete

## Zip contents of directory for upload/download to/from TSD
.PHONY: zip
zip:
	cd ..; \
	zip -r $(PROJECT_NAME).zip \
	$(PROJECT_NAME)/data \
	$(PROJECT_NAME)/src \
	$(PROJECT_NAME)/Makefile \
	$(PROJECT_NAME)/README.md \
	$(PROJECT_NAME)/.gitignore \
	$(PROJECT_NAME)/.git \
	$(PROJECT_NAME)/renv \
	$(PROJECT_NAME)/renv.lock \
	-x $(PROJECT_NAME)/**/*.bed \
	-x $(PROJECT_NAME)/**/*.fam \
	-x $(PROJECT_NAME)/**/*.DS_Store
	
#################################################################################
# PROJECT RULES                                                                 #
#################################################################################

## Download and pre-process sumstats
.PHONY: sumstats
sumstats:
	Rscript src/01_pca_preprocess-sumstats.R

## Run cleansumstats pipeline
.PHONY: cleansumstats
cleansumstats:
	sbatch src/cleansumstats/run_cleansumstats_bulk.job
	sbatch src/cleansumstats/run_cleansumstats_bim_bulk.job

## Identify cases/controls and harmonize in bfiles
.PHONY: metadata
metadata:
	sbatch 03_pca_preprocess-metadata.batch

## Calculate PHS values for each model using matrix multiplication
.PHONY: phs
phs:
	sbatch 06_pca_calculate-phs.batch

