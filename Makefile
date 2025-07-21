#################################################################################
# GLOBALS                                                                       #
#################################################################################

PROJECT_NAME = heritable_cancer

#################################################################################
# COMMANDS                                                                      #
#################################################################################

## Delete all external and interim data files
.PHONY: cleandata
cleandata:
	find data/interim/. -type f -name "*.yaml" -delete
	find data/interim/. -type f -name "*.txt.gz" -delete
#	find data/external/. -type f -name "*.txt.gz" -delete
# removed this because i don't want to accidentally delete datasets shared
# by colleagues

## Delete all slurm log files in the logs/ directory
.PHONY: cleanlogs
cleanlogs:
	find logs/. -type f -name "*.err" -delete
	find logs/. -type f -name "*.out" -delete

## Delete all output in reports/ directory
.PHONY: cleanreports
cleanreports:
	find reports/. -type f -name "*.tsv" -delete
	find reports/. -type f -name "*.png" -delete

## Bundle git repo
.PHONY: bundle
bundle:
	git bundle create main.bundle HEAD main

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
	-x $(PROJECT_NAME)/**/*.bim \
	-x $(PROJECT_NAME)/**/*.bed \
	-x $(PROJECT_NAME)/**/*.fam \
	-x $(PROJECT_NAME)/**/*.DS_Store
	
#################################################################################
# PROJECT RULES                                                                 #
#################################################################################

## Download and pre-process sumstats
.PHONY: sumstats
sumstats:
	Rscript src/01_preprocess-sumstats.R

## Preprocess metadata
.PHONY: metadata
metadata:
	sbatch src/01_p212_preprocess-metadata.batch

## Identify IDs in antegenes to fix bfiles later
.PHONY: antegenes-ids
antegenes-ids:
	Rscript 02_update-ids_antegenes.R

## Run cleansumstats pipeline
.PHONY: cleansumstats
cleansumstats:
	sbatch src/02_run-cleansumstats.batch; \
	bash src/02_record-rsids.sh

.PHONY: cleansumstats-bim
cleansumstats-bim:
	sbatch src/04_run-cleansumstats-bim.batch

## Calculate PHS values for each model using matrix multiplication
.PHONY: phs
phs:
	sbatch src/p212_06_test-phs.batch

