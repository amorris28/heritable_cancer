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

#################################################################################
# PROJECT RULES                                                                 #
#################################################################################

## Download and pre-process sumstats
.PHONY: sumstats
sumstats:
	Rscript src/PCA_sumstats.R

## Run cleansumstats pipeline
.PHONY: cleansumstats
cleansumstats:
	sbatch src/cleansumstats/run_cleansumstats_bulk.job
	sbatch src/cleansumstats/run_cleansumstats_bim_bulk.job

.PHONY: zip
zip:
	cd ..; \
	zip -r $(PROJECT_NAME).zip \
	$(PROJECT_NAME)/data \
	$(PROJECT_NAME)/src \
	$(PROJECT_NAME)/Makefile \
	$(PROJECT_NAME)/.gitignore \
	$(PROJECT_NAME)/.git \
	-x $(PROJECT_NAME)/**/*.bed \
	-x $(PROJECT_NAME)/**/*.fam \