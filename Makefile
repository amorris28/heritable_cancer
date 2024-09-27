
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
