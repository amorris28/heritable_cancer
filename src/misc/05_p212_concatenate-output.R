library(tidyverse)
source("src/00_config.R")

#########################################
# First argument: list of model/sumstat basenames

args = commandArgs(trailingOnly=TRUE)

model_list = args[1] # sumstats_basenames

report = tibble()

for (model in model_list$X1) {
	print(paste0("Print report", report))
	performance = readr::read_tsv(paste0("reports/", model, "_performance.tsv"))
	print(paste0("Print ", model, " performance:", performance))
	diagnostics = readr::read_tsv(paste0("reports/", model, "_diagnostics.tsv"))
	print(paste0("Print ", model, " diagnostics:", diagnostics))
	nmissing = nrow(diagnostics)
	print(paste0("Print ", model, " nmissing:", nmissing))
	performance = cbind(performance, nmissing)
	print(performance)
	report = rbind(report, performance)
}

write_tsv(report, "reports/report.tsv")
