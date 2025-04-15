#library(hazrd)
library(survival)
# library(timeROC)
library(tidyverse)
# library(survivalROC)
source("src/00_config.R")

phenotypes = read_tsv("data/processed/heritable_cancer_phenotypes.tsv")
phs = read_tsv(paste0(PROCESSED_DATA, "p212_PHS601_phs"), col_names = c("bfile_id", "phs"))
phenotypes %>%
    full_join(phs, by = "bfile_id") %>%
    filter(!is.na(phs)) %>%
    filter(!is.na(age)) ->
    full_data

cox_model <- coxph(Surv(age, status) ~ phs, data = full_data)
risk_scores <- predict(cox_model)

# Specify the time point for evaluation
time_point <- 70

# Define a threshold (example: median risk score)
threshold <- median(risk_scores)

# Create predictions based on threshold
predicted_events <- ifelse(risk_scores > threshold, 1, 0)

# Create a confusion matrix at the specified time point
actual_events <- ifelse(full_data$status == 1 & full_data$age <= time_point, 1, 0)

# Calculate True Positives (TP), False Positives (FP), True Negatives (TN), False Negatives (FN)
TP <- sum(predicted_events == 1 & actual_events == 1)
TN <- sum(predicted_events == 0 & actual_events == 0)
FP <- sum(predicted_events == 1 & actual_events == 0)
FN <- sum(predicted_events == 0 & actual_events == 1)

# Calculate sensitivity and specificity
sensitivity <- TP / (TP + FN)
specificity <- TN / (TN + FP)

# Print results
cat("Sensitivity:", round(sensitivity, 3), "\n")
cat("Specificity:", round(specificity, 3), "\n")
