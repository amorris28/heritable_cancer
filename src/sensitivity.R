library(survival)
library(tidyverse)
source("src/00_config.R")

phenotypes_p212 = read_tsv("data/processed/heritable_cancer_phenotypes.tsv")
phs_p212 = read_tsv(paste0(PROCESSED_DATA, "p212_PHS601_phs"), col_names = c("bfile_id", "phs"))
phenotypes_p2521 = read_tsv(paste0(PROCESSED_DATA, "husk_phenotypes.tsv"))
phs_p2521 = read_tsv(paste0(PROCESSED_DATA, "p2521_PHS601_phs"), col_names = c("subjectid", "phs"))

phenotypes_p212 %>%
    full_join(phs_p212, by = "bfile_id") %>%
    filter(!is.na(phs)) %>%
    filter(!is.na(age)) ->
    full_data

phenotypes_p2521 %>%
    full_join(phs_p2521, by = "subjectid") %>%
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

# Calculate True Positives (TP), False Positives (FP), 
# True Negatives (TN), False Negatives (FN)
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

# Calculate True Positives (TP), False Positives (FP), 
# True Negatives (TN), False Negatives (FN)
TP <- sum(predicted_events == 1 & actual_events == 1)
TN <- sum(predicted_events == 0 & actual_events == 0)
FP <- sum(predicted_events == 1 & actual_events == 0)
FN <- sum(predicted_events == 0 & actual_events == 1)

# Calculate Positive Predictive Value (PPV) and Negative Predictive Value (NPV)
PPV <- TP / (TP + FP)  # Proportion of positive predictions that are true positives
NPV <- TN / (TN + FN)  # Proportion of negative predictions that are true negatives

# Print results
cat("PPV:", round(PPV, 3), "\n")
cat("NPV:", round(NPV, 3), "\n")

