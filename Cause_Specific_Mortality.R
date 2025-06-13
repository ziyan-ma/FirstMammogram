####--------------------------------------------------
#### Title: Cause-specific mortality
#### Author: Ziyan Ma
#### Last Updated: 2025-06-02
####--------------------------------------------------
# Note: "LopNr" is the unique identifier used consistently across datasets.

# Clear environment
rm(list = ls())

# Set working directory
setwd("")

# Load necessary libraries
library(tidyverse)
library(survival)

# Load data: screening records and cause of death information
load("mammography_screening_first_screening_defined.RData")
load("cause_death_since1989.RData")

# ------------------------------------------------------------------------------
# Merge data, define follow-up, outcomes, and recode variables
# ------------------------------------------------------------------------------

dat <- mammography_screening_first_screening_defined %>%
  # Merge with cause of death data
  left_join(cause_death_since1989) %>%
  group_by(LopNr) %>%
  # Define end of follow-up: death, emigration, or end of 2023
  mutate(end_of_followup = min(
    as.Date("2023-12-31", tryFormats = "%Y-%m-%d"), 
    date_death, EmigrationDate, 
    na.rm = TRUE)) %>%
  ungroup() %>%
  # Calculate follow-up time in years
  mutate(
    end_of_followup = as.Date(end_of_followup),
    followup_time = as.numeric(difftime(end_of_followup, RequestDate_1, units = "days")) / 365.25
  ) %>%
  # Define outcome based on cause of death
  mutate(
    outcome = ifelse(is.na(date_death) | date_death != end_of_followup, "censor", cause_death_group)
  ) %>%
  # Apply follow-up limit of 25 years
  filter(followup_time > 0) %>%
  mutate(
    outcome = ifelse(followup_time > 25, "censor", outcome),
    followup_time = pmin(followup_time, 25),
    # Binary outcome variables for Cox models
    outcome_all_cause_death = ifelse(outcome == "censor", 1, 2),
    outcome_nonBC = ifelse(outcome == "censor", 1, ifelse(outcome != "Breast_Cancer", 2, 1)),
    outcome_BC_death = ifelse(outcome == "Breast_Cancer", 2, 1),
    # Relevel participation factor (1 = Participant, 0 = Non-participant)
    participate_1 = factor(participate_1, levels = c(1, 0))
  )

# ------------------------------------------------------------------------------
# Handle missing values for common covariates
# ------------------------------------------------------------------------------

common_factors <- c(
  "Education_regroup_at_screening_1",
  "Income_individual_standardized_category_at_screening_1",
  "born_Sweden",
  "n_child_cate",
  "first_child_cate",
  "alcohol_related_diasease",
  "obesity_related_diasease",
  "BC_FH",
  "AllC_FH",
  "CCIw_Regroup",
  "Marital"
)

# Replace NAs with "missing"
dat[common_factors] <- lapply(dat[common_factors], function(x) ifelse(is.na(x), "missing", x))

# ------------------------------------------------------------------------------
# Fit Cox models for 3 outcomes: all-cause death, breast cancer death, non-BC death
# ------------------------------------------------------------------------------

Outcome <- c("outcome_all_cause_death", "outcome_BC_death", "outcome_nonBC")
output <- data.frame()

for (i in Outcome) {
  for (j in c(FALSE, TRUE)) {
    
    # Construct model formula
    if (j == FALSE) {
      formula <- paste0("Surv(followup_time, ", i, ") ~ participate_1 + age_at_invitation_1 + RequestYear_1")
    } else {
      formula <- paste0("Surv(followup_time, ", i, ") ~ participate_1 + age_at_invitation_1 + RequestYear_1 + ",
                        paste(common_factors, collapse = " + "))
    }
    
    # Fit Cox model
    cox <- coxph(as.formula(formula), data = dat)
    
    # Extract and annotate results
    temp <- data.frame(summary(cox)$coefficients)
    temp$outcome <- i
    temp$adjust_for_common_factor <- j
    temp$variables <- rownames(temp)
    
    # Combine results
    output <- rbind(output, temp)
  }
}

# ------------------------------------------------------------------------------
# Extract and format results
# ------------------------------------------------------------------------------

final_results <- output %>%
  filter(variables == "participate_10") %>%
  mutate(
    HR = exp(coef),
    HR_lb = exp(coef - 1.96 * se.coef.),
    HR_ub = exp(coef + 1.96 * se.coef.),
    
    # Create HR (95% CI) string
    HR_ub_lb = paste0(
      round(HR, 2), " (", round(HR_lb, 2), ", ", round(HR_ub, 2), ")"
    )
  ) %>%
  select(outcome, adjust_for_common_factor, HR_ub_lb)

# Print formatted results
print(final_results)