####--------------------------------------------------
#### Title: Define Study Population
#### Author: Ziyan Ma
#### Last Updated: 2025-06-04
####--------------------------------------------------
# Note: "LopNr" is the unique identifier used consistently across datasets.

# Clear workspace
rm(list = ls())  

# Load required libraries
library(tidyverse)
library(ggplot2)
library(survival)
library(survminer)

# Set working directory
setwd("")

# Source custom ggplot2 theme
source("ggplot_theme.R")

# Load data
load("mammography_screening_first_screening_defined.RData")  # Mammography screening register
load("cause_death_since1989.RData")                          # Cause of death register

# ------------------------------------------------------------------------------
# Prepare analysis dataset with follow-up time and breast cancer mortality outcome
# ------------------------------------------------------------------------------

dat <- mammography_screening_first_screening_defined %>%
  left_join(cause_death_since1989, by = "LopNr") %>% # Merge cause of death data
  group_by(LopNr) %>%
  mutate(
    # Define end of follow-up as the earliest of study end date, death, or emigration
    end_of_followup = min(
      as.Date("2023-12-31", tryFormats = "%Y-%m-%d"), 
      date_death, 
      EmigrationDate, 
      na.rm = TRUE
    )
  ) %>%
  ungroup() %>%
  mutate(
    end_of_followup = as.Date(end_of_followup, tryFormats = "%Y-%m-%d"),
    # Calculate follow-up time in years since first invitation
    followup_time = as.numeric(difftime(end_of_followup, RequestDate_1, units = "days")) / 365.25,
    # Define outcome
    outcome = case_when(
      is.na(date_death) | date_death != end_of_followup ~ "censor",
      date_death == end_of_followup & death_from_BC == TRUE ~ "death_BC",
      date_death == end_of_followup & death_from_BC == FALSE ~ "death_other",
      TRUE ~ NA
    ),
    # Censor follow-up times exceeding 25 years
    outcome = ifelse(followup_time > 25, "censor", outcome),
    followup_time = ifelse(followup_time > 25, 25, followup_time),
    # For survival analysis: 2 = event (breast cancer death), 1 = censored/other
    outcome_BC_death = ifelse(outcome == "death_BC", 2, 1),
    # Convert participation at first invited screening to factor
    participate_1 = factor(participate_1, levels = c(0, 1))
  )

# ------------------------------------------------------------------------------
# Model: Breast cancer death by first screening participation
# ------------------------------------------------------------------------------

fit <- survfit(Surv(followup_time, outcome_BC_death) ~ participate_1, data = dat)

# Plot
p <- ggsurvplot(
  fit,
  conf.int = TRUE,
  censor = FALSE,
  risk.table = TRUE,
  ylim = c(0, 0.0175),
  xlim = c(0, 25),
  break.time.by = 5,
  legend.labs = c("Non-participant", "Participant"),
  legend.title = "First screening participation",
  xlab = "Time from first screening invitation (years)",
  ylab = "Breast cancer deaths per 10,000 women",
  fun = "event",  # Plot cumulative incidence
  risk.table.col = "strata"
)

# Apply themes
p$table <- p$table + theme_cleantable()
p$plot <- p$plot + theme_mine +
  scale_x_continuous(limits = c(0, 25), breaks = seq(0, 25, 5)) +
  scale_y_continuous(limits = c(0, 0.0120), breaks = seq(0, 0.0120, 0.0020))

