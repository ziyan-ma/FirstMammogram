####--------------------------------------------------
#### Title: Generate Figure 3
#### Author: Ziyan Ma
#### Last Updated: 2025-06-07
####--------------------------------------------------
# Note: "LopNr" is the unique identifier used consistently across datasets.

# Clear workspace
rm(list = ls())

# Load required libraries
library(tidyverse)
library(ggplot2)
library(ggsurvfit)
library(survival)
library(survminer)
library(ggsci)

# Set working directory
setwd("")

# Load custom ggplot2 theme
source("ggplot_theme.R")

# Load mammography screening data
load("mammography_screening_first_screening_defined.RData")

# ------------------------------------------------------------------------------
# Define end of follow-up and classify breast cancer detection outcomes
# ------------------------------------------------------------------------------

dat <- mammography_screening_first_screening_defined %>%
  group_by(LopNr) %>%
  mutate(
    # Define end of follow-up: earliest of cancer diagnosis, last request + 2 years, death, emigration, or Dec 31, 2020
    end_of_followup = min(
      DiagnoseDate_swecan_1st_cancer_ever, 
      last_RequestDate_plus_2y, 
      deaddate, 
      EmigrationDate, 
      as.Date("2020-12-31", tryFormats = "%Y-%m-%d"), 
      na.rm = TRUE
    )
  ) %>%
  ungroup() %>%
  mutate(
    end_of_followup = as.Date(end_of_followup, tryFormats = "%Y-%m-%d"),
    # Calculate follow-up time in years
    followup_time = as.numeric(difftime(end_of_followup, RequestDate_1, units = "days")) / 365.25
  ) %>%
  filter(followup_time > 0) %>%
  mutate(
    # Define competing risks outcome variable based on diagnosis timing and mode of detection
    outcome = case_when(
      is.na(diagnosis_date) ~ "censor",
      diagnosis_date == end_of_followup & mode_of_detection_clinical == TRUE ~ "clinical_cancer",
      diagnosis_date == end_of_followup & mode_of_detection == "SDC" ~ "SDC",
      diagnosis_date == end_of_followup & mode_of_detection == "IC" ~ "IC",
      TRUE ~ "censor"
    )
  )

# Set outcome as factor for competing risks
dat$outcome <- factor(dat$outcome, levels = c("censor", "clinical_cancer", "SDC", "IC"))

# ------------------------------------------------------------------------------
# Fit a competing risks survival model
# ------------------------------------------------------------------------------

fitCR <- survfit(Surv(followup_time, outcome) ~ participate_1, data = dat, conf.int = FALSE)

# Extract summary of cumulative hazards
df <- summary(fitCR)

# Select relevant columns from summary
cols <- lapply(c(2:6, 8, 12), function(x) df[x])
tbl <- do.call(data.frame, cols) %>%
  # Calculate cumulative incidence per 100 women for each event type
  mutate(
    cumhaz_total = (cumhaz.1.2 + cumhaz.1.3 + cumhaz.1.4) * 100,
    cumhaz_clinical_cancer = cumhaz.1.2 * 100,
    cumhaz_other_SDC = cumhaz.1.3 * 100,
    cumhaz_other_IC = cumhaz.1.4 * 100,
    # Combine clinical and interval cancers for grouped area plotting
    cumhaz_clinical_IC = cumhaz_clinical_cancer + cumhaz_other_IC,
    # Rename strata for clarity in plot
    strata = ifelse(strata == "participate_1=0", "Non-participant", "Participant")
  )

# ------------------------------------------------------------------------------
# Plot: Cumulative incidence of breast cancer by detection mode and participation
# ------------------------------------------------------------------------------

p <- ggplot(tbl) +
  xlim(0, 25) +
  scale_y_continuous(limits = c(0, 9), breaks = seq(0, 9, 1)) +
  
  # Layered area plot (stacking from bottom to top)
  geom_area(aes(x = time, y = cumhaz_total, fill = "SDC")) +
  geom_area(aes(x = time, y = cumhaz_clinical_IC, fill = "IC")) +
  geom_area(aes(x = time, y = cumhaz_clinical_cancer, fill = "CC")) +
  
  # Stratify by screening participation
  facet_wrap(~ strata) +
  
  # Labels
  labs(
    title = "",
    subtitle = "",
    y = "Breast cancer incidence by mode of detection (%)",
    x = "Time from first screening invitation (years)",
    fill = "Detection mode"
  ) +
  # Optional: use your custom theme
  theme_mine
