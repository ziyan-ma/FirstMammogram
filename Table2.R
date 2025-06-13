####--------------------------------------------------
#### Title: Generate Table 2
#####Logistic regression on tumor characteristics, 
#####by first screening participation
#### Author: Ziyan Ma
#### Last Updated: 2025-05-30
####--------------------------------------------------
# Note: "LopNr" is the unique identifier used consistently across datasets.

# Clear environment and load packages
rm(list = ls())
library(tidyverse)
library(tidycmprsk)
library(eoffice)
library(nnet)

# Set working directory
setwd("")

# Source custom theme
source("ggplot_theme.R")

# ------------------------------
# Load and filter screening data
# ------------------------------
load("mammography_screening_first_screening_defined.RData")

screening_data <- mammography_screening_first_screening_defined %>%
  filter(!is.na(mode_of_detection)) %>%
  filter(as.numeric(substr(diagnosis_date, 1, 4)) - as.numeric(substr(last_RequestDate, 1, 4)) <= 3)

# ------------------------------
# Load tumor characteristic data
# ------------------------------
load("nkbc_rcc_tumorcha_1stdiag.RData")

tumor_data <- nkbc_rcc_tumorcha_1stdiag %>%
  select(LopNr, DiagnoseDate, insitu, T, N, M, TNM_stage) %>%
  mutate(
    DiagnoseYear = as.numeric(substr(DiagnoseDate, 1, 4))
  )

# ------------------------------
# Merge and create derived variables
# ------------------------------
dat <- screening_data %>%
  filter(!duplicated(LopNr)) %>%
  inner_join(tumor_data, by = "LopNr") %>%
  mutate(
    age_at_diagnosis = floor(as.numeric(difftime(DiagnoseDate, Birthdate, units = "days") / 365.25)),
    invasive = case_when(insitu == 1 ~ 0, insitu == 0 ~ 1, TRUE ~ NA_real_)
  )

# Convert to factors
dat <- dat %>%
  mutate(
    T = factor(T, levels = 0:2),
    N = factor(N, levels = 0:1),
    invasive = factor(invasive, levels = 0:1),
    M = factor(M, levels = 0:1),
    TNM_stage = factor(TNM_stage, levels = 1:4),
    participate_1 = factor(participate_1, levels = c(1, 0)),
  )

# ------------------------------
# Define analysis function
# ------------------------------
tumor_char_fun <- function(outcome, adjust = FALSE) {
  formula_base <- paste0(outcome, " ~ participate_1 + DiagnoseYear + age_at_diagnosis")
  formula_full <- paste0(formula_base, " + ",
                         paste(common_factors, collapse = " + "))
  formula_use <- ifelse(adjust, formula_full, formula_base)
  
  model <- if (length(levels(dat_sub[[outcome]])) >= 3) {
    multinom(as.formula(formula_use), data = dat_sub)
  } else {
    glm(as.formula(formula_use), data = dat_sub, family = binomial())
  }
  
  if (length(levels(dat_sub[[outcome]])) >= 3) {
    coefs <- summary(model)$coefficients
    se <- summary(model)$standard.errors
    est <- paste(coefs[, "participate_10"], collapse = "+")
    stderr <- paste(se[, "participate_10"], collapse = "+")
  } else {
    smry <- summary(model)$coefficients
    est <- smry["participate_10", "Estimate"]
    stderr <- smry["participate_10", "Std. Error"]
  }
  
  return(list(est, stderr))
}

# ------------------------------
# Handle missing for common covariates
# ------------------------------
common_factors <- c("Education_regroup_at_screening_1", "Income_individual_standardized_category_at_screening_1",
                    "born_Sweden", "n_child_cate", "first_child_cate", "alcohol_related_diasease",
                    "obesity_related_diasease", "BC_FH", "AllC_FH", "CCIw_Regroup", "Marital")

dat[common_factors] <- lapply(dat[common_factors], function(x) ifelse(is.na(x), "missing", x))

# ------------------------------
# Run models for tumor characteristics
# ------------------------------
tumor_outcomes <- c("T", "invasive", "M", "TNM_stage", "N")

glm_output <- data.frame()

for (i in tumor_outcomes) {
  for (j in c(FALSE, TRUE)) {
    if (i == "invasive") {
      dat_sub <- dat
    } else {
      dat_sub <- dat %>%
        filter(invasive != 0)
    }
    temp <- data.frame(
      Outcome = NA, 
      beta = NA, 
      se = NA, 
      adjust_for_common_factors = NA)
    temp$Outcome <- i
    temp$beta <- tumor_char_fun(i, j)[[1]]
    temp$se <- tumor_char_fun(i, j)[[2]]
    temp$adjust_for_common_factors <- j
    glm_output <- rbind(glm_output, temp)
  }}

# ------------------------------
# Format regression results
# ------------------------------
glm_output <- glm_output %>%
  separate_wider_delim(beta, delim = "+", names = paste0("beta_", 1:3), 
                       too_few = "debug",
                       too_many = "debug") %>%
  separate_wider_delim(se, delim = "+", names = paste0("se_", 1:3), 
                       too_few = "debug",
                       too_many = "debug") %>%
  rename(name = Outcome) %>%
  dplyr::select(name, beta_1, beta_2, beta_3, se_1, se_2, se_3, adjust_for_common_factors)

glm_output <- glm_output %>%
  pivot_longer(cols=-c("name", "adjust_for_common_factors"), names_pattern = "(.*)_(.)", names_to = c(".value", "level")) %>%
  filter(!is.na(beta)) %>%
  mutate(adjust_for_common_factors = factor(adjust_for_common_factors, levels = c(TRUE, FALSE)))


glm_output <- glm_output %>%
  mutate(beta = as.numeric(beta)) %>%
  mutate(se = as.numeric(se)) %>%
  mutate(OR = exp(beta)) %>%
  mutate(OR_lb = exp(beta - 1.96*se)) %>%
  mutate(OR_ub = exp(beta + 1.96*se)) %>%
  mutate(OR_ub_lb = paste0(round(OR, digits = 2), "(", round(OR_lb, digits = 2), ",", round(OR_ub, digits = 2), ")"))
