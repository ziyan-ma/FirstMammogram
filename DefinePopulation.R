#####################################################
#### Title: Define Study Population
#### Author: Ziyan Ma
#### Last Updated: 2025-06-08
#####################################################
# Note: "LopNr" is the unique identifier used consistently across datasets.

# Environment setup
rm(list = ls())
library(tidyverse)
library(lubridate)

setwd("")  #Set working directory

####################
### Dataset 1: Mammography Screening
####################

load("mammography_screening.RData")
paste0("Remaining population size: ", length(unique(mammography_screening$LopNr)))

# Step 1: Basic cleaning
# - Keep relevant columns
# - Remove duplicates
# - Remove rows with missing screening invitation dates
# - Keep only invitations from 1991 onwards
# - Remove duplicate screening dates per person

mammography_screening <- mammography_screening %>%
  dplyr::select(LopNr, RequestDate, Screening_ExaminationCancelCode, participate, Birthdate, age_at_invitation) %>%
  mutate(RequestYear = as.numeric(substr(RequestDate, 1, 4))) %>%
  distinct() %>%
  drop_na(RequestDate) %>%
  filter(RequestYear >= 1991) %>%
  group_by(LopNr) %>%
  distinct(RequestDate, .keep_all = TRUE) %>%
  ungroup()
paste0("Remaining population size: ", length(unique(mammography_screening$LopNr)))

# Step 2: Define valid screening intervals
# - Exclude screening intervals <1 year apart
# - Exclude all if interval >5 years
# - Retain only 1st to 10th screening rounds
# - Define follow-up end date as 2 years after last invitation

mammography_screening <- mammography_screening %>%
  group_by(LopNr) %>%
  arrange(RequestDate) %>%
  mutate(
    year_of_screening_interval = ifelse(RequestDate == min(RequestDate), 0.5,
                                        RequestYear - lag(RequestYear)),
    longer_than_5yr_interval = any(year_of_screening_interval > 5),
    keep = case_when(
      is.na(year_of_screening_interval) ~ TRUE,
      year_of_screening_interval == 0 ~ FALSE,
      !longer_than_5yr_interval ~ TRUE,
      RequestYear >= min(RequestYear[year_of_screening_interval > 5]) ~ FALSE,
      TRUE ~ TRUE
    )
  ) %>%
  filter(keep == TRUE) %>%
  mutate(
    Screening_Round_redefine = 1:n(),
    last_RequestDate = max(RequestDate),
    last_RequestDate_plus_2y = last_RequestDate %m+% years(2)
  ) %>%
  filter(Screening_Round_redefine %in% 1:10) %>%
  ungroup()

paste0("Remaining population size: ", length(unique(mammography_screening$LopNr)))

# Step 3: Define cohort based on age at first invitation
# - Women aged 49–51 invited between 1991–2010 (start from 50)
# - Women aged 39–41 invited from 2005 onwards (start from 40)

mammography_screening <- mammography_screening %>%
  group_by(LopNr) %>%
  mutate(first_invitation_age = min(age_at_invitation)) %>%
  ungroup() %>%
  filter(
    (RequestYear %in% 1991:2010 & first_invitation_age %in% 49:51) |
      (RequestYear >= 2005 & first_invitation_age %in% 39:41)
  ) %>%
  mutate(first_invitation_age_cate = case_when(
    first_invitation_age %in% 49:51 ~ "start from 50",
    first_invitation_age %in% 39:41 ~ "start from 40",
    TRUE ~ NA_character_
  ))
paste0("Remaining population size: ", length(unique(mammography_screening$LopNr)))

# Step 4: Flag women with more than one screening
mammography_screening <- mammography_screening %>%
  group_by(LopNr) %>%
  mutate(Max_Screen_More_Than_One = max(Screening_Round_redefine) > 1) %>%
  ungroup() %>%
  dplyr::select(LopNr, RequestDate, RequestYear, participate, Birthdate, age_at_invitation,
                Screening_Round_redefine, last_RequestDate, last_RequestDate_plus_2y,
                Screening_ExaminationCancelCode, Max_Screen_More_Than_One)

paste0("Remaining population size: ", length(unique(mammography_screening$LopNr)))

####################
### Dataset 2: Mode of Detection of Breast Cancer
####################

load("mode_of_detection.RData")

mode_of_detection <- mode_of_detection %>%
  arrange(diagnosis_date) %>%
  filter(!duplicated(LopNr)) %>%
  dplyr::select(LopNr, mode_of_detection, diagnosis_date, RequestDate) %>%
  mutate(mode_of_detection_clinical = (mode_of_detection == "CC_non_participate")) %>%
  rename(RequestDate_bfBC = RequestDate)

# Merge with screening data
mammography_screening_first_screening_defined <- mammography_screening %>%
  left_join(mode_of_detection)

# Remove screenings after breast cancer diagnosis
mammography_screening_first_screening_defined <- mammography_screening_first_screening_defined %>%
  filter(is.na(diagnosis_date) | diagnosis_date > RequestDate)

# Reshape: wide format by screening round
mammography_screening_first_screening_defined <- mammography_screening_first_screening_defined %>%
  pivot_wider(names_from = Screening_Round_redefine,
              values_from = c(participate, RequestDate, RequestYear, age_at_invitation, Screening_ExaminationCancelCode)) %>%
  filter(participate_1 %in% c(0, 1))

paste0("Remaining population size: ", length(unique(mammography_screening$LopNr)))

####################
### Dataset 3: Swedish Cancer Register (Swecan)
####################

load("swecan.RData")

swecan <- swecan %>%
  arrange(DiagnoseDate) %>%
  distinct(LopNr, .keep_all = TRUE) %>%
  dplyr::select(LopNr, DiagnoseDate, BC_diagnosis) %>%
  rename(DiagnoseDate_swecan_1st_cancer_ever = DiagnoseDate)

# Merge and exclude women with cancer before first screening
mammography_screening_first_screening_defined <- mammography_screening_first_screening_defined %>%
  left_join(swecan) %>%
  filter(is.na(DiagnoseDate_swecan_1st_cancer_ever) | DiagnoseDate_swecan_1st_cancer_ever > RequestDate_1)

paste0("Remaining population size: ", length(unique(mammography_screening_first_screening_defined$LopNr)))

#####################################################
### Following sections add background variables (not for cohort filtering)
#####################################################

####################
### Dataset 4: LISA Education & Income
####################

load("lisa_index_edu_income.RData")

lisa_index_edu_income <- lisa_index_edu_income %>%
  inner_join(mammography_screening_first_screening_defined[, c("LopNr", "RequestYear_1")]) %>%
  group_by(LopNr) %>%
  filter(abs(Year_lisa - RequestYear_1) == min(abs(Year_lisa - RequestYear_1))) %>%
  ungroup() %>%
  arrange(Year_lisa) %>%
  filter(!duplicated(LopNr)) %>%
  dplyr::select(LopNr, Education_regroup, Income_individual_standardized_category, Income_family_standardized_category) %>%
  rename(
    Education_regroup_at_screening_1 = Education_regroup,
    Income_individual_standardized_category_at_screening_1 = Income_individual_standardized_category,
    Income_family_standardized_category_at_screening_1 = Income_family_standardized_category
  )

# Merge
mammography_screening_first_screening_defined <- mammography_screening_first_screening_defined %>%
  left_join(lisa_index_edu_income)

####################
### Dataset 5: Key Variables (Demographics)
####################

load("key.RData")

key <- key %>%
  filter(screening == 1) %>%
  mutate(have_child = !is.na(LopNr_Child1)) %>%
  dplyr::select(LopNr, born_Sweden, Birthdate, deaddate, have_child)

# Merge
mammography_screening_first_screening_defined <- mammography_screening_first_screening_defined %>%
  left_join(key)

####################
### Dataset 6: Migration History
####################

load("migration.RData")

# Merge
mammography_screening_first_screening_defined <- mammography_screening_first_screening_defined %>%
  left_join(migration[, c("LopNr", "EmigrationDate")])