####--------------------------------------------------
#### Title: Generate Table 1
#### Author: Ziyan Ma
#### Last Updated: 2025-05-20
####--------------------------------------------------
# Note: "LopNr" is the unique identifier used consistently across datasets.

# Clear environment
rm(list = ls())

# Load libraries
library(tidyverse)
library(table1)

# Set working directory
setwd("")

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

load("mammography_screening_first_screening_defined.RData")
load("number_of_children_first_child_age.RData")
load("NPR_inpatient_outpatient_index_first_screening_alcohol.RData")
load("NPR_inpatient_outpatient_index_first_screening_obesity.RData")
load("FH_bf_1stScreen.RData")
load("NPR_Charlson_Comorbidity_Index.RData")
load("lisa_index_marital.RData")

# ------------------------------------------------------------------------------
# Marital status: Select closest year to screening for each woman
# ------------------------------------------------------------------------------

lisa_index_marital <- lisa_index_marital %>%
  inner_join(mammography_screening_first_screening_defined %>% select(LopNr, RequestYear_1)) %>%
  mutate(DiffYear = as.numeric(Year_lisa) - RequestYear_1) %>%
  group_by(LopNr) %>%
  filter(abs(DiffYear) == min(abs(DiffYear))) %>%
  ungroup() %>%
  distinct(LopNr, .keep_all = TRUE)

# ------------------------------------------------------------------------------
# Prepare Charlson comorbidity index
# ------------------------------------------------------------------------------

Charlson_Comorbidity_Index <- Charlson_Comorbidity_Index %>%
  rename(LopNr = group) %>%
  select(LopNr, CCIw_Regroup)

# ------------------------------------------------------------------------------
# Merge all covariates into main dataset
# ------------------------------------------------------------------------------

mammography_screening_first_screening_defined <- mammography_screening_first_screening_defined %>%
  left_join(child) %>%
  left_join(NPR_inpatient_outpatient_index_first_screening_alcohol) %>%
  left_join(NPR_inpatient_outpatient_index_first_screening_obesity) %>%
  left_join(FH_bf_1stScreen) %>%
  left_join(Charlson_Comorbidity_Index) %>%
  left_join(lisa_index_marital[, c("LopNr", "Marital")])

# ------------------------------------------------------------------------------
# Generate variables for Table 1
# ------------------------------------------------------------------------------

mammography_screening_first_screening_defined <- mammography_screening_first_screening_defined %>%
  mutate(
    alcohol_related_diasease = replace_na(alcohol_related_diasease, FALSE),
    obesity_related_diasease = replace_na(obesity_related_diasease, FALSE),
    CCIw_Regroup = replace_na(CCIw_Regroup, 0)
  )

# ------------------------------------------------------------------------------
# Define labels and formats for Table 1
# ------------------------------------------------------------------------------

df <- mammography_screening_first_screening_defined %>%
  mutate(
    born_Sweden = factor(born_Sweden, levels = c(1, 0), labels = c("Yes", "No")),
    participate_1 = factor(participate_1, levels = c(1, 0), labels = c("Participants", "Non-participants")),
    Education_regroup_at_screening_1 = factor(Education_regroup_at_screening_1, labels = c("<10", "10–12", ">12")),
    Income_individual_standardized_category_at_screening_1 = factor(Income_individual_standardized_category_at_screening_1, labels = c("Q1", "Q2", "Q3", "Q4")),
    obesity_related_diasease = factor(obesity_related_diasease, levels = c(TRUE, FALSE), labels = c("Yes", "No")),
    alcohol_related_diasease = factor(alcohol_related_diasease, levels = c(TRUE, FALSE), labels = c("Yes", "No")),
    first_child_cate = factor(first_child_cate, levels = c("<=30", ">30")),
    n_child_cate = factor(n_child_cate, levels = c("None", "1-2", ">=3")),
    RequestYear_1_cate = factor(RequestYear_1_cate, levels = c("1991-2000", "2001-2010", "2011-2020")),
    AllC_FH = factor(AllC_FH, levels = c("Parent", "Sibling", "Parent_Sibling", "None"), labels = c("Parent", "Sibling", "Parent and sibling", "None")),
    BC_FH = factor(BC_FH, levels = c("Mother", "Sister", "Mother_Sister", "None"), labels = c("Parent", "Sibling", "Mother and sibling", "None")),
    CCIw_Regroup = factor(CCIw_Regroup, levels = c("0", "1-2", ">=3")),
    Marital = factor(Marital, levels = c("Single", "Married_Registered", "Divorced_Separated"), labels = c("Single", "Married", "Divorced"))
  )

label(df$first_invitation_age_cate) <- "Age of 1st invitation"
label(df$Education_regroup_at_screening_1) <- "Educational attainment"
label(df$Income_individual_standardized_category_at_screening_1) <- "Income level"
label(df$born_Sweden) <- "Born in Sweden"
label(df$have_child) <- "Have children"
label(df$participate_1) <- "Participation in 1st screening"
label(df$obesity_related_diasease) <- "Obesity-related diseases"
label(df$alcohol_related_diasease) <- "Alcohol-related diseases"
label(df$n_child_cate) <- "Parity"
label(df$first_child_cate) <- "Age at first child"
label(df$RequestYear_1_cate) <- "Year of 1st invitation"
label(df$AllC_FH) <- "Cancer family history"
label(df$BC_FH) <- "Breast cancer family history"
label(df$CCIw_Regroup) <- "Charlson comorbidity index"
label(df$Marital) <- "Marital status"

units(df$Education_regroup_at_screening_1) <- "years"

# ------------------------------------------------------------------------------
# Generate and display Table 1
# ------------------------------------------------------------------------------

table1(
  ~ first_invitation_age_cate + RequestYear_1_cate + Education_regroup_at_screening_1 +
    Income_individual_standardized_category_at_screening_1 + born_Sweden + Marital +
    n_child_cate + first_child_cate + alcohol_related_diasease + obesity_related_diasease +
    CCIw_Regroup + BC_FH + AllC_FH | participate_1,
  data = df,
  overall = TRUE,
  render.categorical = "FREQ (PCTnoNA%)"
)