####--------------------------------------------------
#### Title: Generate Figure 2
#####Participation trends (overall vs. first screen) and 
#####cumulative screening rounds by initial participation status
#### Author: Ziyan Ma
#### Last Updated: 2025-06-06
####--------------------------------------------------

# Clear environment
rm(list = ls())

# Load libraries
library(tidyverse)
library(ggplot2)
library(ggsci)
library(patchwork)
library(boot)

# Set working directory
setwd("")

# Load custom ggplot theme
source("ggplot_theme.R")

# Load data
load("mammography_screening_first_screening_defined.RData")
load("mammography_screening.RData")

# ------------------------------------------------------------------------------
# Figure 2a: Participation trend over time (first screening vs. any screening)
# ------------------------------------------------------------------------------

# First screening participation rate per year
first_part <- mammography_screening_first_screening_defined %>%
  group_by(RequestYear_1) %>%
  summarise(participation_rate_1 = mean(participate_1)) %>%
  ungroup() %>%
  rename(RequestYear = RequestYear_1)

# Overall screening participation rate (among those in the first screening cohort)
overall_part <- mammography_screening %>%
  filter(LopNr %in% mammography_screening_first_screening_defined$LopNr) %>%
  mutate(RequestYear = as.numeric(substr(RequestDate, 1, 4))) %>%
  group_by(RequestYear) %>%
  summarise(participation_rate = mean(participate)) %>%
  ungroup()

# Combine and convert to percentages
dat <- first_part %>%
  left_join(overall_part, by = "RequestYear") %>%
  mutate(
    participation_rate_1 = participation_rate_1 * 100,
    participation_rate = participation_rate * 100
  )

# Plot Figure 2a
p_2a <- ggplot(data = dat, aes(x = RequestYear)) +
  geom_line(aes(y = participation_rate), size = 1, color = "darkgrey") +
  geom_point(aes(y = participation_rate_1), size = 0.8, color = "black") +
  labs(x = "Year", y = "Screening participation rate (%)") +
  ylim(0, 100) +
  scale_x_continuous(limits = c(1991, 2018), breaks = c(1991, 1995, 2000, 2005, 2010, 2015)) +
  theme_minimal()

# ------------------------------------------------------------------------------
# Figure 2b: Cumulative number of screening rounds (by initial participation)
# ------------------------------------------------------------------------------

# Reshape screening data to long format and calculate cumulative participation
df <- mammography_screening_first_screening_defined %>%
  select(LopNr, paste0("participate_", 1:10)) %>%
  mutate(GroupScreen = participate_1) %>%
  pivot_longer(
    cols = starts_with("participate_"),
    names_to = "round",
    values_to = "participate"
  ) %>%
  drop_na(participate) %>%
  mutate(
    round = gsub("participate_", "", round),
    round = as.integer(round)
  ) %>%
  group_by(LopNr) %>%
  arrange(LopNr, round) %>%
  mutate(cum_participate = cumsum(participate)) %>%
  ungroup()

# Bootstrap function
bootstrap_mean <- function(data, indices) {
  mean(data[indices])
}

# Calculate mean cumulative participation with 95% bootstrap CI
cum_participation_bootstrap <- df %>%
  group_by(round, GroupScreen) %>%
  summarise(
    n_women = n(),
    mean_cum_participate = mean(cum_participate),
    bootstrap_ci = list({
      if (n() >= 10) {
        boot_result <- boot(cum_participate, bootstrap_mean, R = 2000)
        ci_result <- boot.ci(boot_result, type = "perc")
        if (!is.null(ci_result)) ci_result$percent[4:5] else c(NA, NA)
      } else {
        c(NA, NA)
      }
    }),
    .groups = "drop"
  ) %>%
  mutate(
    ci_lower_boot = map_dbl(bootstrap_ci, ~ ifelse(is.na(.x[1]), NA, .x[1])),
    ci_upper_boot = map_dbl(bootstrap_ci, ~ ifelse(is.na(.x[2]), NA, .x[2])),
    round = factor(round, levels = 1:10),
    GroupScreen = ifelse(GroupScreen == 0, "Non-participant", "Participant")
  ) %>%
  select(-bootstrap_ci) %>%
  mutate(across(where(is.numeric), round, 3))

# Plot Figure 2b
p_2b <- ggplot(data = cum_participation_bootstrap) +
  geom_step(aes(x = round, y = mean_cum_participate, group = GroupScreen, colour = GroupScreen), size = 1) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 1)) +
  labs(x = "Screening round", y = "Cumulative participated\nscreening rounds") +
  theme_minimal()

# Optionally combine plots
# combined_plot <- p_2a / p_2b + plot_layout(heights = c(1, 1.3))
# print(combined_plot)

library(patchwork)
combined_plot <- p_2a / p_2b

