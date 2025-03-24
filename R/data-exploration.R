#!/usr/bin/env Rscript

## Setup ----------------------------------------------------------------------------

# Specify options for saving the plots to files
figures_dir = "results/figures/"
tables_dir = "results/tables/"

# Load required packages. 
library(tidyverse)
library(patchwork)

## Load data and Prepare for Analysis -----------------------------------

# Load the data set. Wstratum refers to the strata that determine the
# probability of being sampled for measuring S. 
df = read.csv("data/processed_data.csv") 


df = df %>%
  mutate(trial = factor(trial.lbl)) 

df = df %>%
  mutate(BMI_factor = ifelse(
    BMI_underweight,
    "< 18.5",
    ifelse(
      BMI_normal,
      "[18.5, 25)",
      ifelse(BMI_overweight, "[25, 30)", ">= 30")
    )
  ))

# For computing the limits of detection, we need to use the original data. In
# the processed data, titer values have been modified to reflect the limits of
# detection.
df_original = read.csv("data/CrossProtocolData_moderna_tgt_FINAL.csv")

# Add label for the trials.
df_original = df_original %>%
  mutate(trial.lbl = sapply(
    X = trial,
    FUN = function(x) {
      switch (
        as.character(x),
        "1" = "Moderna",
        "2" = "AstraZeneca",
        "3" = "Novavax",
        "4" = "J&J (Brazil)",
        "5" = "J&J (Colombia)",
        "6" = "J&J (S. America)",
        "7" = "J&J (S. Africa)",
        "8" = "J&J (USA)"
      )
    }
  ))

lod_original = df_original %>%
  filter(A != 0, Delta == 1) %>%
  group_by(protocol) %>%
  summarise(
    "LLOQ IgG Spike" = 10 ** min(bindSpike),
    "ULOQ IgG Spike" = 10 ** max(bindSpike),
    "LLOQ nAb ID50" = 10 ** min(pseudoneutid50),
    "ULOQ nAb ID50" = 10 ** max(pseudoneutid50)
  )
  
# Figures -------------------------------------------------------------

## Distribution of Risk Score  ----------------------------------------

df %>%
  ggplot(aes(x = risk_score, y = trial)) +
  geom_boxplot(position = "identity", color = "black", fill = "gray") +
  ylab("Trial") +
  xlab("Risk Score")

ggsave(
  filename = "distribution-risk-score.pdf",
  device = "pdf",
  path = figures_dir,
  width = double_width,
  height = double_height, 
  units = unit
)

## Distribution of Surrogates -----------------------------------------

# We plot the distribution of the Ab titers before truncation using the LLOD and
# ULODs. We also add the cut-offs based on LLOD and ULOD to the plots as
# vertical dashed lines. In the data used in the final analyses, titers below
# and above these cut-offs are truncated to the respective cut-off values.

# Compute LODs for each trial and type of surrogate.
LOD <- df_original %>% group_by(protocol) %>% dplyr::summarise(
  LLOD.Spike = min(bindSpike, na.rm = T),
  LLOD.Neut = min(pseudoneutid50, na.rm = T)
)

# Maximum LOQs across all trials. 
maxLLODSpike <- max(LOD$LLOD.Spike)
maxLLODNeut <- max(LOD$LLOD.Neut)

# Plot for binding Ab.
bAb_plot = df_original %>%
  filter(A != 0, Delta == 1) %>%
  ggplot(aes(x = bindSpike, y = trial.lbl)) +
  geom_boxplot(position = "identity",
               color = "black",
               fill = "gray") +
  geom_vline(
    xintercept = c(maxLLODSpike),
    color = "gray",
    linetype = "dashed"
  ) +
  ylab("Trial") +
  xlab(expression(paste("Spike Protein IgG (", log[10], "BAU/ml)"))) 

# Plot for neutralizing Ab.
nAb_plot = df_original %>%
  filter(A != 0, Delta == 1) %>%
  ggplot(aes(x = pseudoneutid50, y = trial.lbl)) +
  geom_boxplot(position = "identity",
               color = "black",
               fill = "gray") +
  geom_vline(
    xintercept = c(maxLLODNeut),
    color = "gray",
    linetype = "dashed"
  ) +
  ylab("Trial") +
  xlab(expression(paste("Neutralizing Antibody (", log[10], "ID50)"))) 

bAb_plot / nAb_plot

ggsave(
  filename = "distribution-titers.pdf",
  device = "pdf",
  path = figures_dir,
  width = double_width,
  height = double_height, 
  units = unit
)


# Tables ------------------------------------------------------------------

## Case-Cohort Sampling and Events ----------------------------------------

# Compute the number of infections by trial and treatment arm.
table_infections = df %>%
  mutate(
    treatment = ifelse(
      vax,
      "Vaccine (number of infected subjects/total number of subjects)",
      "Placebo (number of infected subjects/total number of subjects)"
    )
  ) %>%
  group_by(trial, treatment) %>%
  summarise(
    n = n(),
    n_infected = sum(Y),
    results = paste0(n_infected, "/", n)
  ) %>%
  select(-n, -n_infected) %>%
  pivot_wider(names_from = "treatment", values_from = c("results"))

# Compute the number of subjects sampled in the case-cohort sampling by trial
# and infection status. We only select the vaccine recipients here because the
# titer is not measured for placebo recipients.
table_case_cohort = df %>%
  filter(vax == 1) %>%
  mutate(
    infected = ifelse(
      Y,
      "Cases (number of titers measures/total number of subjects)",
      "Non-Cases (number of titers measures/total number of subjects)"
    )
  ) %>%
  group_by(trial, infected) %>%
  summarise(
    n = n(),
    n_Delta = sum(Delta),
    results = paste0(n_Delta, "/", n)
  ) %>%
  select(-n, -n_Delta) %>%
  pivot_wider(names_from = "infected", values_from = "results")


# Join to tables together and save.
write.csv(
  table_infections %>%
    left_join(table_case_cohort),
  file = paste0(tables_dir, "infections-case-cohort.csv")
)

## Limits of Detection ----------------------------------------------------



lod_original %>%
  write.csv(
    file = paste0(tables_dir, "limits-of-detection.csv")
  )

         