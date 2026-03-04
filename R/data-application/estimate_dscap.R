#!/usr/bin/env Rscript

# load data and setup ---------------------------------------------------------------
t1 = Sys.time()

## Setup ----------------------------------------------------------------------------

set.seed(2)

# Extract arguments for analysis. 
args = commandArgs(trailingOnly=TRUE)
antibody_type <- args[1]
target <- args[2]
vars = args[3]
mnum <- as.integer(args[4])
n_boot <- as.integer(args[5])
n_perm <- as.integer(args[6])
estimate_weights = as.integer(args[7])
truncation = as.integer(args[8])
data_location = args[9]

if (estimate_weights) {
  outfile_wts = "estwts"
} else {
  outfile_wts = "kwnwts"
}

population = "full"
if (truncation) population = "truncated"

# Load required packages. 
library(tidyverse)

# Load functions that will be used later on. 
source("R/helper-functions/treatment-effect-estimators.R")
source("R/helper-functions/DSCAP-estimators.R")
source("R/helper-functions/permutation-LRT.R")

# Define formulas for the outcome models for (i) clinical outcome Y and (ii)
# surrogate outcome S. The same linear predictors are used for both models. 
formula_Y = as.formula(paste0("Y ~ ", vars))
formula_S = as.formula(paste0("S ~ ", vars))
formula_T = as.formula(paste0("trial ~ ", vars))
formula_CC = as.formula(paste0("Delta ~ ", "CC_stratum"))

## Load data and Prepare for Analysis -----------------------------------

# Load the data set. Wstratum refers to the strata that determine the
# probability of being sampled for measuring S. 
df = read.csv(data_location)

target_trial = switch(
  target,
  jjsa = "J&J (S. America)",
  AZ = "AstraZeneca",
  jjbr = "J&J (Brazil)",
  jjusa  = "J&J (USA)",
  jjcol = "J&J (Colombia)",
  nov = "Novavax"
)


df = df %>%
  mutate(trial = factor(trial.lbl)) 

# Depending on the value of `mnum`, we drop some trials from the analysis.
if (mnum == 1) {
  # No trials are dropped.
  df <- df
} else if (mnum == 2) {
  # Drop two trials.
  df <- df %>% filter(!(
    trial.lbl %in% c("J&J (Brazil)", "J&J (Colombia)", "Novavax")
  ))  %>%
    mutate(trial = droplevels(trial)) #drop 
}

# If target trial is not in the set of trials in the selected data, raise an
# error.
if (!(target_trial %in% levels(df$trial))) {
  stop("Target trial not in set of selected trials for analysis.")
}

# Truncate the risk scores in the target trial is asked. 
if (truncation) {
  # Compute the minimum of the trial-specific maxima of the risk score. 
  min_max_risk_score = df %>%
    group_by(trial) %>%
    summarise(max = max(risk_score)) %>%
    ungroup() %>%
    group_by(1) %>%
    summarise(min = min(max)) %>%
    pull(min)
  # Compute the maximum of the trial-specific minima of the risk score.
  max_min_risk_score = df %>%
    group_by(trial) %>%
    summarise(min = min(risk_score)) %>%
    ungroup() %>%
    group_by(1) %>%
    summarise(max = max(min)) %>%
    pull(max)
  
  # We add an epsilon to the interval limits. Otherwise, we will only select few
  # subjects in the target trial.
  min_max_risk_score = min_max_risk_score
  max_min_risk_score = max_min_risk_score
  
  print(paste0("Subjects in the target trial with risk scores outside of [", max_min_risk_score, ", ", min_max_risk_score, "] are dropped."))
  
  # Exclude subjects from the target trial with risk scores outside the above
  # determined bounds.
  df = df %>%
    filter((trial != target_trial) | (risk_score >= max_min_risk_score & min_max_risk_score >= risk_score))
}

# Variable defining the population (i.e., trial + placebo group). Subjects from
# the placebo group across different trials are assumed to be exchangeable; they
# all get A = 0.
df$A <- ifelse(df$A == 0, 0, as.integer(df$trial))

# Compute interaction covariate between the risk score and age. 
df$riskxage <- df$risk_score*df$age.geq.65

# Define S outcome based on inputs
if (antibody_type == "spike") {
  df = df %>%
    dplyr::select(-pseudoneutid50) %>%
    rename(S = bindSpike) %>%
    droplevels()
} else if (antibody_type == "neut") {
  df = df %>%
    dplyr::select(-bindSpike) %>%
    rename(S= pseudoneutid50) %>%
    droplevels()
}

# Recode Delta as equal to one in the placebo groups because the immune response
# is known to be the lower bound for this subjects. 
df <- df %>% mutate(Delta = ifelse(A == 0, 1, Delta))


# Data Analysis -----------------------------------------------------------

## Estimation -------------------------------------------------------------

# Use the doubly robust, IPW, standardized, and naive estimators to estimate the
# treatment effects and corresponding measures of association.

results_DR = inference_DSCAP(
  data = df,
  formula_S = formula_S,
  formula_Y = formula_Y,
  formula_T = formula_T,
  formula_CC = formula_CC,
  target_trial = target_trial,
  estimate_weights = estimate_weights,
  alpha = 0.05,
  B = n_boot,
  type = "doubly robust",
  treatment_var = "A",
  trial_var = "trial",
  corrected_target_trial = FALSE
)

results_ipw = inference_DSCAP(
  data = df,
  formula_S = formula_S,
  formula_Y = formula_Y,
  formula_T = formula_T,
  formula_CC = formula_CC,
  target_trial = target_trial,
  estimate_weights = estimate_weights,
  alpha = 0.05,
  B = n_boot,
  type = "ipw",
  treatment_var = "A",
  trial_var = "trial",
  corrected_target_trial = FALSE
)

results_st = inference_DSCAP(
  data = df,
  formula_S = formula_S,
  formula_Y = formula_Y,
  formula_T = formula_T,
  formula_CC = formula_CC,
  target_trial = target_trial,
  estimate_weights = estimate_weights,
  alpha = 0.05,
  B = n_boot,
  type = "standardized",
  treatment_var = "A",
  trial_var = "trial",
  corrected_target_trial = FALSE
)

results_naive = inference_DSCAP(
  data = df,
  formula_S = formula_S,
  formula_Y = formula_Y,
  formula_T = formula_T,
  formula_CC = formula_CC,
  target_trial = target_trial,
  estimate_weights = estimate_weights,
  alpha = 0.05,
  B = n_boot,
  type = "naive",
  treatment_var = "A",
  trial_var = "trial",
  corrected_target_trial = FALSE
)

# Determine to which trial each treatment group corresponded.
treatment_trial_mapping <- df %>%
  filter(A != 0) %>%
  group_by(A) %>%
  summarise(trial = unique(trial)) %>%
  ungroup()

results_DR <- results_DR %>% 
  left_join(treatment_trial_mapping, by = c("treatment" = "A"))
results_ipw <- results_ipw %>%
  left_join(treatment_trial_mapping, by = c("treatment" = "A"))
results_st <- results_st %>%
  left_join(treatment_trial_mapping, by = c("treatment" = "A"))
results_naive <- results_naive %>%
  left_join(treatment_trial_mapping, by = c("treatment" = "A"))

# Merge the tibbles with estimates and save as an rds file.
all_results_tbl = bind_rows(
  results_DR %>% mutate(type = "doubly robust"),
  results_ipw %>% mutate(type = "ipw"),
  results_st %>% mutate(type = "standardized"),
  results_naive %>% mutate(type = "naive")
)
saveRDS(
  all_results_tbl,
  file = paste0(
    "results/data-application/raw-results/DSCAP_estimates_",
    population,
    "_",
    target,
    "_",
    outfile_wts,
    "_",
    mnum,
    "_",
    antibody_type,
    ".rds"
  )
)


