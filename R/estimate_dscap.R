#!/usr/bin/env Rscript

# load data and setup ---------------------------------------------------------------
t1 = Sys.time()

## Setup ----------------------------------------------------------------------------

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
# A custom version of the geex package is used to avoid issues with large
# objects.
# install.packages("geex_1.1.1.tar.gz", repos = NULL)
library(geex)

# Load functions that will be used later on. 
source("R/helper-functions/estFUN.R")
source("R/helper-functions/runDSCAP.R")

# Define formulas for the outcome models for (i) clinical outcome Y and (ii)
# surrogate outcome S. The same linear predictors are used for both models. 
Ymod = as.formula(paste0("Y ~ ",vars))
Smod = as.formula(paste0("S ~ ",vars))

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

# Get number of trials
n_trials <- length(unique(df$trial))

# Add variable that defines the target population. 0: for target populations; 1
# otherwise.
df$R <- ifelse(df$trial == target_trial , 0, 1)


# Data Analysis -----------------------------------------------------------

## Estimation -------------------------------------------------------------

# Estimate all models.
RESULT <- RunDSCAP(
  data = df,
  formula_S = Smod,
  formula_Y = Ymod,
  target_trial = target_trial,
  sim = F, 
  estimate_weights = TRUE
)

# If some some strata have an estimated weight of Inf (i.e., zero probability of
# being sampled in the case-cohort sampling, they are excluded from the further
# analyses).
df = df %>%
  filter(!(CC_stratum %in% RESULT$excluded_CC_strata))
# Include estimated weights to data.
df = df %>%
  left_join(
    RESULT$weights_df
  )

## Conditional Independence Tests --------------------------------------------------------

# Fit null model (which does not contain trial as covariate).
glm_null <- glm(as.formula(paste0("Y ~ ", vars)), binomial, df %>%
                  filter(A == 0))
# Fit alternative model (which contains the interaction between trial and all
# covariates which were  already in the null model).
glm_w_trial <-update(glm_null, ~ . * trial)
# Perform likelihood ratio test. 
glm_lr_test = lmtest::lrtest(glm_null, glm_w_trial)

# Extract LRT statistic and corresponding p-value.
lrt_test_statistic <- glm_lr_test[2,4]
lrt_p_value <- glm_lr_test[2,5]

# Initialize data set in which the values of `trial` will be permuted. 
df_placebo_perm = df %>%
  filter(A == 0)
# Extract trial variable which will be permuted.
trial_vec = df_placebo_perm$trial

lrt_test_statistic_permutations = rep(NA, n_perm)
set.seed(1)
for(i in 1:n_perm){
  df_placebo_perm$trial <- sample(trial_vec, replace = FALSE)
  
  glm_w_trial <-update(glm_null, ~ . * trial, data = df_placebo_perm)
  lrt_test_statistic_permutations[i] <- lmtest::lrtest(glm_null, glm_w_trial)[2,4] #extract LRT stat
}

lrt_permuted_p_value <- mean(lrt_test_statistic_permutations > lrt_test_statistic)
lrt_out <- c("p-value" = lrt_p_value, "p-value (permutations)" = lrt_permuted_p_value, "LRT statistic" = lrt_test_statistic)

outfile_lrt = paste0(
  "results/raw-results/lrt_",
  antibody_type,
  "_",
  target,
  "_",
  population,
  "_M",
  mnum,
  "_",
  outfile_wts,
  ".csv"
)
write.csv(lrt_out, outfile_lrt, row.names=FALSE)


## Variance Estimation -----------------------------------------------------


outfile = paste0(
  "results/raw-results/",
  c("trt_effects_", "cor_est_"),
  antibody_type,
  "_",
  target,
  "_",
  population,
  "_M",
  mnum,
  "_",
  outfile_wts,
  ".csv"
)

write.csv(
  bind_rows(
    RESULT$standardized_trt_effects_df %>%
      mutate(type = "standardized"),
    RESULT$naive_trt_effects_df %>%
      mutate(type = "naive")
  ),
  outfile[1],
  row.names = FALSE
)
write.csv(
  bind_rows(
    RESULT$cor_standardized_df %>%
      mutate(type = "standardized"),
    RESULT$cor_naive_df %>%
      mutate(type = "naive")
  ),
  outfile[2],
  row.names = FALSE
)

## Bootstrap --------------------------------------------------------------

# Initialize dataframe in which the bootstrap replicates will be saved. Note the
# the results of a single bootstrap replication are saved over many rows. This
# facilities further processing.
bootstrap_replicates_df = data.frame(
  estimand = character(0),
  type = character(0),
  trial = character(0),
  estimate = numeric(0)
)

# Perform bootstrap by resampling within each trial.
set.seed(1)
for (i in 1:n_boot) {
  # Resample within each trial
  booti <- df %>%
    group_by(trial) %>%
    slice_sample(prop = 1.0, replace = TRUE)
  
  try(expr =  {
    temp <- RunDSCAP(
      data = booti,
      formula_S = Smod,
      formula_Y = Ymod,
      target_trial = target_trial,
      sim = F,
      estimate_weights = estimate_weights
    )
    bootstrap_replicates_df = bootstrap_replicates_df %>%
      bind_rows(
        tibble(
          estimand = c("cor_p", "cor_s", "beta"),
          type = "standardized",
          trial = NA,
          estimate = temp$cor_standardized_df[1, ] %>% as.numeric()
        ),
        tibble(
          estimand = c("cor_p", "cor_s", "beta"),
          type = "naive",
          trial = NA,
          estimate = temp$cor_naive_df[1, ] %>% as.numeric()
        ),
        tibble(
          estimand = "VE",
          type = "standardized",
          trial = temp$standardized_trt_effects_df %>% pull(trial),
          estimate = temp$standardized_trt_effects_df %>% pull(VE_est)
        ),
        tibble(
          estimand = "mean_diff_S",
          type = "standardized",
          trial = temp$standardized_trt_effects_df %>% pull(trial),
          estimate = temp$standardized_trt_effects_df %>% pull(mean_diff_S_est)
        ),
        tibble(
          estimand = "VE",
          type = "naive",
          trial = temp$naive_trt_effects_df %>% pull(trial),
          estimate = temp$naive_trt_effects_df %>% pull(VE_est)
        ),
        tibble(
          estimand = "mean_diff_S",
          type = "naive",
          trial = temp$naive_trt_effects_df %>% pull(trial),
          estimate = temp$naive_trt_effects_df %>% pull(mean_diff_S_est)
        )
      )
  }
  )

  
  
}

bootstrap_inferences_df = bootstrap_replicates_df %>%
  group_by(estimand, type, trial) %>%
  summarise(
    CI_lower = quantile(estimate, 0.025, na.rm = TRUE),
    CI_upper = quantile(estimate, 0.975, na.rm = TRUE),
    mean = mean(estimate, na.rm = TRUE),
    se = sd(estimate, na.rm = TRUE)
  )

outfile_bootstrap = paste0(
  "results/raw-results/bootstrap_",
  antibody_type,
  "_",
  target,
  "_",
  population,
  "_M",
  mnum,
  "_",
  outfile_wts,
  ".csv"
)
write.csv(bootstrap_inferences_df, outfile_bootstrap, row.names=FALSE)

## Sandwich Standard Errors ----------------------------------------------

# If some of the estimated weights are equal to 1, this will break the sandwich
# estimator because the corresponding bread matrix will not be invertible. 
problem_weights = RESULT$weights_df %>%
  filter(CC_stratum != "Placebo", weight == 1) %>%
  summarize(n() >= 1) %>%
  pull() %>% # Returns TRUE if there is an estimated weight of exactly 1.
  any()

if (estimate_weights & problem_weights) {
  warning(
    "At least one weight stratum has an estimated weight equal to one. The original sandwich estimator fails because the bread matrix is not invertible. A modified sandwich estimator is obtained by treating the problematic weight stratum's weight as fixed."
  )
  # The problematic stratum is joined with the placebo stratum. This forces the
  # code to treat the corresponding (estimated) weight as fixed. 
  problematic_strata = RESULT$weights_df %>%
    filter(CC_stratum != "Placebo", weight == 1) %>%
    pull(CC_stratum) %>%
    unique()
  df = df %>%
    mutate(CC_stratum = ifelse(CC_stratum %in% problematic_strata, "Placebo", CC_stratum))
  
  # The models and weights are re-estimated, now with the modified weight
  # strata. Note that all results remain unchanged, except that there is now
  # one weight stratum fewer.
  RESULT <- RunDSCAP(
    data = df,
    formula_S = Smod,
    formula_Y = Ymod,
    target_trial = target_trial,
    sim = F, 
    estimate_weights = TRUE
  )
}

# Extract estimated parameter vector that corresponds to the set of stacked
# estimating equations.
theta = extract_coefs(RESULT, estimate_weights)

# Attach estimated weights to original data frame. These weights will be used if
# the weights are treated as known for the sandwich variance estimator.
df = df %>%
  select(-weight) %>% left_join(RESULT$weights_df, by = "CC_stratum")

# Compute sandwich estimate
m_est = m_estimate(
  estFUN = estFUN_taudelta,
  data = df,
  outer_args = list(
    models_tbl = RESULT$glm_fits_df %>%
      select(trial_modified, glm_fit_Y, glm_fit_S),
    weights_df = RESULT$weights_df,
    target_trial = target_trial
  ),
  inner_args = list(
    estimate_weights = estimate_weights
  ),
  roots = theta,
  compute_roots = FALSE,
  deriv_control = setup_deriv_control(method = "simple"),
  compute_vcov = TRUE
)

vcov_m_est <- data.frame(vcov(m_est))
colnames(vcov_m_est) = names(theta)



outfile_vcov = paste0("results/raw-results/vcov_",
                      antibody_type,
                      "_",
                      target,
                      "_",
                      population,
                      "_M",
                      mnum,
                      "_", 
                      outfile_wts,
                      ".csv")
write.csv(vcov_m_est, outfile_vcov, row.names = FALSE)

t2 = Sys.time()
print(t2 - t1)