# Preliminaries -----------------------------------------------------------
t1 <- Sys.time()
 
# Load all R packages
library(tidyverse)
library(furrr) # for parallel computing

# Set seed for reproducibility.
set.seed(2)

# Set up parallel computing
if (parallelly::supportsMulticore()) {
  plan("multicore", workers = parallel::detectCores() - 1)
} else {
  plan(multisession, workers = parallel::detectCores() - 1)
}
# Extract arguments for analysis.
args = commandArgs(trailingOnly = TRUE)

# Number of MC IPD replicates (per trial) to approximate true values.
n_MC <- as.numeric(args[1])

# Source parameter values that also includes a tibble containing all scenarios
# studied.
source("R/simulations/parameter_values.R")

# Change number of subjects per trials (which is not relevant for approximating
# the true values) to a large number `n_MC`.
dgm_param_tbl$n_t = n_MC

# We need one row per setting for generating data.  
dgm_param_tbl = dgm_param_tbl %>%
  group_by(setting) %>%
  slice(1) %>%
  ungroup() %>%
  # The target parameters (trial-level treatment effects etc) do not depend on
  # case-cohort sampling, so we can drop this variable to avoid redundant
  # computations.
  select(-CC_sampling)

# For scenario X (i.e., violation of exchangeability), we need to set
# `different_placebo_model` to TRUE.
dgm_param_tbl <- dgm_param_tbl %>%
  mutate(different_placebo_model = ifelse(setting == "X", TRUE, FALSE))

# Define formulas for the outcome models for (i) clinical outcome Y and (ii)
# surrogate outcome S. The same linear predictors are used for both models. In
# the simulations, the corresponding models are automatically stratified by
# treatment.
formula_Y = as.formula(Y ~ X1 + X2)
formula_S = as.formula(S ~ X1 + X2)

# Define formula for trial participation model given covariates X (but not
# treatment).
formula_T = as.formula(trial ~ X1 + X2)

# Helper Functions --------------------------------------------------------

# Source helper functions.
source("R/helper-functions/DSCAP-estimators.R")
source("R/helper-functions/treatment-effect-estimators.R")
source("R/simulations/generate_simulated_data.R")

## Simulation Function --------------------------------------------------

# Function to simulate and analyze one data set.
simulate_and_analyze <- function(n_trials,
                                 n_t,
                                 gamma,
                                 theta,
                                 zeta,
                                 CC_sampling,
                                 formula_S,
                                 formula_Y,
                                 formula_T,
                                 target_trial,
                                 estimate_weights = FALSE,
                                 alpha = 0.05,
                                 B,
                                 different_placebo_model,
                                 bridge_error) {
  # Generate data.
  simulated_data <- generate_simulated_data(
    n_trials = n_trials,
    n_t = n_t,
    formula_S = formula_S,
    formula_Y = formula_Y,
    gamma = gamma,
    theta = theta,
    zeta = zeta,
    CC_sampling = CC_sampling, 
    target_trial = target_trial,
    different_placebo_model = different_placebo_model,
    bridge_error = bridge_error
  )
  
  # Compute proportion of events in each trial-treatment arm.
  prop_events_tbl <- simulated_data %>%
    group_by(trial, A) %>%
    summarise(prop_events = mean(Y))

  # Analyze data.
  inferences_tbl <- estimate_DSCAP(
    data = simulated_data,
    type = "naive",
    formula_S = formula_S,
    formula_Y = formula_Y,
    formula_T = formula_T,
    trial_var = "trial",
    treatment_var = "A",
    target_trial = target_trial,
    estimate_weights = estimate_weights,
    corrected_target_trial = FALSE,
    CC_weight_var = NULL
  ) %>%
    do.call(what = estimates_table_join)
  
  # Combines estimates with event proportions.
  inferences_tbl <- bind_rows(
    inferences_tbl,
    prop_events_tbl %>%
      rename(treatment = A, estimate = prop_events) %>%
      mutate(measure = "prop_events", type = NA)
  )
  
  return(inferences_tbl)
}

# data analysis -----------------------------------------------------------

## Simulate data and compute estimates  --------------------

# Generate `N_MC` data sets and compute the trial-level treatment effect
# estimates and association measures for each data set.
data_set_indicator = 1

true_values_tbl <- expand_grid(data_set_indicator, dgm_param_tbl) %>%
  # The estimands are compute both the naive and the standardized ones. If
  # `target_trial` = NULL, the data are generated according to the "default" DGP
  # and MC approximation returns naive values. Setting `target_trial` = 1 means
  # that the data are generated according to the covariate distribution of trial
  # 1 and MC approximation returns standardized values.
  expand_grid(tibble(
    target_trial = list(NULL, 1)
  ))
  

true_values_tbl$inferences_tbl = future_pmap(
  .l = list(
    n_trials = true_values_tbl$n_trials,
    n_t = true_values_tbl$n_t,
    gamma = true_values_tbl$gamma,
    theta = true_values_tbl$theta,
    zeta = true_values_tbl$zeta,
    target_trial = true_values_tbl$target_trial,
    different_placebo_model = true_values_tbl$different_placebo_model,
    bridge_error = true_values_tbl$bridge_error
  ),
  .f = simulate_and_analyze,
  formula_S = formula_S,
  formula_Y = formula_Y,
  formula_T = formula_T,
  CC_sampling = FALSE,
  estimate_weights = FALSE,
  .options = furrr_options(
    seed = TRUE,
    stdout = FALSE,
    conditions = character()
  )
)

# Remove redundant information and convert to a long format where each row
# contains one estimate for a given parameter. Hence, the estimates and
# inferences for a single simulated data set span many rows.
true_values_tbl = true_values_tbl %>%
  select(-theta, -gamma, -zeta, -data_set_indicator, -n_t) %>%
  # Recode the `target_trial` variable to be more interpretable.
  rowwise(everything()) %>%
  dplyr::summarize(target_trial = ifelse(is.null(target_trial), "naive", "standardized")) %>%
  ungroup() %>%
  unnest(inferences_tbl)

# Save Results -------------------------------------------------------------

# Save the results in an RDS file.
saveRDS(
  true_values_tbl,
  "results/simulations/raw-results/true_values_tbl.rds"
)

print(Sys.time() - t1)