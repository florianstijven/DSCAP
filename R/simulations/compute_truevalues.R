# Preliminaries -----------------------------------------------------------
t1 <- Sys.time()

# Load all R packages
library(tidyverse)
library(furrr) # for parallel computing

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

# Define formulas for the outcome models for (i) clinical outcome Y and (ii)
# surrogate outcome S. The same linear predictors are used for both models. In
# the simulations, the corresponding models are automatically stratified by
# treatment.
formula_Y = as.formula(Y ~ X1 + X2)
formula_S = as.formula(S ~ X1 + X2)

# Define formula for trial participation model given covariates X (but not
# treatment).
formula_T = as.formula(trial ~ X1 + X2)
# Formula for case-cohort sampling model (i.e., model for probability of being
# sampled in the case-cohort sampling design given covariates X, treatment A,
# outcome Y, and trial).
formula_CC = as.formula(Delta ~ CC_stratum * trial * as.factor(A))

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
                                 ZOPT = NA,
                                 formula_S,
                                 formula_Y,
                                 formula_T,
                                 target_trial,
                                 estimate_weights = FALSE,
                                 alpha = 0.05,
                                 B) {
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
    target_trial = 1
  )
  
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
    corrected_target_trial = FALSE
  ) %>%
    do.call(what = estimates_table_join)
  
  return(inferences_tbl)
}

# data analysis -----------------------------------------------------------

## Simulate data and compute estimates  --------------------

# Generate `N_MC` data sets and compute the trial-level treatment effect
# estimates and association measures for each data set.
data_set_indicator = 1

true_values_tbl <- expand_grid(data_set_indicator, dgm_param_tbl)

true_values_tbl$inferences_tbl = future_pmap(
  .l = list(
    n_trials = simulations_results_tbl$n_trials,
    n_t = simulations_results_tbl$n_t,
    gamma = simulations_results_tbl$gamma,
    theta = simulations_results_tbl$theta,
    zeta = simulations_results_tbl$zeta,
    CC_sampling = simulations_results_tbl$CC_sampling
  ),
  .f = simulate_and_analyze,
  formula_S = formula_S,
  formula_Y = formula_Y,
  formula_T = formula_T,
  target_trial = 1,
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
  select(-theta, -gamma, -data_set_indicator) %>%
  unnest(inferences_tbl)

# Save Results -------------------------------------------------------------

# Save the results in an RDS file.
saveRDS(
  true_values_tbl,
  "results/simulations/raw-results/true_values_tbl.rds"
)

print(Sys.time() - t1)