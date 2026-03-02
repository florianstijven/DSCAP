# Preliminaries -----------------------------------------------------------
t1 <- Sys.time()

# Load all R packages
library(tidyverse)
library(nnet) # multinom() function for multinomial logistic regression
library(Hmisc)
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

# Number of MC replicates in the simulations.
N_MC <- as.numeric(args[1])
# Number of bootstrap replications for inference.
n_boot <- as.numeric(args[2])
# Number of replications for permutation test for conditional independence.r
n_perm <- as.numeric(args[3])

# Source parameter values that also includes a tibble containing all scenarios
# studied.
source("R/simulations/parameter_values.R")

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
formula_CC = as.formula(Delta ~ X1 * trial * as.factor(A) * as.factor(Y))

# If CC sampling is required for the simulations, we set the formula_CC in the
# corresponding row of the tibble containing the DGP parameters. Otherwise, this
# element is null.
dgm_param_tbl$formula_CC <- ifelse(
  dgm_param_tbl$CC_sampling,
  list(formula_CC),
  list(NULL)
)

# For the setting "X", we have a different DGP for the different trial's control
# groups. We need to suuply the `different_placebo_model` argument to
# `different_placebo_model()` later on to indicate this.
dgm_param_tbl$different_placebo_model <- ifelse(
  dgm_param_tbl$setting == "X",
  TRUE,
  FALSE
)

# Helper Functions --------------------------------------------------------

# Source helper functions.
source("R/helper-functions/DSCAP-estimators.R")
source("R/helper-functions/treatment-effect-estimators.R")
source("R/simulations/generate_simulated_data.R")
source("R/helper-functions/permutation-LRT.R")

## Simulation Function --------------------------------------------------

# Function to analyze one data set.
analyze <- function(data,
                    formula_S,
                    formula_Y,
                    formula_T,
                    formula_CC,
                    target_trial,
                    estimate_weights = FALSE,
                    alpha = 0.05,
                    B) {
  # Analyze the simulated data using the naive, standardization, ipw, and doubly
  # robust estimators.
  inferences_naive_tbl <- inference_DSCAP(
    data = data,
    type = "naive",
    formula_S = formula_S,
    formula_Y = formula_Y,
    formula_T = formula_T,
    formula_CC = formula_CC,
    target_trial = target_trial,
    estimate_weights = estimate_weights,
    trial_var = "trial",
    treatment_var = "A",
    alpha = alpha,
    B = B
  )
  inferences_standardization_tbl <- inference_DSCAP(
    data = data,
    type = "standardized",
    formula_S = formula_S,
    formula_Y = formula_Y,
    formula_T = formula_T,
    formula_CC = formula_CC,
    target_trial = target_trial,
    estimate_weights = estimate_weights,
    trial_var = "trial",
    treatment_var = "A",
    alpha = alpha,
    B = B
  )
  inferences_ipw_tbl <- inference_DSCAP(
    data = data,
    type = "ipw",
    formula_S = formula_S,
    formula_Y = formula_Y,
    formula_T = formula_T,
    formula_CC = formula_CC,
    target_trial = target_trial,
    estimate_weights = estimate_weights,
    trial_var = "trial",
    treatment_var = "A",
    alpha = alpha,
    B = B
  )
  inferences_DR_tbl <- inference_DSCAP(
    data = data,
    type = "doubly robust",
    formula_S = formula_S,
    formula_Y = formula_Y,
    formula_T = formula_T,
    formula_CC = formula_CC,
    target_trial = target_trial,
    estimate_weights = estimate_weights,
    trial_var = "trial",
    treatment_var = "A",
    alpha = alpha,
    B = B
  )
  # Join the inference results for the different estimators into one table and
  # add a column indicating the estimator type.
  inferences_tbl <- bind_rows(
    inferences_naive_tbl %>% mutate(type = "naive"),
    inferences_standardization_tbl %>% mutate(type = "standardization"),
    inferences_ipw_tbl %>% mutate(type = "ipw"),
    inferences_DR_tbl %>% mutate(type = "doubly robust")
  )
  
  # Compute permutation p-value for testing the trial exchangeability assumption
  # for control treatment.
  p_value_permutation_Y <- permutation_LRT(
    data = data,
    formula_O = formula_Y,
    formula_CC = formula_CC,
    family = binomial(),
    treatment_var = "A",
    trial_var = "trial",
    n_perm = n_perm,
    estimate_weights = estimate_weights
  )
  
  p_value_permutation_S <- permutation_LRT(
    data = data,
    formula_O = formula_S,
    formula_CC = formula_CC,
    family = gaussian(),
    treatment_var = "A",
    trial_var = "trial",
    n_perm = n_perm,
    estimate_weights = estimate_weights
  )
  
  # Add p-value to inferences table.
  inferences_tbl = inferences_tbl %>%
    bind_rows(
      tibble(
        measure = c("permutation_LRT_p_value_S", "permutation_LRT_p_value_Y"),
        estimate = c(p_value_permutation_S, p_value_permutation_Y),
        type = NA,
        treatment = NA,
        CI_lower_bs = NA,
        CI_upper_bs = NA,
        SE = NA
      )
    )
  
  return(inferences_tbl)
}

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
                                 formula_CC,
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
    CC_sampling = CC_sampling
  )
  
  # Analyze data.
  inferences_tbl <- analyze(
    data = simulated_data,
    formula_S = formula_S,
    formula_Y = formula_Y,
    formula_T = formula_T,
    formula_CC = formula_CC,
    target_trial = target_trial,
    estimate_weights = estimate_weights,
    alpha = alpha,
    B = B
  )
  
  return(inferences_tbl)
}

# data analysis -----------------------------------------------------------

## Simulate data and compute estimates  --------------------

# Generate `N_MC` data sets and compute the trial-level treatment effect
# estimates and association measures for each data set.
data_set_indicator = 1:N_MC

simulations_results_tbl <- expand_grid(data_set_indicator, dgm_param_tbl)

simulations_results_tbl$inferences_tbl = future_pmap(
  .l = list(
    n_trials = simulations_results_tbl$n_trials,
    n_t = simulations_results_tbl$n_t,
    gamma = simulations_results_tbl$gamma,
    theta = simulations_results_tbl$theta,
    zeta = simulations_results_tbl$zeta,
    CC_sampling = simulations_results_tbl$CC_sampling,
    formula_CC = simulations_results_tbl$formula_CC
  ),
  .f = simulate_and_analyze,
  formula_S = formula_S,
  formula_Y = formula_Y,
  formula_T = formula_T,
  target_trial = 1,
  estimate_weights = FALSE,
  alpha = 0.05,
  B = n_boot,
  .options = furrr_options(
    seed = TRUE,
    stdout = FALSE,
    conditions = character()
  )
)

# Remove redundant information and convert to a long format where each row
# contains one estimate for a given parameter. Hence, the estimates and
# inferences for a single simulated data set span many rows.
simulations_results_tbl = simulations_results_tbl %>%
  select(-theta, -gamma, -zeta) %>%
  unnest(inferences_tbl)


# Save Results -------------------------------------------------------------

# Save the results in an RDS file.
saveRDS(
  simulations_results_tbl,
  "results/simulations/raw-results/simulations_results_tbl.rds"
)

print(Sys.time() - t1)