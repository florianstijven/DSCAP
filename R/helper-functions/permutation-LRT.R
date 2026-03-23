# Function to compute the LRT test statistic.
LRT_statistic <- function(data_subset, formula_O, family) {
  # Fit null model, pooling data from the control groups of all trials and not
  # including the trial variable as predictor.
  if (family$family == "binomial") {
    # If the outcome is binomial, there is not case-cohort sampling and all
    # observations should be used.
    data_subset$CC_weight = 1
    data_subset$Delta = 1
  }
  glm_null <- glm(formula = formula_O,
                  family = family,
                  data = data_subset,
                  subset = (Delta == 1),
                  weights = CC_weight,
                  x = FALSE,
                  y = FALSE)
  # Fit alternative model that contains the interaction between all predictors
  # in the null model and the trial variable.
  glm_w_trial <-update(glm_null, ~ . * trial)
  # Compute the likelihood ratio test statistic. 
  glm_lr_test = lmtest::lrtest(glm_null, glm_w_trial)

  # Extract LRT statistic.
  lrt_test <- c("statistic" = glm_lr_test[2,4], "p-value" = glm_lr_test[2,5])
  
  return(lrt_test)
}


# Function to compute the LRT test permutation p-value.
permutation_LRT <- function(data,
                            formula_O,
                            family,
                            formula_CC,
                            treatment_var,
                            trial_var,
                            n_permutations,
                            estimate_weights,
                            CC_weight_var = NULL) {
  df <- data_preparation(
    data = data,
    formula_CC = formula_CC,
    trial_var = trial_var,
    treatment_var = treatment_var,
    target_trial = 1,
    estimate_weights = estimate_weights,
    CC_weight_var = CC_weight_var
  )
  
  df_subset = df %>%
    filter(treatment == 0)
  # Compute observed LRT test statistic for the original data.
  lrt_test_statistic <- LRT_statistic(data_subset = df_subset, formula_O = formula_O, family = family)["statistic"]
  
  # Initialize vector to store LRT test statistics from the permutations.
  lrt_test_statistic_permuted = rep(NA, n_permutations)
  # Vector of trial labels to permute.
  trial_vec = df_subset$trial
  # Permute the data `n_permutations` times, and compute the LRT test statistic
  # for each permutation.
  for (i in 1:n_permutations) {
    # Permute trial labels.
    df_subset$trial <- sample(trial_vec, replace = FALSE)
    # Compute LRT test statistic for the permuted data.
    lrt_test_statistic_permuted[i] <- LRT_statistic(data_subset = df_subset,
                                                    formula_O = formula_O,
                                                    family = family)["statistic"]
  }
  # Compute the permutation p-value.
  p_value <- mean(lrt_test_statistic_permuted > lrt_test_statistic)
  
  return(p_value)
}

# Function to compute the LRT p-value based on the chi-squared asymptotic
# distribution.
asymptotic_LRT <- function(data,
                           formula_O,
                           family,
                           formula_CC,
                           treatment_var,
                           trial_var,
                           estimate_weights,
                           CC_weight_var = NULL) {
  df <- data_preparation(
    data = data,
    formula_CC = formula_CC,
    trial_var = trial_var,
    treatment_var = treatment_var,
    target_trial = 1,
    estimate_weights = estimate_weights,
    CC_weight_var = CC_weight_var
  )
  
  df_subset = df %>%
    filter(treatment == 0)
  
  # Compute LRT test statistic and p-value based on the chi-squared asymptotic distribution.
  lrt_test <- LRT_statistic(data_subset = df_subset, formula_O = formula_O, family = family)
  
  return(lrt_test["p-value"])
}