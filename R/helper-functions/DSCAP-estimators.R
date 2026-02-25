# Function that computes the DSCAPs (i.e., the association measures) for a given
# data set from multiple trials. The `type` argument specifies the type of
# trial-level treatment effect estimators (standardized, ipw, doubly-robust,
# naive).
estimate_DSCAP <- function(data,
                           type,
                           formula_Y,
                           formula_S,
                           formula_T,
                           formula_CC = NULL,
                           trial_var = "trial",
                           treatment_var,
                           target_trial,
                           estimate_weights = FALSE,
                           CC_weight_var) {
  # Estimate the trial-level treatment effects with the requested estimator.
  if (type == "standardized") {
    trt_effects_df = data %>%
      standardization_estimator(
        formula_Y = formula_Y,
        formula_S = formula_S,
        formula_CC = formula_CC,
        trial_var = trial_var,
        treatment_var = treatment_var,
        target_trial = target_trial,
        estimate_weights = estimate_weights,
        CC_weight_var = CC_weight_var
      ) %>%
      trt_effects(by_trial = FALSE)
  } else if (type == "ipw") {
    trt_effects_df = data %>%
      ipw_estimator(
        formula_T = formula_T,
        formula_CC = formula_CC,
        trial_var = trial_var,
        treatment_var = treatment_var,
        target_trial = target_trial,
        estimate_weights = estimate_weights,
        CC_weight_var = CC_weight_var
      ) %>%
      trt_effects(by_trial = FALSE)
  } else if (type == "doubly-robust") {
    trt_effects_df = data %>%
      doubly_robust_estimator(
        formula_Y = formula_Y,
        formula_S = formula_S,
        formula_T = formula_T,
        formula_CC = formula_CC,
        trial_var = trial_var,
        treatment_var = treatment_var,
        target_trial = target_trial,
        estimate_weights = estimate_weights,
        CC_weight_var = CC_weight_var
      ) %>%
      trt_effects(by_trial = FALSE)
  } else if (type == "naive") {
    trt_effects_df = data %>%
      naive_estimator(
        formula_CC = formula_CC,
        trial_var = trial_var,
        treatment_var = treatment_var,
        estimate_weights = estimate_weights,
        CC_weight_var = CC_weight_var
      ) %>%
      trt_effects(by_trial = TRUE)
  }
  
  # Estimate the association measures given the trial-level treatment effect
  # estimates.
  association_measures_df = association_measures(trt_effects_df)
  
  return(association_measures_df)
}

# Function to compute association measures (i.e., the DSCAPs) given the table of
# treatment effect estimates.
association_measures <- function(trt_effects_df) {
  # Pearson correlation
  cor_p_est = cor(trt_effects_df$VE_est, trt_effects_df$mean_diff_S_est)
  
  # Spearman correlation
  cor_s_est = cor(trt_effects_df$VE_est,
                  trt_effects_df$mean_diff_S_est,
                  method = "spearman")
  
  # Linear regression slope
  beta_est = coef(lm(VE_est ~ mean_diff_S_est, data = trt_effects_df))[2]
  
  association_measures_df = tibble(
    cor_p_est = cor_p_est,
    cor_s_est = cor_s_est,
    beta_est = beta_est
  )
  
  return(association_measures_df)
}