# Function that computes the trial-level treatment effects and DSCAPs (i.e., the
# association measures) for a given data set from multiple trials. The `type`
# argument specifies the type of trial-level treatment effect estimators
# (standardized, ipw, doubly-robust, naive).
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
  
  return(
    list(
      association_measures_df = association_measures_df,
      trt_effects_df = trt_effects_df
    )
  )
}

# Function that computes the trial-level treatment effect estimates and
# association measures, and performs the bootstrap for inferences.
inference_DSCAP <- function(data,
                            type,
                            formula_Y,
                            formula_S,
                            formula_T,
                            formula_CC = NULL,
                            trial_var = "trial",
                            treatment_var,
                            target_trial,
                            estimate_weights = FALSE,
                            CC_weight_var,
                            B,
                            alpha = 0.05) {
  # Estimate the trial-level treatment effects and association measures for the
  # original data.
  estimates_df <- estimate_DSCAP(
    data = data,
    type = type,
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
    do.call(what = estimates_table_join)

  # Function that takes in data and indices, and return the estimates.
  boot_fun <- function(data, indices) {
    data_resampled <- data[indices, ] # allows boot to select sample
    estimates_list = estimate_DSCAP(
      data = data_resampled,
      type = type,
      formula_Y = formula_Y,
      formula_S = formula_S,
      formula_T = formula_T,
      formula_CC = formula_CC,
      trial_var = trial_var,
      treatment_var = treatment_var,
      target_trial = target_trial,
      estimate_weights = estimate_weights,
      CC_weight_var = CC_weight_var
    )
    estimates_table_join(estimates_list[[1]], estimates_list[[2]]) %>%
      pull(estimate) %>%
      return()
  }
  
  # Perform the bootstrap to obtain standard errors and confidence intervals for
  # the association measures.
  boot_out <- boot::boot(data = df, statistic = boot_fun, R = B)
  # Compute bootstrap CIs and SEs for all estimates in estimates_df.
  estimates_df = estimates_df %>%
    mutate(
      CI_lower_bs = purrr::map_dbl(
        1:nrow(estimates_df),
        ~ boot::boot.ci(boot_out, type = "perc", index = .x, conf = 1 - alpha)$percent[4]
      ),
      CI_upper_bs = purrr::map_dbl(
        1:nrow(estimates_df),
        ~ boot::boot.ci(boot_out, type = "perc", index = .x, conf = 1 - alpha)$percent[5]
      ),
      SE_bs = apply(boot_out$t, 2, sd)
    )
  
  return(estimates_df)
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

# Function that maps the association measures and treatment effect measures
# tables into a single table with one estimate per row.
estimates_table_join <- function(association_measures_df, trt_effects_df) {
  # Reshape the association measures table to have one row per measure.
  association_measures_long_df = association_measures_df %>%
    pivot_longer(cols = everything(),
                 names_to = "measure",
                 values_to = "estimate") %>%
    mutate(type = "association_measure")
  
  # Reshape the treatment effects table to have one row per trial and measure.
  trt_effects_long_df = trt_effects_df %>%
    pivot_longer(cols = c(VE_est, mean_diff_S_est),
                 names_to = "measure",
                 values_to = "estimate") %>%
    mutate(type = "treatment_effect")
  
  # Join the two tables together.
  estimates_df = bind_rows(association_measures_long_df, trt_effects_long_df)
  
  return(estimates_df)
}