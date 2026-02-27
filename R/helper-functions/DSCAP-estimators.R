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
                           CC_weight_var,
                           corrected_target_trial = TRUE) {
  # Estimate the trial-level treatment effects with the requested estimator.
  if (type == "standardized") {
    means_df = data %>%
      standardization_estimator(
        formula_Y = formula_Y,
        formula_S = formula_S,
        formula_CC = formula_CC,
        trial_var = trial_var,
        treatment_var = treatment_var,
        target_trial = target_trial,
        estimate_weights = estimate_weights,
        CC_weight_var = CC_weight_var
      ) 
    
    trt_effects_df = means_df %>%
      trt_effects(by_trial = FALSE)
  } else if (type == "ipw") {
    means_df = data %>%
      ipw_estimator(
        formula_T = formula_T,
        formula_CC = formula_CC,
        trial_var = trial_var,
        treatment_var = treatment_var,
        target_trial = target_trial,
        estimate_weights = estimate_weights,
        CC_weight_var = CC_weight_var
      ) 
    
    trt_effects_df = means_df %>%
      trt_effects(by_trial = FALSE)
  } else if (type == "doubly robust") {
    means_df = data %>%
      DR_estimator(
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
    
    trt_effects_df = means_df %>%
      trt_effects(by_trial = FALSE)
  } else if (type == "naive") {
    means_df = data %>%
      naive_estimator(
        formula_CC = formula_CC,
        trial_var = trial_var,
        treatment_var = treatment_var,
        target_trial = target_trial,
        estimate_weights = estimate_weights,
        CC_weight_var = CC_weight_var
      ) 
    
    trt_effects_df = means_df %>%
      trt_effects(by_trial = TRUE)
  }
  
  # If any estimated proportions are zero, add 0.5 divided by that treatment
  # group's sample size to the estimated proportion.
  if (any(means_df$mean_Y == 0)) {
    if (type == "naive") {
      means_df = means_df %>%
        group_by(trial, treatment) %>%
        mutate(mean_Y = ifelse(mean_Y == 0, 0.5 / n(), mean_Y)) %>%
        ungroup()
    } else {
      means_df = means_df %>%
        group_by(treatment) %>%
        mutate(mean_Y = ifelse(mean_Y == 0, 0.5 / n(), mean_Y)) %>%
        ungroup()
    }
  }
  
  # If `corrected_target_trial` is TRUE, then we will use the treatment effect
  # estimates based on corrected estimators (ipw, standardized, or
  # doubly-robust) for the target trial. Otherwise, the treatment effect
  # estimate on the target trial is based on the uncorrected mean in the target
  # trial's active arm. 
  if (!corrected_target_trial & type != "naive") {
    # Estimate uncorrected mean in the target trial's active arm and use the
    # control group's estimated mean (based on whatever method) to compute the
    # treatment effects in the target trial.
    naive_trt_effect_target_df = data %>%
      naive_estimator(
        formula_CC = formula_CC,
        trial_var = trial_var,
        treatment_var = treatment_var,
        target_trial = target_trial,
        estimate_weights = estimate_weights,
        CC_weight_var = CC_weight_var
      ) %>%
      filter(trial == target_trial, treatment != 0) %>%
      bind_rows(means_df %>%
                  filter(treatment == 0)) %>%
      trt_effects(by_trial = FALSE)
    # Replace the treatment effect estimates in the target trial with the
    # uncorrected estimates.
    trt_effects_df = trt_effects_df %>%
      filter(treatment != naive_trt_effect_target_df$treatment) %>%
      bind_rows(naive_trt_effect_target_df)
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
                            CC_weight_var = NULL,
                            B,
                            stratified_bs = TRUE,
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
  # the association measures. If `stratified_bs` is TRUE, then we will resample
  # within each trial.
  if (stratified_bs) {
    strata = data[[trial_var]]
  } else {
    strata = rep(1, nrow(data))
  }
  
  boot_out <- boot::boot(
    data = data,
    statistic = boot_fun,
    R = B,
    strata = strata
  )
  # Compute bootstrap CIs and SEs for all estimates in estimates_df.
  estimates_df = estimates_df %>%
    mutate(
      CI_lower_bs = purrr::map_dbl(1:nrow(estimates_df), function(.x) {
        ci_limit <- tryCatch(
          expr = {
            boot::boot.ci(boot_out,
                          type = "perc",
                          index = .x,
                          conf = 1 - alpha)$percent[4]
          },
          error = function(e) {
            NA
          }
        ) 
        if (!is.numeric(ci_limit)) {
          ci_limit = NA
        }
        return(ci_limit)
      }),
      CI_upper_bs = purrr::map_dbl(1:nrow(estimates_df), function(.x) {
        ci_limit <- tryCatch(
          expr = {
            boot::boot.ci(boot_out,
                          type = "perc",
                          index = .x,
                          conf = 1 - alpha)$percent[5]
          },
          error = function(e) {
            NA
          }
        ) 
        if (!is.numeric(ci_limit)) {
          ci_limit = NA
        }
        return(ci_limit)
      }),
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
  if (any(is.na(trt_effects_df$VE_est)) | any(is.na(trt_effects_df$mean_diff_S_est))) {
    browser()
  }
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