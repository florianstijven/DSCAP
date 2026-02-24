RunDSCAP <- function(data,
                     formula_Y,
                     formula_S,
                     formula_T, 
                     target_trial,
                     #  sim = F,
                     estimate_weights = FALSE,
                     CC_sampling = FALSE,
                     alpha_level = 0.05)
{
  df = data
  df$trial <- factor(df$trial, levels = 1:5)
  
  # Determine which level of the trial factor corresponds to the target trial. 
  target_trial_level_n = which(levels(df$trial) == target_trial)
  
  n_trials <- length(levels(df$trial))
  
  # Re-estimate weights if required.
  if (estimate_weights) {
    df <- df %>% group_by(CC_stratum) %>% mutate(prop.delta = sum(Delta == 1) /
                                                   n(),
                                                 weight = 1 / prop.delta) %>%
      dplyr::select(-prop.delta) %>%
      # Placebo patients should get a weight of one.
      mutate(weight = ifelse(A == 0, 1, weight))
  }
  
  # Data frame that contains the estimated weight for each category in
  # `CC_stratum`.
  if(CC_sampling == TRUE){
    weights_df = df %>%
      group_by(CC_stratum) %>%
      slice_head(n = 1) %>%
      dplyr::select(CC_stratum, weight)
    
    problematic_strata = NA
    if (any(weights_df$weight == Inf)) {
      problematic_strata = weights_df %>%
        filter(weight == Inf) %>%
        pull(CC_stratum)
      warning(paste0(
        paste(
          c("Estimated weight(s) of Inf for strata", problematic_strata),
          collapse = " "
        ),
        ". Corresponding observations are ignored for this analysis."
      ))
      weights_df = weights_df %>%
        filter(!(CC_stratum %in% problematic_strata))
      df = df  %>%
        filter(!(CC_stratum %in% problematic_strata))
    }
  }
  
  
  # standardization ---------------------------------------------------------
  
  if(CC_sampling == TRUE){
    outcome_model_fits_df = df %>%
      mutate(trial_modified = ifelse(A == 0, "Placebo", as.character(trial)),
             trial_modified = factor(trial_modified,
                                     levels = c("Placebo", "1", "2", "3", "4", "5"))) %>%
      filter(trial != target_trial) %>%
      group_by(trial_modified) %>%
      dplyr::summarize(
        outcome_model_fit_Y = glm(formula_Y, data = pick(everything()), family = binomial, x = FALSE, y = FALSE) %>% list(),
        outcome_model_fit_S = glm(
          formula_S,
          data = pick(everything()),
          family = gaussian,
          subset = (Delta == 1),
          weights = weight, 
          x = FALSE,
          y = FALSE
        ) %>% list()
      ) %>%
      # Extract coefficients vectors.
      ungroup() %>%
      mutate(
        outcome_model_coef_Y = purrr::map(outcome_model_fit_Y, coef),
        outcome_model_coef_S = purrr::map(outcome_model_fit_S, coef)
      )
  } else{
    outcome_model_fits_df = df %>%
      mutate(trial_modified = ifelse(A == 0, "Placebo", as.character(trial)),
             trial_modified = factor(trial_modified,
                                     levels = c("Placebo", "1", "2", "3", "4", "5"))) %>%
      filter(trial != target_trial) %>%
      group_by(trial_modified) %>%
      dplyr::summarize(
        outcome_model_fit_Y = glm(formula_Y, data = pick(everything()), family = binomial, x = FALSE, y = FALSE) %>% list(),
        outcome_model_fit_S = glm(
          formula_S,
          data = pick(everything()),
          family = gaussian,
          # subset = (Delta == 1),
          # weights = weight, 
          x = FALSE,
          y = FALSE
        ) %>% list()
      ) %>%
      # Extract coefficients vectors.
      ungroup() %>%
      mutate(
        outcome_model_coef_Y = purrr::map(outcome_model_fit_Y, coef),
        outcome_model_coef_S = purrr::map(outcome_model_fit_S, coef)
      )
  }
  
  # For each subject in the target trial, we predict the outcome using the model
  # estimated in any of the other trials.
  target_trial_df = df %>% 
    filter(trial == target_trial) 
  
  outcome_model_fits_df = outcome_model_fits_df %>%
    mutate(
      predicted_Y = purrr::map(
        .x = outcome_model_fit_Y,
        .f = function(outcome_model_fit) {
          predict(outcome_model_fit, newdata = target_trial_df, type = "response")
        }
      ),
      predicted_S = purrr::map(
        .x = outcome_model_fit_S,
        .f = function(outcome_model_fit) {
          predict(outcome_model_fit, newdata = target_trial_df, type = "response")
        }
      )
    )
  
  # Compute the standardized means as the means of the predicted outcomes, based
  # on any of the trial-specific models, of subjects from the target trial. 
  standardized_means_df = outcome_model_fits_df %>%
    mutate(
      mean_Y = purrr::map_dbl(
        .x = predicted_Y,
        .f = mean
      ),
      mean_S = purrr::map_dbl(
        .x = predicted_S,
        .f = mean
      )
    ) %>%
    dplyr::select(trial_modified, mean_Y, mean_S)
  
  
  # weighted (IPW) estimator ----------------------------------------------------------------
  
  
  # conditional treatment assignment probability is treated 
  # as a known function of T and X.
  # note that this table does not include X2, just X1 
  conditional_treatment_assignment_table <- expand.grid(
    X1 = unique(df$X1), A = unique(df$A), trial = unique(df$trial)) |> 
    mutate(treatment_assignment_prob = case_when(
      A == 0 ~ 0.5,
      A == 1 & trial == 1 ~ 0.5,
      A == 2 & trial == 2 ~ 0.5,
      A == 3 & trial == 3 ~ 0.5,
      A == 4 & trial == 4 ~ 0.5,
      A == 5 & trial == 5 ~ 0.5,
      .default = 0
    )) |> 
    arrange(trial, A, X1)
  
  df_ipw <- df %>%
    mutate(trial_modified = ifelse(A == 0, "Placebo", as.character(trial)),
           trial_modified = factor(trial_modified,
                                   levels = c("Placebo", "1", "2", "3", "4", "5")),
           trial = relevel(trial, ref = "1")
    )
  
  # model for T | X 
  prob_t_equals_1 <- nrow(df_ipw[df_ipw$trial==1,]) / nrow(df_ipw)
  trial_participation_model_fit <- nnet::multinom(formula_T, df_ipw)
  
  df_ipw$trial_participation_model_preds <- asplit(fitted(trial_participation_model_fit), 1)
  
  calculate_W_hat_a <- function(conditional_treatment_assignment_table,
                                trial_participation_model_fit, 
                                trial_participation_model_preds,
                                prob_t_equals_1,
                                A, 
                                X1, 
                                X2) {
    
    conditional_treatment_assignment_probs_by_trial <-  
      conditional_treatment_assignment_table[(conditional_treatment_assignment_table$A == A) & 
                                               (conditional_treatment_assignment_table$X1 == X1),]$treatment_assignment_prob
    
    W_hat_a <- trial_participation_model_preds[1] /
      (prob_t_equals_1 * sum(
        conditional_treatment_assignment_probs_by_trial * 
          trial_participation_model_preds
      ))
    
    return(W_hat_a)
  }
  
  df_ipw <- df_ipw %>%
    mutate(W_hat_a = calculate_W_hat_a(conditional_treatment_assignment_table,
                                       trial_participation_model_fit, 
                                       trial_participation_model_preds,
                                       prob_t_equals_1,
                                       A, 
                                       X1, 
                                       X2)) 
  
  ipw_means_df <- df_ipw |> 
    #filter(trial != target_trial) |> 
    group_by(trial_modified) |> 
    dplyr::summarize(mean_Y = sum(Y * W_hat_a) / nrow(df_ipw),
                     mean_S = sum(S * W_hat_a) / nrow(df_ipw)) |> 
    filter(trial_modified != 1) 
  
  # naive estimates and output all estimates --------------------------------
  
  # Compute the naive trial-specific means. 
  if(CC_sampling == TRUE){
    naive_means_df = df %>%
      group_by(trial, vax) %>%
      dplyr::summarize(
        mean_Y = mean(Y),
        n = n(),
        mean_S = mean(weight * S),
        var_S = var(weight * S),
        .groups = "drop"
      )
  } else{
    naive_means_df = df %>%
      group_by(trial, vax) %>%
      dplyr::summarize(
        mean_Y = mean(Y),
        n = n(),
        mean_S = mean(S),
        var_S = var(S),
        .groups = "drop"
      )
  }
  
  # Compute treatment effects from the estimated means. For the naive
  # trial-specific estimates, we also compute confidence intervals. 
  naive_trt_effects_df = naive_means_df %>%
    pivot_wider(
      names_from = c("vax"),
      values_from = c("mean_Y", "n", "mean_S", "var_S")
    ) %>%
    mutate(
      VE_est = 1 - mean_Y_1 / mean_Y_0,
      log_RR_est = log(1 - VE_est),
      log_RR_se = sqrt(((1 - mean_Y_1) / (mean_Y_1 * n_1)) + ((1 - mean_Y_0) / (mean_Y_0 * n_0))),
      mean_diff_S_est = mean_S_1 - mean_S_0,
      mean_diff_S_se = sqrt((var_S_1 / n_1)  + (var_S_0 / n_0)),
    ) %>%
    mutate(
      VE_CI_lower = 1 - exp(log_RR_est + qnorm(1 - alpha_level / 2) * log_RR_se),
      VE_CI_upper = 1 - exp(log_RR_est - qnorm(1 - alpha_level / 2) * log_RR_se),
      mean_diff_S_CI_lower = mean_diff_S_est - qnorm(1 - alpha_level / 2) * mean_diff_S_se,
      mean_diff_S_CI_upper = mean_diff_S_est + qnorm(1 - alpha_level / 2) * mean_diff_S_se
    )
  
  standardized_trt_effects_df = standardized_means_df %>%
    bind_rows(
      naive_means_df %>%
        filter(trial == target_trial, vax == 1) %>%
        rename(trial_modified = trial) %>%
        mutate(trial_modified = factor(trial_modified,
                                       levels = c("Placebo", "1", "2", "3", "4", "5"))) |> 
        dplyr::select(trial_modified, mean_Y, mean_S)
    ) %>%
    filter(trial_modified != "Placebo") %>%
    rename(mean_Y_1 = mean_Y, mean_S_1 = mean_S) %>%
    cross_join(
      standardized_means_df %>%
        filter(trial_modified == "Placebo") %>%
        rename(mean_Y_0 = mean_Y, mean_S_0 = mean_S) %>%
        dplyr::select(-trial_modified)
    ) %>%
    mutate(VE_est = 1 - mean_Y_1 / mean_Y_0,
           mean_diff_S_est = mean_S_1 - mean_S_0,) %>%
    rename(trial = "trial_modified")
  
  ipw_trt_effects_df = ipw_means_df %>%
    bind_rows(
      naive_means_df %>%
        filter(trial == target_trial, vax == 1) %>%
        rename(trial_modified = trial) %>%
        mutate(trial_modified = factor(trial_modified,
                                       levels = c("Placebo", "1", "2", "3", "4", "5"))) |> 
        dplyr::select(trial_modified, mean_Y, mean_S)
    ) %>%
    filter(trial_modified != "Placebo") %>%
    rename(mean_Y_1 = mean_Y, mean_S_1 = mean_S) %>%
    cross_join(
      ipw_means_df %>%
        filter(trial_modified == "Placebo") %>%
        rename(mean_Y_0 = mean_Y, mean_S_0 = mean_S) %>%
        dplyr::select(-trial_modified)
    ) %>%
    mutate(VE_est = 1 - mean_Y_1 / mean_Y_0,
           mean_diff_S_est = mean_S_1 - mean_S_0,) %>%
    rename(trial = "trial_modified")
  
  # Compute the measures of surrogacy based on the naive trial-specific and
  # standardized treatment effect estimates.
  cor_naive_df = naive_trt_effects_df %>%
    group_by() %>%
    dplyr::summarize(
      cor_p = cor(VE_est, mean_diff_S_est, method = "pearson"),
      cor_s = cor(VE_est, mean_diff_S_est, method = "spearman"),
      beta = coef(lm(VE_est ~ mean_diff_S_est))[2]
    )
  
  cor_standardized_df = standardized_trt_effects_df %>%
    group_by() %>%
    dplyr::summarize(
      cor_p = cor(VE_est, mean_diff_S_est, method = "pearson"),
      cor_s = cor(VE_est, mean_diff_S_est, method = "spearman"),
      beta = coef(lm(VE_est ~ mean_diff_S_est))[2]
    )
  
  cor_ipw_df = ipw_trt_effects_df %>%
    group_by() %>%
    dplyr::summarize(
      cor_p = cor(VE_est, mean_diff_S_est, method = "pearson"),
      cor_s = cor(VE_est, mean_diff_S_est, method = "spearman"),
      beta = coef(lm(VE_est ~ mean_diff_S_est))[2]
    )
  
  
  # Combine everything into a list and return this list. 
  if(CC_sampling == TRUE){
    obj = 
      list(
        outcome_model_fits_df = outcome_model_fits_df %>%
          dplyr::select(-predicted_Y, -predicted_S),
        standardized_means_df = standardized_means_df,
        standardized_trt_effects_df = standardized_trt_effects_df,
        naive_means_df = naive_means_df,
        naive_trt_effects_df = naive_trt_effects_df,
        cor_standardized_df = cor_standardized_df,
        cor_naive_df = cor_naive_df,
        weights_df = weights_df,
        target_trial = target_trial,
        excluded_CC_strata = problematic_strata
      )
  } else{
    obj = 
      list(
        outcome_model_fits_df = outcome_model_fits_df %>%
          dplyr::select(-predicted_Y, -predicted_S),
        standardized_means_df = standardized_means_df,
        standardized_trt_effects_df = standardized_trt_effects_df,
        ipw_means_df = ipw_means_df,
        ipw_trt_effects_df = ipw_trt_effects_df,
        naive_means_df = naive_means_df,
        naive_trt_effects_df = naive_trt_effects_df,
        cor_standardized_df = cor_standardized_df,
        cor_naive_df = cor_naive_df,
        cor_ipw_df = cor_ipw_df,
        #weights_df = weights_df,
        target_trial = target_trial
        # excluded_CC_strata = problematic_strata
      )
  }
  
  gc()
  
  return(obj)
}

extract_coefs = function(obj, estimate_weights, target_trial) {
  # The ordering of the trial-specific estimates is determined by the ordering
  # of the rows in `obj$outcome_model_fits_df`.
  trials_modified_ordering = obj$outcome_model_fits_df %>%
    pull(trial_modified) %>%
    as.character()
  # Ordering excluding placebo.
  trials_ordering = obj$outcome_model_fits_df %>%
    filter(trial_modified != "Placebo") %>%
    pull(trial_modified) %>%
    as.character()
  
  # Extract coefficients of the outcome regression models for Y.
  ncov = length(obj$outcome_model_fits_df$outcome_model_coef_Y[[1]])
  coefs_vec_Y = obj$outcome_model_fits_df %>%
    pull(outcome_model_coef_Y) %>%
    unlist()
  names(coefs_vec_Y) = paste(names(coefs_vec_Y), rep(trials_modified_ordering, each = ncov), sep = " - ")
  # Extract coefficients of the outcome regression models for S. We don't need
  # the model for the placebo group.
  coefs_vec_S = obj$outcome_model_fits_df %>%
    filter(trial_modified != "Placebo") %>%
    pull(outcome_model_coef_S, name = "trial_modified") %>%
    unlist()
  names(coefs_vec_S) = paste(names(coefs_vec_S), rep(trials_ordering, each = ncov), sep = " - ")
  
  # Extract standardized mean estimates for Y and S. 
  standardized_mean_Y = left_join(
    tibble(trial_modified = trials_modified_ordering),
    obj$standardized_means_df) %>%
    pull(mean_Y)
  names(standardized_mean_Y) = paste("standardized_mean_Y", trials_modified_ordering, sep = " - ")
  
  standardized_mean_S = left_join(
    tibble(trial_modified = trials_ordering),
    obj$standardized_means_df) %>%
    pull(mean_S)
  names(standardized_mean_S) = paste("standardized_mean_S", trials_ordering, sep = " - ")
  
  # Extract mean estimates for Y and S for the target trial.
  target_trial_mean_Y = obj$naive_means_df %>%
    filter(trial == obj$target_trial, vax == 1) %>%
    pull(mean_Y)
  names(target_trial_mean_Y) = paste("mean_Y", target_trial, sep = " - ")
  
  target_trial_mean_S = obj$naive_means_df %>%
    filter(trial == obj$target_trial, vax == 1) %>%
    pull(mean_S)
  names(target_trial_mean_S) = paste("mean_S", target_trial, sep = " - ")
  
  # Extract naive treatment effect estimate for the target trial. 
  target_trial_VE_est = obj$standardized_trt_effects_df %>% filter(trial == obj$target_trial) %>%
    pull(VE_est, name = trial)
  names(target_trial_VE_est) = paste("VE_est", target_trial, sep = " - ")
  
  # Extract standardized VE estimates. 
  standardized_VE_est = left_join(tibble(trial = trials_ordering),
                                  obj$standardized_trt_effects_df) %>%
    pull(VE_est)
  names(standardized_VE_est) = paste("standardized_VE_est", trials_ordering, sep = " - ")
  
  # Extract estimated Pearson correlation and regression slope for the
  # standardized estimates.
  standardized_cor_p = obj$cor_standardized_df %>%
    pull(cor_p)
  names(standardized_cor_p) = "standardized_cor_p"
  standardized_beta = obj$cor_standardized_df %>%
    pull(beta)
  names(standardized_beta) = "standardized_beta"
  
  # Join all estimates.
  estimates_vec = c(
    coefs_vec_Y,
    coefs_vec_S,
    target_trial_mean_Y,
    standardized_mean_Y,
    target_trial_mean_S,
    standardized_mean_S,
    target_trial_VE_est,
    standardized_VE_est,
    standardized_cor_p,
    standardized_beta
  )
  # If the weights are treated as being estimated, they are included in the
  # returned vector of estimates.
  if (estimate_weights) {
    estimated_weights = obj$weights_df %>%
      # The weight for placebo patients is known anyhow, the
      # corresponding weight is thus always treated as known.
      filter(CC_stratum != "Placebo") %>% 
      pull(weight, name = CC_stratum)
    
    estimates_vec = c(estimates_vec, 1 / estimated_weights)
  }
  return(estimates_vec)
}


