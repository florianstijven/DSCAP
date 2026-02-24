RunDSCAP <- function(data,
                     formula_Y,
                     formula_S,
                     formula_T,
                     formula_CC,
                     trial_var = "trial",
                     treatment_var,
                     target_trial,
                     estimate_weights = FALSE,
                     CC_weight_var,
                     alpha_level = 0.05)
{
  # Use a new data frame to avoid modifying the input data frame.
  df = data
  # Add unique subject id
  df = df %>%
    ungroup() %>%
    mutate(subject_id = row_number())
  # Recode the trial variable as a factor and add a new variable indicating the 
  # target trial
  df = df %>%
    ungroup() %>%
    mutate(
      trial = as.factor(.data[[trial_var]]),
      target_population = ifelse(.data[[trial_var]] == target_trial, TRUE, FALSE),
      treatment = .data[[treatment_var]]
    )
  
  # Does the analysis require case-cohort sampling?
  CC_sampling = estimate_weights | !is.null(CC_weight_var)
  
  # If the analysis requires case-cohort sampling, we need to check whether the 
  # weights need to be estimated are already given in the data.
  if(CC_sampling == TRUE){
    # Estimate weights if required.
    if (estimate_weights) {
      # From the formula_CC, extract the name of the case-cohort indicator
      # variable and create a variable Delta that indicates whether the
      # surrogate variable has been observed.
      outcome_var_CC = all.vars(formula_CC)[1]
      df <- df %>%
        mutate(Delta = .data[[outcome_var_CC]])
      # Fit logistic regression model for the case-cohort sampling mechanism and
      # use the fitted model to estimate the weights.
      CC_model_fit = glm(formula_CC, data = df, family = binomial)
      df <- df %>%
        mutate(CC_weight = 1 / predict(CC_model_fit, type = "response", newdata = pick(everything())))
    } else {
      df <- df %>%
        mutate(CC_weight = .data[[CC_weight_var]])
    }
    
    # Check for any estimated weights of Inf, which can occur if the fitted
    # model for the case-cohort sampling mechanism predicts a probability of
    # zero for some strata. 
    problematic_subjects = NULL
    if (any(df$CC_weight == Inf)) {
      # Select subjects with Inf weight.
      problematic_subjects = df %>%
        filter(CC_weight == Inf)
      n_problematic_subjects = nrow(problematic_subjects)
      warning(paste0(
        paste(
          c(
            "Estimated weight(s) of Inf for",
            n_problematic_subjects,
            "subject(s)"
          ),
          collapse = " "
        ),
        ". Corresponding observations are ignored for this analysis."
      ))
      df = df  %>%
        filter(CC_weight != Inf)
    }
  } else {
    # If case-cohort sampling is not used, we set the weights to 1 and Delta to
    # 1 for all subjects.
    df = df %>%
      mutate(CC_weight = 1,
             Delta = 1)
  }

  
  
  # standardization ---------------------------------------------------------
  
  # Fit the outcome models for the clinical (Y) and surrogate (S) endpoints.
  # Separate models are fit by `stratum_var` (e.g., trial) to allow for flexible
  # differences in the relationships between covariates and outcomes across
  # strata. The models for S are fit only among subjects with Delta == 1 and
  # using the case-cohort weights if case-cohort sampling is used.
  outcome_model_fits_df <- df %>%
    group_by(treatment) %>%
    group_split() %>%
    purrr::map(.x = ., .f = function(df_stratum) {
      # Fit the model for Y.
      outcome_model_fit_Y = glm(
        formula_Y,
        data = df_stratum,
        family = binomial,
        x = FALSE,
        y = FALSE
      )
      
      # Fit the model for S.
      outcome_model_fit_S = glm(
        formula_S,
        data = df_stratum,
        family = gaussian,
        subset = (Delta == 1),
        weights = CC_weight,
        x = FALSE,
        y = FALSE
      ) 
      
      # Return a data frame with the fitted models.
      tibble(
        treatment = unique(df_stratum$treatment),
        fit_Y = list(outcome_model_fit_Y),
        fit_S = list(outcome_model_fit_S),
        coef_Y = list(coef(outcome_model_fit_Y)),
        coef_S = list(coef(outcome_model_fit_S))
      )
    }) %>%
    bind_rows()

  
  # For each subject in the target trial, we predict the outcome using the model
  # estimated in any of the other trials/treatments. For each regression model
  # in `outcome_model_fits_df`, we thus create a new variable with the predicted
  # Y and S for all subjects from the target trial.
  df_target_population = df %>%
    filter(target_population)
  
  outcome_model_fits_df <- outcome_model_fits_df %>%
    mutate(
      predicted_Y = purrr::map(
        .x = fit_Y,
        .f = function(outcome_model_fit) {
          predict(outcome_model_fit, newdata = df_target_population, type = "response")
        }
      ),
      predicted_S = purrr::map(
        .x = fit_S,
        .f = function(outcome_model_fit) {
          predict(outcome_model_fit, newdata = df_target_population, type = "response")
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
    dplyr::select(treatment, mean_Y, mean_S)
  
  
  # IPW estimator ----------------------------------------------------------------
  
  # Estimate the treatment assignment probability, conditionally on X and trial.
  # Because we have simple randomization in each trial, the treatment assignment
  # probability is 0.5 for control in all trials, and 0.5 for the active
  # treatment in the trial in which the active treatment is given and 0 for
  # other treatments.
  
  # Create tibble in long format that contains for every subject (i.e., for X_i)
  # the treatment assignment probability P(A = a | X, T = t) for all a and t.
  df_weight_helper <- df %>%
    cross_join(expand_grid(treatment_pred = unique(df$treatment), trial_pred = unique(df$trial))) %>%
    mutate(predicted_prob_treatment_XT = case_when(
      treatment_pred == 0 ~ 0.5,
      treatment_pred >= 1 & trial_pred == treatment_pred ~ 0.5,
      .default = 0
    )) %>%
    # We only need to retain id, treatment_pred, and predicted_prob_treatment_XT
    # for the IPW estimator. 
    select(subject_id, treatment_pred, trial_pred, predicted_prob_treatment_XT) %>%
    as_tibble()
    
  
  # Estimate model for trial participation given X. The linear predictor of this
  # model is given in `formula_T`.
  trial_participation_model_fit <- nnet::multinom(formula_T, df)
  # Compute the predicted probabilities of trial participation for each subject
  # and add them to df_ipw.
  df_weight_helper <- df_weight_helper %>%
    left_join(
      fitted(trial_participation_model_fit) %>%
        as_tibble() %>%
        mutate(subject_id = df$subject_id) %>%
        pivot_longer(
          cols = -subject_id,
          names_to = "trial_pred",
          values_to = "trial_participation_prob_X"
        ),
      by = c("subject_id", "trial_pred"),
      relationship = "many-to-many"
    )
  
  # Compute the IPW weights for each subject in the target trial. These weights
  # are defined for treatment A = a and target trial T = t as P(T = t | X) /
  # (P(A = a | X) * P(T = t)) where P(A = a | X) := Sum_{t} P(A = a | T = t, X)
  # * P(T = t | X) and the probabilities in the denominator are estimated using
  # the models fitted above. The probabilities P(A = a | T = t, X) are given in
  # `df_ipw$predicted_prob_treatment` and the probabilities P(T = t | X) are
  # given in `df_ipw$trial_participation_prob_X`. The probability P(T = t) is
  # estimated as the proportion of subjects in the target trial.
  prob_target_trial = mean(df$target_population)
  
  # Compute P(A = a | X) for each subject in the target trial for all a. 
  df_weights <- df_weight_helper %>%
    group_by(subject_id, treatment_pred) %>%
    dplyr::summarize(predicted_prob_treatment_X = sum(predicted_prob_treatment_XT * trial_participation_prob_X))
  
  # Compute P(T = t | X) for t equal to the target trial.
  df_weights <- df_weights %>%
    left_join(
      df_weight_helper %>%
        filter(trial_pred == target_trial) %>%
        # There are multiple rows per subject, but the `trial_participation_prob_X`
        # values are identical. We thus take the first value for each subject.
        group_by(subject_id) %>%
        slice_max(n = 1, order_by = treatment_pred) %>%
        select(subject_id, trial_participation_prob_X),
      by = "subject_id",
      relationship = "many-to-one"
    )
  
  # Compute the overall weight.
  df_weights <- df_weights %>%
    mutate(weight = trial_participation_prob_X / (predicted_prob_treatment_X * prob_target_trial))
  
  # Compute the weighted means for Y and S for each treatment. .
  ipw_means_df <- df_weights %>%
    left_join(
      df %>%
        select(subject_id, treatment, Y, S),
      by = "subject_id",
      relationship = "many-to-one"
    ) %>%
    # Multiply missingness indicator for whether desired treatment does not match observed
    # treatment with the weight.
    mutate(W_hat_a = ifelse(treatment == treatment_pred, weight, 0)) %>%
    group_by(treatment) %>%
    dplyr::summarize(mean_Y = sum(Y * W_hat_a) / sum(W_hat_a),
                     mean_S = sum(S * W_hat_a) / sum(W_hat_a))
  
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
  
  # Because the naive means in the target trial's treatment group is a valid
  # estimator of that treatment group's mean standardized to the target trial,
  # we replace the standardized means with the corresponding naive means. This
  # is not necessary and this step could be omitted, however.
  standardized_means_df <- standardized_means_df %>%
    filter(treatment != target_trial) %>%
    bind_rows(naive_means_df %>%
                filter(trial == target_trial, vax == 1) %>%
                select(mean_Y, mean_S) %>%
                mutate(treatment = target_trial))
  
  ipw_means_df <- ipw_means_df %>%
    filter(treatment != target_trial) %>%
    bind_rows(naive_means_df %>%
                filter(trial == target_trial, vax == 1) %>%
                select(mean_Y, mean_S) %>%
                mutate(treatment = target_trial))
  
  # For computing the standardized treatment effects
  standardized_trt_effects_df = standardized_means_df %>%
    filter(treatment != 0) %>%
    rename(mean_Y_1 = mean_Y, mean_S_1 = mean_S) %>%
    cross_join(
      standardized_means_df %>%
        filter(treatment == 0) %>%
        rename(mean_Y_0 = mean_Y, mean_S_0 = mean_S) %>%
        select(-treatment)
    ) %>%
    mutate(VE_est = 1 - mean_Y_1 / mean_Y_0,
           mean_diff_S_est = mean_S_1 - mean_S_0)
  
  ipw_trt_effects_df = ipw_means_df %>%
    filter(treatment != 0) %>%
    rename(mean_Y_1 = mean_Y, mean_S_1 = mean_S) %>%
    cross_join(
      standardized_means_df %>%
        filter(treatment == 0) %>%
        rename(mean_Y_0 = mean_Y, mean_S_0 = mean_S) %>%
        select(-treatment)
    ) %>%
    mutate(VE_est = 1 - mean_Y_1 / mean_Y_0,
           mean_diff_S_est = mean_S_1 - mean_S_0)
  
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
        ipw_means_df = ipw_means_df,
        ipw_trt_effects_df = ipw_trt_effects_df,
        naive_means_df = naive_means_df,
        naive_trt_effects_df = naive_trt_effects_df,
        cor_standardized_df = cor_standardized_df,
        cor_naive_df = cor_naive_df,
        target_trial = target_trial,
        problematic_subjects = problematic_subjects
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
        target_trial = target_trial
        problematic_subjects = problematic_subjects
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


