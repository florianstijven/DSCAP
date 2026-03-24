# Function to prepare the data and compute case-cohort weights if necessary
data_preparation <- function(data,
                             formula_CC,
                             trial_var = "trial",
                             treatment_var,
                             target_trial,
                             estimate_weights = FALSE,
                             CC_weight_var) {
  # If no target trial is specified, we set the target trial to 1 by default.
  if (is.null(target_trial)) {
    target_trial = 1
  }
  # Use a new data frame to avoid modifying the input data frame.
  df = data
  # Add unique subject ids to the data frame. These ids are needed further on to
  # merge various data frames.
  df = df %>%
    ungroup() %>%
    mutate(subject_id = row_number())
  # Recode the trial variable as a factor and add a new variable indicating the
  # target trial and the treatment variable. Further functions and code assume
  # the existence of the variables `trial`, `target_population`, and `treatment`
  # in the data frame.
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
  # weights need to be estimated or whether they are already given in the data.
  if (CC_sampling == TRUE) {
    # Estimate weights if required.
    if (estimate_weights) {
      # From the formula_CC, extract the name of the case-cohort indicator
      # variable and create a variable Delta that indicates whether the
      # surrogate variable has been observed.
      outcome_var_CC = all.vars(formula_CC)[1]
      df <- df %>%
        mutate(Delta = .data[[outcome_var_CC]])
      # Fit logistic regression model for the case-cohort sampling mechanism and
      # use the fitted model to estimate the weights. Note that the model is not
      # automatically stratified by trial in the current implementation.
      CC_model_fit = glm(
        formula_CC,
        data = df,
        family = binomial,
        x = FALSE,
        y = FALSE
      )
      df <- df %>%
        mutate(CC_weight = 1 / predict(
          CC_model_fit,
          type = "response",
          newdata = pick(everything())
        ))
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
      mutate(CC_weight = 1, Delta = 1)
  }
  
  return(df)
}

# Standardization estimator
standardization_estimator <- function(data,
                                      formula_Y,
                                      formula_S,
                                      formula_CC,
                                      trial_var = "trial",
                                      treatment_var,
                                      target_trial,
                                      estimate_weights = FALSE,
                                      CC_weight_var,
                                      alpha = 0.05) {
  df = data_preparation(
    data,
    formula_CC,
    trial_var = "trial",
    treatment_var,
    target_trial,
    estimate_weights,
    CC_weight_var
  )
  
  # Does the analysis require case-cohort sampling?
  CC_sampling = estimate_weights | !is.null(CC_weight_var)
  
  # Fit the outcome models for the clinical (Y) and surrogate (S) endpoints.
  # Separate models are fit by trial to allow for flexible differences in the
  # relationships between covariates and outcomes across strata. The models for
  # S are fit only among subjects with `Delta` == 1 and using the case-cohort
  # weights `CC_weight` if case-cohort sampling is used.
  outcome_model_fits_df <- fit_outcome_models(df, formula_Y, formula_S)
  
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
      mean_Y = purrr::map_dbl(.x = predicted_Y, .f = mean),
      mean_S = purrr::map_dbl(.x = predicted_S, .f = mean)
    ) %>%
    dplyr::select(treatment, mean_Y, mean_S)
  
  return(standardized_means_df)
  
}

# IPW estimator
ipw_estimator <- function(data,
                          formula_T,
                          formula_CC,
                          trial_var = "trial",
                          treatment_var,
                          target_trial,
                          estimate_weights = FALSE,
                          CC_weight_var,
                          alpha = 0.05) {
  df = data_preparation(
    data,
    formula_CC,
    trial_var = "trial",
    treatment_var,
    target_trial,
    estimate_weights,
    CC_weight_var
  )
  
  # Does the analysis require case-cohort sampling?
  CC_sampling = estimate_weights | !is.null(CC_weight_var)
  
  # Compute the IPW weights for each subject in the target trial.
  df_weights = compute_ipw_weights(df, formula_T, target_trial)
  
  # Compute the weighted means for Y and S for each treatment.
  ipw_means_df <- df_weights %>%
    left_join(df %>%
                select(subject_id, treatment, Y, S, Delta, CC_weight),
              by = "subject_id",
              relationship = "many-to-one") %>%
    # Multiply missingness indicator for whether desired treatment does not match observed
    # treatment with the weight.
    mutate(W_hat_a = ifelse(treatment == treatment_pred, weight, 0)) %>%
    group_by(treatment) %>%
    mutate(S = ifelse(Delta == 1, S, 0)) %>%
    # Divide the weighted sum of the outcomes by the sum of the weights. This is
    # asymptotically equivalent to dividing by the sample size, but generally
    # leads to better finite-sample performance.
    dplyr::summarize(
      mean_Y = sum(Y * W_hat_a) / sum(W_hat_a),
      mean_S = sum(S * W_hat_a * Delta * CC_weight) / sum(W_hat_a * Delta * CC_weight)
    )
  
  return(ipw_means_df)
}

# Doubly robust estimator
DR_estimator <- function(data,
                         formula_Y,
                         formula_S,
                         formula_T,
                         formula_CC,
                         trial_var = "trial",
                         treatment_var,
                         target_trial,
                         estimate_weights = FALSE,
                         CC_weight_var,
                         alpha = 0.05) {
  df = data_preparation(
    data,
    formula_CC,
    trial_var = "trial",
    treatment_var,
    target_trial,
    estimate_weights,
    CC_weight_var
  )
  
  # Does the analysis require case-cohort sampling?
  CC_sampling = estimate_weights | !is.null(CC_weight_var)
  
  # Fit the outcome models for the clinical (Y) and surrogate (S) endpoints.
  # Separate models are fit by trial to allow for flexible differences in the
  # relationships between covariates and outcomes across strata. The models for
  # S are fit only among subjects with `Delta` == 1 and using the case-cohort
  # weights `CC_weight` if case-cohort sampling is used.
  outcome_model_fits_df <- fit_outcome_models(df, formula_Y, formula_S)
  
  
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
  # These estimates will be corrected below using the IPW weights to get the
  # doubly-robust estimates.
  standardized_means_df = outcome_model_fits_df %>%
    mutate(
      mean_Y = purrr::map_dbl(.x = predicted_Y, .f = mean),
      mean_S = purrr::map_dbl(.x = predicted_S, .f = mean)
    ) %>%
    dplyr::select(treatment, mean_Y, mean_S)
  
  # For every subject (in all trials), we need the estimate of E[Y | A = a, X]
  # and E[S | A = a, X] for all a. We can obtain these estimates from the models
  # fitted above by predicting the outcomes for all subjects in the data, but
  # setting the treatment variable to a specific value a. We thus create new
  # variables with the predicted Y and S for all subjects in the data for each
  # treatment a.
  all_predictions_df <- outcome_model_fits_df %>%
    group_by(treatment) %>%
    group_split() %>%
    purrr::map(
      .x = .,
      .f = function(df_stratum) {
        # For each fitted model, we predict the outcomes for all subjects in the
        # data, but setting the treatment variable to a specific value a.
        predicted_Y_a = predict(df_stratum$fit_Y[[1]], newdata = df, type = "response")
        predicted_S_a = predict(df_stratum$fit_S[[1]], newdata = df, type = "response")
        
        # Return a data frame with the predicted outcomes.
        tibble(
          subject_id = df$subject_id,
          treatment = df$treatment,
          treatment_pred = unique(df_stratum$treatment),
          Y = df$Y,
          predicted_Y_a = predicted_Y_a,
          S = df$S,
          predicted_S_a = predicted_S_a
        )
      }
    ) %>%
    bind_rows()
  
  # Compute the IPW weights for each subject in the target trial.
  df_weights = compute_ipw_weights(df, formula_T, target_trial)
  
  # Compute correction terms to get the doubly-robust estimates of the means for
  # Y and S for each treatment.

  DR_corrections_df <- all_predictions_df %>%
    left_join(
      df_weights %>%
        select(subject_id, treatment_pred, weight),
      by = c("subject_id", "treatment_pred"),
      relationship = "one-to-one"
    ) %>%
    left_join(
      df %>%
        select(subject_id, Delta, CC_weight),
      by = "subject_id",
      relationship = "many-to-one"
    ) %>%
    mutate(W_hat_a = ifelse(treatment == treatment_pred, weight, 0)) %>%
    mutate(S = ifelse(Delta == 1, S, 0)) %>%
    group_by(treatment) %>%
    dplyr::summarize(
      mean_Y_correction = sum(W_hat_a * (Y - predicted_Y_a) * Delta * CC_weight) / sum(W_hat_a * Delta * CC_weight),
      mean_S_correction = sum(W_hat_a * (S - predicted_S_a) * Delta * CC_weight) / sum(W_hat_a * Delta * CC_weight)
    )
  
  # Add the correction terms to the standardized means to get the doubly-robust
  # estimates of the means for Y and S (standardized to the target population)
  # for each treatment.
  DR_means_df <- DR_corrections_df %>%
    left_join(standardized_means_df,
              by = "treatment",
              relationship = "many-to-one") %>%
    mutate(mean_Y = mean_Y + mean_Y_correction,
           mean_S = mean_S + mean_S_correction) %>%
    select(treatment, mean_Y, mean_S)
  
  return(DR_means_df)
}

naive_estimator <- function(data,
                            formula_CC,
                            trial_var = "trial",
                            treatment_var,
                            target_trial,
                            estimate_weights = FALSE,
                            CC_weight_var,
                            alpha = 0.05) {
  df = data_preparation(
    data,
    formula_CC,
    trial_var = "trial",
    treatment_var,
    target_trial,
    estimate_weights,
    CC_weight_var
  )
  naive_means_df = df %>%
    mutate(S = ifelse(Delta == 1, S, 0)) %>%
    group_by(trial, treatment) %>%
    dplyr::summarize(
      mean_Y = mean(Y),
      n = n(),
      mean_S = sum(Delta * CC_weight * S) / sum(Delta * CC_weight),
      var_S = var(Delta * CC_weight * S),
      .groups = "drop"
    )
  
  return(naive_means_df)
}

# Compute treatment effects from the estimated means in each trial.
trt_effects <- function(means_df, by_trial) {
  # Compute the treatment effects for each trial. The `by_trial` argument
  # indicates whether one control group mean should be used for all trials
  # (i.e., for the standardized means estimator) or whether trial-specific
  # control group means should be used (i.e., for the naive estimator).
  if (by_trial) {
    # We assume here that the control group (indicated by treatment = 0) is
    # present in all trials and that the trials have one unique active treatment
    # group.
    means_df = means_df %>%
      select(trial, treatment, mean_Y, mean_S) %>%
      filter(treatment == 0) %>%
      rename(mean_Y_0 = mean_Y, mean_S_0 = mean_S) %>%
      select(-treatment) %>%
      left_join(
        means_df %>%
          select(trial, treatment, mean_Y, mean_S) %>%
          filter(treatment != 0) %>%
          rename(mean_Y_1 = mean_Y, mean_S_1 = mean_S),
        by = "trial",
        relationship = "one-to-many"
      ) %>%
      select(-trial)
  } else {
    # We assume here that the control group (indicated by treatment = 0) is
    # present in all trials and that the trials have one unique active treatment
    # group. We also assume that the control group means are identical across
    # trials.
    means_df = means_df %>%
      select(treatment, mean_Y, mean_S) %>%
      filter(treatment == 0) %>%
      rename(mean_Y_0 = mean_Y, mean_S_0 = mean_S) %>%
      select(-treatment) %>%
      left_join(
        means_df %>%
          select(treatment, mean_Y, mean_S) %>%
          filter(treatment != 0) %>%
          rename(mean_Y_1 = mean_Y, mean_S_1 = mean_S),
        by = character(),
        relationship = "one-to-many"
      )
  }
  trt_effects_df = means_df %>%
    mutate(VE_est = 1 - (mean_Y_1 / mean_Y_0),
           mean_diff_S_est = mean_S_1 - mean_S_0) %>%
    select(treatment, VE_est, mean_diff_S_est)
  
  return(trt_effects_df)
}

# Function to fit the outcome models needed in the standardization and
# doubly-robust estimators. Note that the `df` argument should be a data frame
# returned by `data_preparation()`.
fit_outcome_models <- function(df, formula_Y, formula_S) {
  # Fit the outcome models for the clinical (Y) and surrogate (S) endpoints.
  # Separate models are fit by `stratum_var` (e.g., trial) to allow for flexible
  # differences in the relationships between covariates and outcomes across
  # strata. The models for S are fit only among subjects with Delta == 1 and
  # using the case-cohort weights if case-cohort sampling is used.
  outcome_model_fits_df <- df %>%
    group_by(treatment) %>%
    group_split() %>%
    purrr::map(
      .x = .,
      .f = function(df_stratum) {
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
      }
    ) %>%
    bind_rows()
  
  return(outcome_model_fits_df)
}

# Compute the IPW weights for each subject to standardize treatment groups means
# to the target population.
compute_ipw_weights <- function(df, formula_T, target_trial) {
  # Estimate the treatment assignment probability, conditionally on X and trial.
  # Because we have simple randomization in each trial, the treatment assignment
  # probability is 0.5 for control in all trials, and 0.5 for the active
  # treatment in the trial in which the active treatment is given and 0 for
  # other treatments.

  # Create tibble in long format that contains for every subject (i.e., for X_i)
  # the treatment assignment probability P(A = a | X, T = t) for all a and t.
  # Because we assume simple randomized, we need to estimate P(A = a | T = t)
  # for all a, t.
  treatment_prob_by_trial_tbl = df %>%
    cross_join(expand_grid(
      treatment_pred = unique(df$treatment),
    )) %>%
    group_by(trial, treatment_pred) %>%
    dplyr::summarize(predicted_prob_treatment_XT = mean(treatment == treatment_pred)) %>%
    rename(trial_pred = trial) %>%
    ungroup()
  

  df_weight_helper <- df %>%
    cross_join(expand_grid(
      treatment_pred = unique(df$treatment),
      trial_pred = unique(df$trial)
    )) %>%
    left_join(
      treatment_prob_by_trial_tbl,
      by = c("trial_pred", "treatment_pred"),
      relationship = "many-to-one"
    ) %>%
    # We only need to retain id, treatment_pred, and predicted_prob_treatment_XT
    # for the IPW estimator.
    select(subject_id,
           treatment_pred,
           trial_pred,
           predicted_prob_treatment_XT) %>%
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
  
  df_weights <- df_weight_helper %>%
    group_by(subject_id, treatment_pred) %>%
    dplyr::summarize(
      predicted_prob_treatment_X = sum(predicted_prob_treatment_XT * trial_participation_prob_X),
      target_trial_participation_prob_X = trial_participation_prob_X[trial_pred == target_trial],
      weight = target_trial_participation_prob_X / (predicted_prob_treatment_X * prob_target_trial)
    )
  
  return(df_weights)
}
