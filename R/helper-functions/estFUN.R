estFUN_taudelta <- function(data, models_tbl, target_trial, estimate_weights, weights_df){
  n_trials = nrow(models_tbl)
  
  A = data$A
  if (A == 0) {
    trial_subject = "Placebo"
  } else {
    trial_subject = as.character(data$trial)
  }
  belongs_to_target_trial = trial_subject == target_trial
  trial_number_S = 0
  trial_number = 0
  placebo_number = which(models_tbl$trial_modified == "Placebo")
  
  if (!belongs_to_target_trial) {
    if (trial_subject != "Placebo") {
      trial_number_S = which(trial_subject == (
        models_tbl %>%
          filter(as.character(trial_modified) != "Placebo") %>%
          pull(trial_modified)
      ))
    }
    trial_number = which(trial_subject == as.character(models_tbl$trial_modified))
  }
  R = data$R
  Y = data$Y
  S = data$S
  S = ifelse(is.na(S), 0, S)
  weight_original = data$weight
  weight = weight_original
  CC_stratum = data$CC_stratum
  Delta = data$Delta
  
  ncov <- length(names(models_tbl$glm_fit_Y[[1]]$coefficients))
  covnames <- names(models_tbl$glm_fit_Y[[1]]$coefficients)[2:ncov]
  covs <- subset(data, select = covnames)
  X = as.matrix(cbind(1, covs))
  nX = ncol(X)
  
  # Estimating function for outcome regression model for the clinical outcome
  # for the patient's trial.
  clin_or_psiFUN_trial_A = if (belongs_to_target_trial) {
    function(theta_or) {
      return(rep(0, ncov))
    }
  } else {
    grab_psiFUN(models_tbl %>%
                  filter(as.character(trial_modified) == trial_subject) %>%
                  pull(glm_fit_Y) %>%
                  `[[`(1),
                data)
  }
  
  
  # Estimating function for outcome regression model, for the clinical outcome,
  # for all trial combined. Since each patient belongs to a unique trial, this
  # corresponds to a concatenation of estimating functions, all of which are
  # zero except the one corresponding to the trial to which the patient belongs.
  clin_or_psiFUN = function(theta_or){
    return_vec = rep(0, ncov * n_trials)
    
    if (trial_subject %in% target_trial) return(return_vec)
    
    # Determine position of elements corresponding to trial A
    trial_pos = (1:ncov) # position for trial on first row of `models_tbl`
    trial_pos = trial_pos + (trial_number - 1) * ncov
    
    # Set elements corresponding to trial A to the corresponding trial's
    # estimating function.
    # tryCatch({return_vec[trial_pos] = clin_or_psiFUN_trial_A(theta_or)}, error = function(e) browser())
    return_vec[trial_pos] = clin_or_psiFUN_trial_A(theta_or)
    
    return(return_vec)
  }
  
  # Determine position of parameters for the clinical outcome regression model.
  theta_pos_clin_or = 1:ncov
  theta_pos_clin_or = theta_pos_clin_or + (trial_number - 1) * ncov
  
  # Estimating function for outcome regression model for the surrogate outcome
  # for the patient's trial.
  surr_or_psiFUN_trial_A = if (Delta == 0 |
                               trial_subject %in% c("Placebo", target_trial)) {
    # If the surrogate is missing, the estimating function is zero.
    function(theta)
      return(rep(0, ncov))
  } else {
    grab_psiFUN(models_tbl %>%
                  filter(as.character(trial_modified) == trial_subject) %>%
                  pull(glm_fit_S) %>%
                  `[[`(1),
                data)
  }
  
  # Estimating function for outcome regression model, for the surrogate outcome,
  # for all trials combined. Since each patient belongs to a unique trial, this
  # corresponds to a concatenation of estimating functions, all of which are
  # zero except the one corresponding to the trial to which the patient belongs.
  surr_or_psiFUN = function(theta_or){
    return_vec = rep(0, ncov * (n_trials - 1))
    
    if (trial_subject %in% c("Placebo", target_trial)) return(return_vec)
    
    # Determine position of elements corresponding to trial A
    trial_pos = 1:ncov  # position for trial on first row of `models_tbl`
    trial_pos = trial_pos + (trial_number_S - 1) * ncov
    
    # Set elements corresponding to trial A to the corresponding trial's
    # estimating function.
    # browser()
    return_vec[trial_pos] = surr_or_psiFUN_trial_A(theta_or)
    
    return(return_vec)
  }
  
  # Determine position of parameters for the clinical outcome regression model.
  theta_pos_surr_or = NULL
  if (!(trial_subject %in% c("Placebo", target_trial))) {
    theta_pos_surr_or = 1:ncov + n_trials * ncov + (trial_number_S - 1) * ncov
  }
  
  # Matrix with in each column the positions of the regression parameters (for
  # the clinical outcome) where each column corresponds to a different trial.
  clin_or_parm_pos = matrix(1:(n_trials * nX), nrow = nX, ncol = n_trials, byrow = FALSE)
  
  # First position of the treatment effect parameters.
  trt_effect_start = (2 * n_trials * nX - nX) 
  
  # Position of parameters corresponding to the standardized mean clinical
  # outcomes.
  clin_stand_parm_pos = (trt_effect_start + 1):(trt_effect_start + 1 + n_trials)
  
  # Matrix with in each column the positions of the regression parameters (for
  # the surrogate outcome) where each column corresponds to a different trial.
  # Trial A == 1 or 0 are skipped because they correspond to the reference
  # population and the placebo group.
  surr_or_parm_pos = matrix((n_trials*nX + 1):((2 * n_trials - 1)*nX), nrow = nX, ncol = n_trials - 1, byrow = FALSE)
  
  # Position of parameters corresponding to the standardized mean surrogate
  # outcomes.
  surr_stand_parm_pos = (trt_effect_start + (1 + n_trials) + 1):(trt_effect_start  + (1 + n_trials) * 2 - 1)
  
  # Position of parameters corresponding to tau 1 to `n_trials`.
  tau_1_to_8_pos = (trt_effect_start + (1 + n_trials) * 2):(trt_effect_start + (1 + n_trials) * 3 - 2)
  
  # Position of correlation parameters (Pearson correlation and beta; in this order).
  corr_pos = (trt_effect_start + (1 + n_trials) * 3 - 1):(trt_effect_start + (1 + n_trials) * 3)
  
  # Number of weight strata
  n_CC_strata = weights_df %>%
    filter(CC_stratum != "Placebo") %>%
    nrow()
  
  # Position of the weight parameter in the theta vector. The theta-vector
  # actually contains the probability of being sampled parameters.
  CC_stratum_vec = weights_df %>%
    filter(CC_stratum != "Placebo") %>%
    pull(CC_stratum)
  weight_pos = which(CC_stratum == CC_stratum_vec) + corr_pos[2]
  
  function(theta, estimate_weights){
    # Determine subject-specific weight (which is a parameter in theta).
    if (estimate_weights) {
      if (CC_stratum != "Placebo") {
        weight = 1 / theta[weight_pos] %>% as.numeric()
      }
      else {
        weight = 1
      }
    }

    # Compute the predicted outcomes if this patient belongs to the target trial.
    if (R == 0) {
      # Predict Y given the models estimated in each trial (except the target trial)
      predicted_Y_all_trials = plogis(X %*% matrix(theta[clin_or_parm_pos], nrow = nX, byrow = FALSE))
      # Predict S given the models estimated in each trial (except the target trial)
      predicted_S_all_trials = X %*%  matrix(theta[surr_or_parm_pos], nrow = nX, byrow = FALSE)
    }
    else {
      # If a patient does not belong to the target population, their covariates
      # values do not contribute to the standardized estimates directly.
      predicted_Y_all_trials = rep(0, n_trials)
      predicted_S_all_trials = rep(0, n_trials - 1)
    }
    
    # Estimating equations for the standardized mean clinical outcome
    # parameters.
    clin_stand_parm_ee = c(
      (1 - R) * (A != 0) * (Y - theta[clin_stand_parm_pos[1]]),
      (1 - R) * (predicted_Y_all_trials - theta[clin_stand_parm_pos[-1]])
    )
    
    # Estimating equations for the standardized mean surrogate outcome
    # parameters.
    surr_stand_parm_ee = c(
      (1 - R) * (A != 0) * Delta * weight * (S - theta[surr_stand_parm_pos[1]]),
      (1 - R) * (predicted_S_all_trials - theta[surr_stand_parm_pos[-1]])
    )
    
    if (any(is.na(weight * surr_or_psiFUN(theta[theta_pos_surr_or])))) simpleError("NAs produced.")
    
    stacked_ee <- c(
      clin_or_psiFUN(theta[theta_pos_clin_or]),
      weight * surr_or_psiFUN(theta[theta_pos_surr_or]),
      clin_stand_parm_ee,
      surr_stand_parm_ee,
      (1 - theta[clin_stand_parm_pos[-c(placebo_number + 1)]] / theta[clin_stand_parm_pos[placebo_number + 1]]) - theta[tau_1_to_8_pos], #tau 1 to 8
      # causal association pearson rho
      cor(theta[tau_1_to_8_pos], theta[surr_stand_parm_pos]) - theta[corr_pos[1]],
      # causal association parameter beta -- linear regression slope
      (cov(theta[tau_1_to_8_pos], theta[surr_stand_parm_pos]) / 
        var(theta[surr_stand_parm_pos])) - theta[corr_pos[2]]
    )

    # If weights are treated as unknown, their corresponding estimating
    # equations are added to the set of stacked estimating equations.
    if (estimate_weights) {
      inv_weights_ee = rep(0, n_CC_strata)
      stacked_ee = c(stacked_ee, inv_weights_ee)
      if (CC_stratum != "Placebo") {
        stacked_ee[weight_pos] = Delta - theta[weight_pos]
      }
    }

    return(stacked_ee)
  }
}

