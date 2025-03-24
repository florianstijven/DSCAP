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

estFUN_taudelta.fh.8trial <- function(data, models){
  # print(data$Ptid)
  
  A = data$A
  R = data$R
  Y = data$Y
  S = data$S
  W = data$weight
  Delta = data$Delta
  
  covnames <- names(models[[1]]$coefficients)[2:length(names(models[[1]]$coefficients))]
  ncov <- length(names(models[[1]]$coefficients))
  covs <- subset(data, select=covnames)
  X = as.matrix(cbind(1,covs))
  nX = ncol(X)
  
  # Estimating function for outcome regression model for the clinical outcome
  # for the patient's trial.
  clin_or_psiFUN_trial_A = if (A == 0) {
    grab_psiFUN(models[[1]], data)
  } else if (A == 1) {
    function(theta_or) {
      return(rep(0, ncov))
    }
  }
    else {
    grab_psiFUN(models[[A]], data)
  }
  
  # Estimating function for outcome regression model, for the clinical outcome,
  # for all trial combined. Since each patient belongs to a unique trial, this
  # corresponds to a concatenation of estimating functions, all of which are
  # zero except the one corresponding to the trial to which the patient belongs.
  clin_or_psiFUN = function(theta_or){
    return_vec = rep(0, ncov * 8)
    
    # Determine position of elements corresponding to trial A
    trial_pos = (1:ncov) # position for trial A == 0
    if (A != 0) trial_pos = trial_pos + (A - 1) * ncov
    
    # Set elements corresponding to trial A to the corresponding trial's
    # estimating function.
    return_vec[trial_pos] = clin_or_psiFUN_trial_A(theta_or)
    
    return(return_vec)
  }
  
  # Determine position of parameters for the clinical outcome regression model.
  theta_pos_clin_or = 1:ncov
  if (A != 0) theta_pos_clin_or = theta_pos_clin_or + (A - 1) * ncov
  
  # Estimating function for outcome regression model for the surrogate outcome
  # for the patient's trial.
  surr_or_psiFUN_trial_A = if (Delta == 0 | A <= 1) {
    # If the surrogate is missing, the estimating function is zero.
    function(theta) return(rep(0, ncov))
  } else {
    grab_psiFUN(models[[A + 8]], data)
  }
  
  # Estimating function for outcome regression model, for the surrogate outcome,
  # for all trials combined. Since each patient belongs to a unique trial, this
  # corresponds to a concatenation of estimating functions, all of which are
  # zero except the one corresponding to the trial to which the patient belongs.
  surr_or_psiFUN = function(theta_or){
    return_vec = rep(0, ncov * 7)
    
    if (A <= 1) return(return_vec)
    
    # Determine position of elements corresponding to trial A
    trial_pos = 1:ncov # position for trial A == 2
    trial_pos = trial_pos + max(c(0, A - 2)) * ncov
    
    # Set elements corresponding to trial A to the corresponding trial's
    # estimating function.
    # if (length(return_vec[trial_pos]) != length(surr_or_psiFUN_trial_A(theta_or)))  browser()
    return_vec[trial_pos] = surr_or_psiFUN_trial_A(theta_or)
    
    
    return(return_vec)
  }
  
  # Determine position of parameters for the clinical outcome regression model.
  theta_pos_surr_or = NULL
  if (A >= 2) {
    theta_pos_surr_or = 1:ncov + 8 * ncov + (A - 2) * ncov
  }
  
  # weights
  # NOTE: 1 weights category has 100% coverage, nZ-1 is the proper
  # number of categories

  # position of the weight parameter in the theta vector
  Zind <- data$Zind2
  wtvec <- rep(0,(data$nZ))
  
  py1 <- ifelse(A == 1, Y, 0)
  ps0 <- data$ps0
  ps1 <- ifelse(A == 1, S, 0)
  
  # Matrix with in each column the positions of the regression parameters (for
  # the clinical outcome) where each column corresponds to a different trial.
  # Trial A == 1 is skipped because this corresponds to the reference
  # population.
  clin_or_parm_pos = matrix(1:(8 * nX), nrow = nX, ncol = 8, byrow = FALSE)
  
  # Position of parameters corresponding to the standardized mean clinical
  # outcomes.
  clin_stand_parm_pos = (15*nX + 1):(15*nX + 9)
  
  # Matrix with in each column the positions of the regression parameters (for
  # the surrogate outcome) where each column corresponds to a different trial.
  # Trial A == 1 or 0 are skipped because they correspond to the reference
  # population and the placebo group.
  surr_or_parm_pos = matrix((8*nX + 1):(15*nX), nrow = nX, ncol = 7, byrow = FALSE)
  
  # Position of parameters corresponding to the standardized mean surrogate
  # outcomes.
  surr_stand_parm_pos = (15*nX + 11):(15*nX + 17)
  
  # Position of parameters corresponding to tau 1 to 8.
  tau_1_to_8_pos = (15*nX + 18):(15*nX + 25)

  function(theta){
    # py0 <- plogis(X %*% theta[1:nX])
    # # py1 <- ifelse(A == 1, Y, 0)
    # py2 <- plogis(X %*% theta[(nX+1):(2*nX)])
    # py3 <- plogis(X %*% theta[(2*nX + 1):(3*nX)])
    # py4 <- plogis(X %*% theta[(3*nX + 1):(4*nX)])
    # py5 <- plogis(X %*% theta[(4*nX + 1):(5*nX)])
    # py6 <- plogis(X %*% theta[(5*nX + 1):(6*nX)])
    # py7 <- plogis(X %*% theta[(6*nX + 1):(7*nX)])
    # py8 <- plogis(X %*% theta[(7*nX + 1):(8*nX)])
    
    # The following code computes the above set of values more efficiently.
    py0_to_8 = plogis(X %*% matrix(theta[clin_or_parm_pos], nrow = nX, byrow = FALSE))
    py0_to_8 = c(py0_to_8[1], Y, py0_to_8[2:8])
    
    
    # #ps0 <- X %*% theta[16:18]
    # # ps0 <- data$ps0
    # # ps1 <- ifelse(A == 1, S, 0)
    # ps2 <- X %*% theta[(8*nX + 1):(9*nX)]
    # ps3 <- X %*% theta[(9*nX + 1):(10*nX)]
    # ps4 <- X %*% theta[(10*nX + 1):(11*nX)]
    # ps5 <- X %*% theta[(11*nX + 1):(12*nX)]
    # ps6 <- X %*% theta[(12*nX + 1):(13*nX)]
    # ps7 <- X %*% theta[(13*nX + 1):(14*nX)]
    # ps8 <- X %*% theta[(14*nX + 1):(15*nX)]
    
    # The following code computes the above set of values more efficiently.
    ps2_to_8 = X %*%  matrix(theta[surr_or_parm_pos], nrow = nX, byrow = FALSE)
    
    # Estimating equations for the standardized mean clinical outcome
    # parameters.
    clin_stand_parm_ee = py0_to_8 - theta[clin_stand_parm_pos]
    clin_stand_parm_ee[2] = clin_stand_parm_ee[2] * (A == 1)
    # browser()
    
    rtrn <- c(
      clin_or_psiFUN(theta[theta_pos_clin_or]),
      surr_or_psiFUN(theta[theta_pos_surr_or]),
      clin_stand_parm_ee,
      if(A == 1 & Delta == 1) (S - theta[15*nX + 10])*W else 0,
      (1-R)*(ps2_to_8 - theta[surr_stand_parm_pos]),
      (1 - theta[clin_stand_parm_pos[-1]] / theta[clin_stand_parm_pos[1]]) - theta[tau_1_to_8_pos], #tau 1 to 8
      # causal association pearson rho
      # cov(c(theta[(15*nX + 18):(15*nX + 25)]),
      #     c(theta[(15*nX + 10):(15*nX + 17)])) / 
      #   ( sd(c(theta[(15*nX + 18):(15*nX + 25)])) * 
      #       sd(c(theta[(15*nX + 10):(15*nX + 17)]))) - 
      #   theta[15*nX + 26],
      cor(theta[tau_1_to_8_pos], theta[surr_or_parm_pos]) - theta[15*nX + 26],
      # causal association parameter beta -- linear regression slope
      cov(c(theta[(15*nX + 18):(15*nX + 25)]),
          c(theta[(15*nX + 10):(15*nX + 17)]-ps0)) / 
        var(c(theta[(15*nX + 10):(15*nX + 17)]))  -
        theta[15*nX + 27]
    )
    

    wtvec[Zind] <- ifelse(!(W == 1), Delta - theta[102 + Zind], 0)
    # browser()
    return(c(rtrn, wtvec))
  }
}


estFUN_taudelta.fh.6trial <- function(data){
  A = data$A
  R = data$R
  Y = data$Y
  S = data$S
  Delta = data$Delta

  function(theta, models){
    
    covnames <- names(models[[1]]$coefficients)[2:length(names(models[[1]]$coefficients))]
    ncov <- length(names(models[[1]]$coefficients))
    covs <- subset(data, select=covnames)
    X = as.matrix(cbind(1,covs))
    
    py0 <- plogis(X %*% theta[1:6])
    py1 <- ifelse(data$A==1, data$Y, 0)
    py2 <- plogis(X %*% theta[7:12])
    py3 <- plogis(X %*% theta[13:18])
    py4 <- plogis(X %*% theta[19:24])
    py5 <- plogis(X %*% theta[25:30])
    py6 <- plogis(X %*% theta[31:36])
    
    #ps0 <- X %*% theta[16:18]
    ps0 <- data$ps0
    ps1 <- ifelse(A==1, data$S, 0)
    ps2 <- X %*% theta[37:42]
    ps3 <- X %*% theta[43:48]
    ps4 <- X %*% theta[49:54]
    ps5 <- X %*% theta[67:60]
    ps6 <- X %*% theta[73:66]
    
    rtrn <- c(
      if(A==0) grab_psiFUN(models[[1]], data)(theta[1:6]) else rep(0, ncov),
      #if(A==1) (data$Y - theta[4]) else 0,
      if(A==2) grab_psiFUN(models[[2]], data)(theta[7:12]) else rep(0, ncov),
      if(A==3) grab_psiFUN(models[[3]], data)(theta[13:18]) else rep(0, ncov),
      if(A==4) grab_psiFUN(models[[4]], data)(theta[19:24]) else rep(0, ncov),
      if(A==5) grab_psiFUN(models[[5]], data)(theta[25:30]) else rep(0, ncov),
      if(A==6) grab_psiFUN(models[[6]], data)(theta[31:36]) else rep(0, ncov),
      
      #if(A==0) grab_psiFUN(models[[7]], data)(theta[16:18]) else rep(0, 3),
      #if(A==1) (data$S - theta[20]) else 0,
      if(A==2 & Delta==1) grab_psiFUN(models[[8]], data)(theta[37:42]) else rep(0, ncov),
      if(A==3 & Delta==1) grab_psiFUN(models[[9]], data)(theta[43:48]) else rep(0, ncov),
      if(A==4 & Delta==1) grab_psiFUN(models[[10]], data)(theta[49:54]) else rep(0, ncov),
      if(A==5 & Delta==1) grab_psiFUN(models[[11]], data)(theta[55:60]) else rep(0, ncov),
      if(A==6 & Delta==1) grab_psiFUN(models[[12]], data)(theta[61:66]) else rep(0, ncov),
      
      (1-R)*(py0 - theta[67]),
      if(A==1) (data$Y - theta[68]) else 0,
      (1-R)*(py2 - theta[69]),
      (1-R)*(py3 - theta[70]),
      (1-R)*(py4 - theta[71]),
      (1-R)*(py5 - theta[72]),
      (1-R)*(py6 - theta[73]),
      
      #(1-R)*(ps0 - theta[37]),
      if(A==1 & Delta==1) (data$S - theta[74])*data$weight else 0,
      (1-R)*(ps2 - theta[75]),
      (1-R)*(ps3 - theta[76]),
      (1-R)*(ps4 - theta[77]),
      (1-R)*(ps5 - theta[78]),
      (1-R)*(ps6 - theta[79]),
      
      (1 - theta[68] / theta[67]) - theta[80], #tau 1
      (1 - theta[69] / theta[67]) - theta[81], #tau 2
      (1 - theta[70] / theta[67]) - theta[82], #tau 3
      (1 - theta[71] / theta[67]) - theta[83], #tau 4
      (1 - theta[72] / theta[67]) - theta[84], #tau 5
      (1 - theta[73] / theta[67]) - theta[85], #tau 6
      
      
      # Delta and psM are collinear with ps0 as a constant
      # (theta[38] - theta[37]) - theta[48], #delta1
      # (theta[39] - theta[37]) - theta[49], #delta2
      # (theta[40] - theta[37]) - theta[50], #delta3
      # (theta[41] - theta[37]) - theta[51], #delta4
      # (theta[42] - theta[37]) - theta[52], #delta5
      
      
      # causal association parameter beta -- linear regression slope
      cov(c(theta[80:85]),
          c(theta[74:79]-ps0)) / 
        var(c(theta[74:79]))  -
        theta[81],
      
      # causal association pearson rho
      cov(c(theta[80:85]),
          c(theta[75:79])) / 
        ( sd(c(theta[80:85])) * 
            sd(c(theta[75:79]))) - 
        theta[80])
    
    # weights
    # NOTE: 1 weights category has 100% coverage, nZ-1 is the proper
    # number of categories
    Zind <- ifelse(data$Zind <= 33, data$Zind, data$Zind-1)
    wtvec <- rep(0,(data$nZ-1))
    wtvec[Zind] <- ifelse(!(data$weight==1),
                          Delta-theta[80+Zind],0)
    
    return(c(rtrn,wtvec))
  }
}

estFUN_taudelta.fh.5trial <- function(data){
  A = data$A
  R = data$R
  Y = data$Y
  S = data$S
  Delta = data$Delta
  
  function(theta, models){
    
    covnames <- names(models[[1]]$coefficients)[2:length(names(models[[1]]$coefficients))]
    ncov <- length(names(models[[1]]$coefficients))
    covs <- subset(data, select=covnames)
    X = as.matrix(cbind(1,covs))
    nX = ncol(X)
    
    py0 <- plogis(X %*% theta[1:nX])
    py1 <- ifelse(data$A==1, data$Y, 0)
    py2 <- plogis(X %*% theta[(nX+1):(2*nX)])
    py3 <- plogis(X %*% theta[(2*nX + 1):(3*nX)])
    py4 <- plogis(X %*% theta[(3*nX + 1):(4*nX)])
    py5 <- plogis(X %*% theta[(4*nX + 1):(5*nX)])
    
    ps0 <- data$ps0
    ps1 <- ifelse(A==1, data$S, 0)
    ps2 <- X %*% theta[(5*nX + 1):(6*nX)]
    ps3 <- X %*% theta[(6*nX + 1):(7*nX)]
    ps4 <- X %*% theta[(7*nX + 1):(8*nX)]
    ps5 <- X %*% theta[(8*nX + 1):(9*nX)]
  
    rtrn <- c(
      if(A==0) grab_psiFUN(models[[1]], data)(theta[1:nX]) else rep(0, ncov),
      if(A==2) grab_psiFUN(models[[2]], data)(theta[(nX+1):(2*nX)]) else rep(0, ncov),
      if(A==3) grab_psiFUN(models[[3]], data)(theta[(2*nX + 1):(3*nX)]) else rep(0, ncov),
      if(A==4) grab_psiFUN(models[[4]], data)(theta[(3*nX + 1):(4*nX)]) else rep(0, ncov),
      if(A==5) grab_psiFUN(models[[5]], data)(theta[(4*nX + 1):(5*nX)]) else rep(0, ncov),
     

      if(A==2 & Delta==1) grab_psiFUN(models[[7]], data)(theta[(5*nX + 1):(6*nX)]) else rep(0, ncov),
      if(A==3 & Delta==1) grab_psiFUN(models[[8]], data)(theta[(6*nX + 1):(7*nX)]) else rep(0, ncov),
      if(A==4 & Delta==1) grab_psiFUN(models[[9]], data)(theta[(7*nX + 1):(8*nX)]) else rep(0, ncov),
      if(A==5 & Delta==1) grab_psiFUN(models[[10]], data)(theta[(8*nX + 1):(9*nX)]) else rep(0, ncov),
      
      (1-R)*(py0 - theta[9*nX + 1]),
      if(A==1) (data$Y - theta[9*nX + 2]) else 0,
      (1-R)*(py2 - theta[9*nX + 3]),
      (1-R)*(py3 - theta[9*nX + 4]),
      (1-R)*(py4 - theta[9*nX + 5]),
      (1-R)*(py5 - theta[9*nX + 6]),
      
      if(A==1 & Delta==1) (data$S - theta[9*nX + 7])*data$weight else 0,
      (1-R)*(ps2 - theta[9*nX + 8]),
      (1-R)*(ps3 - theta[9*nX + 9]),
      (1-R)*(ps4 - theta[9*nX + 10]),
      (1-R)*(ps5 - theta[9*nX + 11]),
 
      
      (1 - theta[9*nX + 2] / theta[9*nX + 1]) - theta[9*nX + 12], #tau 1
      (1 - theta[9*nX + 3] / theta[9*nX + 1]) - theta[9*nX + 13], #tau 2
      (1 - theta[9*nX + 4] / theta[9*nX + 1]) - theta[9*nX + 14], #tau 3
      (1 - theta[9*nX + 5] / theta[9*nX + 1]) - theta[9*nX + 15], #tau 4
      (1 - theta[9*nX + 6] / theta[9*nX + 1]) - theta[9*nX + 16], #tau 5
      
      
      # causal association pearson rho
      cov(c(theta[(9*nX + 12):(9*nX + 16)]),
          c(theta[(9*nX + 7):(9*nX + 11)])) / 
        ( sd(c(theta[(9*nX + 12):(9*nX + 16)])) * 
            sd(c(theta[(9*nX + 7):(9*nX + 11)]))) - 
        theta[9*nX + 17],
      
      # causal association parameter beta -- linear regression slope
      cov(c(theta[(9*nX + 12):(9*nX + 16)]),
          c(theta[(9*nX + 7):(9*nX + 11)]-ps0)) / 
        var(c(theta[(9*nX + 7):(9*nX + 11)]))  -
        theta[9*nX + 18]
    )

    # weights
    Zind <- data$Zind2
    wtvec <- rep(0,data$nZ)
    if (data$weight!=1) {
      wtvec[Zind] <- Delta-theta[(9*nX + 18)+Zind]
    }
    
    return(c(rtrn,wtvec))
    
  }
}
