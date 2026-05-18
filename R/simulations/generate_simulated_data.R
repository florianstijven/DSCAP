#######################################################
# Programmer: K. Cross (adapted from code by S. Rosin) 
# Puropose: Generate a dataset 
# Inputs: n_trials = number of trials
#         n_t = number of observations per trial (must be same for all trials)
#         Ymod = Model for Y outcome
#         Smod = Model for S outcome
#         gamma.set = list of settings for gamma (Y models)
#         theta.set = list of settings for theta (S models)
# Output: dataset of simulated data
#######################################################
generate_simulated_data <- function(
    n_trials,
    n_t,
    formula_Y,
    formula_S,
    gamma,
    theta,
    zeta,
    CC_sampling = FALSE,
    target_trial = NULL,
    different_placebo_model = FALSE,
    bridge_error = FALSE
) {
  
  N <- n_trials * n_t
  
  inv_logit <- function(x) plogis(x)
  
  # Generate covariates and trial indicators. 
  if (is.null(target_trial)) {
    # Sample covariates and trial assignment according to the "default" DGP.
     sample_df <- sample_covariates_trials(zeta = zeta, N = N, target_trial = NULL)
  } else {
    # Sample covariates from the target trial covariate distribution and sample the trial indicator
    # from a uniform distribution on the trial indicators.
    sample_df <- sample_covariates_trials(zeta = zeta, N = N, target_trial = target_trial)
  }
  
  # Sample treatment indicators.
  sample_df$vax <- rbinom(N, 1, 0.5)
  sample_df$A   <- ifelse(sample_df$vax == 1, sample_df$trial, 0)
  
  sample_df$Y <- -1
  sample_df$S <- -1
  
  # ----------------------------
  # Generate potential outcomes Treatments are 0,1,...,n_trials
  # ----------------------------
  # Make deisgn matrix based on formula_Y and formula_S. 
  model_matrix_Y = model.matrix(formula_Y, data = sample_df)
  model_matrix_S = model.matrix(formula_S, data = sample_df)
  
  # If `different_placebo_model` is TRUE, we use different models for the
  # placebo arm in each trial. Otherwise, we use the same model for the placebo
  # arm in each trial.
  if (different_placebo_model) {
    # The first `n_trials` rows of theta and gamma contain the model parameters
    # for the control groups. The last `n_trials` rows contain the model
    # parameters for the active treatment groups.
    for (t in 1:n_trials) {
      for (a in c(0, t)) {
        idx <- which(sample_df$A == a & sample_df$trial == t)
        
        param_row <- ifelse(a == 0, t, t + n_trials)
        
        theta_a <- theta[param_row, ]
        gamma_a <- gamma[param_row, ]
        
        mu_Y <- inv_logit(model_matrix_Y[idx, ] %*% theta_a)
        mu_S <- model_matrix_S[idx, ] %*% gamma_a
        
        
        # If `bridge_error` is TRUE, we sample the residual for the linear model
        # for S from a bridge distribution instead of a normal distribution, and
        # we sample Y from a distribution conditional on S.
        if (bridge_error) {
          # Define phi parameter that gives unit variance.
          phi = 1 / sqrt(1 + 3 / pi^2)
          residuals_S <- ifelse(
            a == 0,
            0.28 * bridgedist::rbridge(length(idx), scale = phi),
            0.64 * bridgedist::rbridge(length(idx), scale = phi)
          )
          sample_df$S[idx] <- mu_S + residuals_S
          
          mu_Y = mu_Y - residuals_S
          sample_df$Y[idx] <- rbinom(length(idx), 1, mu_Y)
        } else {
          sample_df$S[idx] <- rnorm(length(idx), mu_S, ifelse(a == 0, 0.28, 0.64))
          
          sample_df$Y[idx] <- rbinom(length(idx), 1, mu_Y)
        }
      }
    }
  } else {
    # The first row of theta and gamma contain the model parameters for the
    # control groups. The other `n_trials` rows contain the model parameters for
    # the active treatment groups.
    for (a in 0:n_trials) {
      idx <- which(sample_df$A == a)
      
      theta_a <- theta[a + 1, ]
      gamma_a <- gamma[a + 1, ]
      
      mu_Y <- inv_logit(model_matrix_Y[idx, ] %*% theta_a)
      mu_S <- model_matrix_S[idx, ] %*% gamma_a
      
      sample_df$Y[idx] <- rbinom(length(idx), 1, mu_Y)
      sample_df$S[idx] <- rnorm(length(idx), mu_S, ifelse(a == 0, 0.28, 0.64))
    }
  }
 
  
  # ==================================================
  # CASE–COHORT DESIGN
  # ==================================================
  if (CC_sampling) {
    # For subjects with Y = 1, the sampling probability is 0.75. For subjects
    # with Y = 0, the sampling probability depends on X1: 0.10 if X1 = 0 and
    # 0.15 if X1 = 1.
    sample_df$Delta <- ifelse(sample_df$Y == 1,
                              rbinom(N, 1, 0.75),
                              ifelse(sample_df$X1 == 0, rbinom(N, 1, 0.10), rbinom(N, 1, 0.15)))
  }
  
  # ----------------------------
  # Final selection
  # ----------------------------
  xvy  <- unlist(strsplit(as.character(formula_Y)[3], " "))
  cols <- c("trial", "A", "Y", "S", xvy, "Z", "Delta", "vax")
  
  dplyr::select(sample_df, any_of(cols))
}

sample_covariates_trials <- function(zeta, N, target_trial) {
  # If `target_trial` is NULL, sample covariates and trial assignment
  # according to the "default" DGP.
  if (is.null(target_trial)) {
    df_list <- replicate(N, sample_row_covariates_trials(zeta), simplify = FALSE)
  } else {
    # If `target_trial` contains an integer, sample covariates from the target
    # trial covariate distribution and sample the trial indicator from a uniform
    # distribution on the trial indicators.
    df_list <- replicate(N, sample_one_covariates_trials_target(zeta, target_trial), simplify = FALSE)
  }
  df <- do.call(rbind, df_list) %>%
    as_tibble()
  colnames(df) <- c("X1", "X2", "trial")
  return(df)
}

# Function to sample covariates and trial assignment for one subject. 
sample_row_covariates_trials <- function(zeta) {
  # Sample covariates for one subject.
  X1 <- rbinom(1, 1, 0.3)
  X2 <- pmin(pmax(rnorm(1, 55, 22.3), 18), 85)
  X <- c(1, X1, X2)

  # Sample the trial indicator.
  denom <- 1 +
    exp(X %*% zeta[1, ]) +
    exp(X %*% zeta[2, ]) +
    exp(X %*% zeta[3, ]) +
    exp(X %*% zeta[4, ])
  
  probs <- c(
    1 / denom,
    exp(X %*% zeta[4, ]) / denom,
    exp(X %*% zeta[3, ]) / denom,
    exp(X %*% zeta[2, ]) / denom,
    exp(X %*% zeta[1, ]) / denom
  )
  
  trial <- which.max(rmultinom(1, 1, probs))
  
  return(c(X1, X2, trial))
  # return(data.frame(X1 = X1, X2 = X2, trial = trial))
}

# Sample covariates from the target trial covariate distribution and sample the
# trial indicator from a uniform distribution on the trial indicators.
sample_one_covariates_trials_target <- function(zeta, target_trial) {
  sample_vec <- sample_row_covariates_trials(zeta)
  while (sample_vec[3] != target_trial) {
    # We can sample covariate from the target trial covariate distribution by
    # sampling from the default DGP until the trial indicator is equal to the
    # target trial. 
    sample_vec <- sample_row_covariates_trials(zeta)
  }
  # Number of trial.
  n_trials <- nrow(zeta) + 1
  # Sample trial indicator from a uniform distribution on the trial indicators.
  sample_vec[3] <- sample(1:n_trials, 1)
  
  return(sample_vec)
}
