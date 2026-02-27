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
    ZOPT = NULL,
    target_trial = NULL
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
  
  df$Y <- NA_real_
  df$S <- NA_real_
  
  # ----------------------------
  # Generate potential outcomes
  # Treatments are 0,1,...,n_trials
  # ----------------------------
  for (a in 0:n_trials) {
    
    idx <- which(A == a)
    if (length(idx) == 0) next
    
    theta_a <- theta[a + 1, ]
    gamma_a <- gamma[a + 1, ]
    
    mu_Y <- inv_logit(X[idx, ] %*% theta_a)
    mu_S <- X[idx, ] %*% gamma_a
    
    df$Y[idx] <- rbinom(length(idx), 1, mu_Y)
    df$S[idx] <- rnorm(length(idx), mu_S, ifelse(a == 0, 0.28, 0.64))
  }
  
  # ==================================================
  # CASE–COHORT DESIGN
  # ==================================================
  if (isTRUE(CC_sampling)) {
    
    if (ZOPT == 1) {
      Z <- interaction(X1, X2 >= 65, drop = TRUE)
    } else if (ZOPT == 2) {
      Z <- factor(X1)
    } else {
      stop("ZOPT must be 1 or 2 when CC_sampling=TRUE")
    }
    
    df$Z <- Z
    
    # Trial-specific sampling probability
    D1pr <- rnorm(n_trials, mean = 0.1, sd = 0.03)
    df$D1pr <- D1pr[df$trial]
    
    df$D1 <- 0
    df$D2 <- 0
    
    for (i in seq_len(n_trials)) {
      
      idx_trial <- which(df$trial == i)
      
      # Stage 1 sample
      n1 <- ceiling(length(idx_trial) * D1pr[i])
      s1 <- sample(idx_trial, n1)
      df$D1[s1] <- 1
      
      # Stage 2 sample (cases)
      idx_cases <- idx_trial[df$Y[idx_trial] == 1]
      if (length(idx_cases) > 0) {
        n2 <- ceiling(length(idx_cases) * 0.75)
        s2 <- sample(idx_cases, n2)
        df$D2[s2] <- 1
      }
    }
    
    df$Delta <- as.integer(df$D1 == 1 | df$D2 == 1)
  }
  
  # ----------------------------
  # Final selection
  # ----------------------------
  xvy  <- unlist(strsplit(as.character(formula_Y)[3], " "))
  cols <- c("trial", "A", "Y", "S", xvy, "Z", "Delta", "vax")
  
  dplyr::select(df, any_of(cols))
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
  df <- do.call(rbind, df_list)
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
  
  return(data.frame(X1 = X1, X2 = X2, trial = trial))
}

# Sample covariates from the target trial covariate distribution and sample the
# trial indicator from a uniform distribution on the trial indicators.
sample_one_covariates_trials_target <- function(zeta, target_trial) {
  while (TRUE) {
    # We can sample covariate from the target trial covariate distribution by
    # sampling from the default DGP until the trial indicator is equal to the
    # target trial. 
    sample_df <- sample_one_covariates_trials(zeta)
    if (sample_df$trial == target_trial) {
      break
    }
  }
  # Number of trial.
  n_trials <- length(zeta) + 1
  # Sample trial indicator from a uniform distribution on the trial indicators.
  sample_df$trial <- sample(1:n_trials, 1)
  
  return(sample_df)
}
