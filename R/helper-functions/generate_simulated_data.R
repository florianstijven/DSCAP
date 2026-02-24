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
    Ymod,
    Smod,
    gamma.set,
    theta.set,
    CC = FALSE,
    ZOPT = NULL
) {
  
  N <- n_trials * n_t
  
  inv_logit <- function(x) plogis(x)
  
  # ----------------------------
  # Generate covariates
  # ----------------------------
  X1 <- rbinom(N, 1, 0.3)
  X2 <- pmin(pmax(rnorm(N, 55, 22.3), 18), 85)
  
  X <- cbind(1, X1, X2)
  
  # ----------------------------
  # Multinomial trial assignment
  # ----------------------------
  denom <- 1 +
    exp(X %*% zeta) +
    exp(X %*% eta) +
    exp(X %*% lambda) +
    exp(X %*% xi)
  
  probs <- cbind(
    1 / denom,
    exp(X %*% xi)     / denom,
    exp(X %*% lambda) / denom,
    exp(X %*% eta)    / denom,
    exp(X %*% zeta)   / denom
  )
  
  trial <- max.col(
    t(apply(probs, 1, function(p) rmultinom(1, 1, p))),
    ties.method = "first"
  )
  
  # ----------------------------
  # Treatment assignment
  # ----------------------------
  # Trial a now has treatments {0, a}
  vax <- rbinom(N, 1, 0.5)
  A   <- ifelse(vax == 1, trial, 0)
  
  df <- data.frame(
    ID    = seq_len(N),
    trial = trial,
    X1    = X1,
    X2    = X2,
    vax   = vax,
    A     = A
  )
  
  df$Y <- NA_real_
  df$S <- NA_real_
  
  # ----------------------------
  # Generate potential outcomes
  # Treatments are 0,1,...,n_trials
  # ----------------------------
  for (a in 0:n_trials) {
    
    idx <- which(A == a)
    if (length(idx) == 0) next
    
    theta_a <- matrix(theta.set[[a + 1]], ncol = 1)
    gamma_a <- matrix(gamma.set[[a + 1]], ncol = 1)
    
    mu_Y <- inv_logit(X[idx, ] %*% theta_a)
    mu_S <- X[idx, ] %*% gamma_a
    
    df$Y[idx] <- rbinom(length(idx), 1, mu_Y)
    df$S[idx] <- rnorm(length(idx), mu_S, ifelse(a == 0, 0.28, 0.64))
  }
  
  # ==================================================
  # CASE–COHORT DESIGN
  # ==================================================
  if (isTRUE(CC)) {
    
    if (ZOPT == 1) {
      Z <- interaction(X1, X2 >= 65, drop = TRUE)
    } else if (ZOPT == 2) {
      Z <- factor(X1)
    } else {
      stop("ZOPT must be 1 or 2 when CC=TRUE")
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
  xvy  <- unlist(strsplit(as.character(Ymod)[3], " "))
  cols <- c("trial", "A", "Y", "S", xvy, "Z", "Delta", "vax")
  
  dplyr::select(df, any_of(cols))
}
