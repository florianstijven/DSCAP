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
generate_simulated_data <- function(n_trials,n_t,Ymod,Smod,gamma.set,theta.set,CC=FALSE,ZOPT) {
  
  xvy <- unlist(strsplit(as.character(Ymod)[3],split=" "))
  
  # generate covariates
  X1 <- rbinom(prob = 0.3, n = n_trials * n_t, size = 1)
  X2 <- rnorm(n_trials * n_t, mean = 55, sd = 22.3)
  
  # truncate age (X2) into [18, 85]
  X2[X2 < 18] <- 18
  X2[X2 > 85] <- 85
  
  # sample T from T | X 
  X <- cbind(1, X1, X2)
  
  # parameters - higher intercept corresponds to higher n
  # zeta <- c(-0.17, 1, 0) #p5
  # eta <- c(1.2, 1.5, -0.025) #p4
  # lambda <- c(2.35, 2, -0.05) #p3
  # xi <- c(3.02, 3, -.075) #p2
  
  multinomial_probs_denom <- (1 + exp(X %*% zeta) + exp(X %*% eta) + 
                                exp(X %*% lambda) + exp(X %*% xi))
  
  p5 <- exp(X %*% zeta) / multinomial_probs_denom 
  p4 <- exp(X %*% eta) / multinomial_probs_denom 
  p3 <- exp(X %*% lambda) / multinomial_probs_denom 
  p2 <- exp(X %*% xi) / multinomial_probs_denom 
  p1 <- 1 - p2 - p3 - p4 - p5
  
  probs <- matrix(c(p1, p2, p3, p4, p5), ncol = 5)
  trial <- rMultinom(probs = probs, m = 1)
  
  df <- data.frame(trial, X1, X2) |> 
    rowwise() |> 
    mutate(vax = rbinom(1,1, .5), 
           A = ifelse(vax == 1, trial, 0)) 
  
  df$ID <- 1:nrow(df)
  df$Y <- NA
  df$S <- NA
  
  for(a in 1:(n_trials+1)){
    Y_a <- paste("Y_", a-1, sep="")
    S_a <- paste("S_", a-1, sep="")
    #NOTE: Theta and Gamma lists are defined in the simulation settings
    theta_a <- matrix(theta.set[[a]], ncol=1)
    gamma_a <- matrix(gamma.set[[a]], ncol=1)
    
    # inv_logit(X %*% theta_a) defines the TRUE probability of response if the subject got treatment A==a
    # rbinom introduces random error (draw from bernoulli(p))
    # inv_logit(x):= exp(x) / (1 + exp(x))
    df[[Y_a]] <- rbinom(nrow(df), 1, 
                        prob = (exp(X %*% theta_a)/(1 + exp(X %*% theta_a))))
    
    # true mean of Surrogate Outcome is given by X %*% gamma(a). True variance defined here.
    # Draw from normal(mean,sd) introduces random error to simulation
    ## Case Cohort design assigns probability of OBSERVING 
    if(a == 1){
      df[[S_a]] <- rnorm(nrow(df), mean = (X %*% gamma_a), sd = 0.28)
    } else {
      df[[S_a]] <- rnorm(nrow(df), mean = (X %*% gamma_a), sd = 0.64)
    }
    
    # Create a label for the potential outcome for clinical outcome Y for
    # treatment a-1 (i.e. 0,1,2,...,ntrial)
    Y_a <- paste("Y_", a-1, sep="")
    # Set the Y value for those with actual treatment A=a to the potential
    # outcome Y_(a-1)
    df[df$A==(a-1),]$Y <- unlist(df[df$A==(a-1),Y_a])
    # Create a label for the potential outcome for surrogate outcome S for
    # treatment a-1 (i.e. 0,1,2,...,ntrial)
    S_a <- paste("S_", a-1, sep="")
    # Set the S value for those with actual treatment A=a to the potential
    # outcome S_(a-1)
    df[df$A==(a-1),]$S <- unlist(df[df$A==(a-1),S_a])
  }
  
  if (CC==T) {
    # CASE-COHORT UPDATE
    if (ZOPT==1) {
      # Define strata Z (1 through 6) based on A (treatment), X1 (0/1), and X2 (<65,>=65)
      XA <- as.data.frame(cbind(X,A))
      XA$X2. <- (XA$X2 >= 65)
      Z <- factor(with(XA, paste(X1, X2., sep = "_")))
      levels(Z) <- c(1:4)
    }
    else if (ZOPT==2) {
      XA <- as.data.frame(cbind(X,A))
      Z <- factor(XA$X1)
      levels(Z) <- c(1:2)
    }
    
    ## Define delta: 
    ## stage 1: (D1) equal proportion sampled from each category of Z within each trial
    ## stage 2: (D2) 75% of cases sampled within each trial
    ## Delta = 1 if D1==1 OR D2==1
    df$Z <- Z
    df <- df %>% mutate(ID = cur_group_id(),
                        TZ = paste0(trial,"_",Z),
                        TY = paste0(trial,"_",Y),
                        TZY = paste0(trial,"_",Z,"_",Y))
    
    #Now define sampling fractions
    D1pr <- rnorm(n=n_trials,mean=0.1,sd=0.03)
    
    df$D1pr <- ifelse(df$trial==1,D1pr[1],
                      ifelse(df$trial==2,D1pr[2],
                             ifelse(df$trial==3,D1pr[3],
                                    ifelse(df$trial==4,D1pr[4],
                                           ifelse(df$trial==5,D1pr[5],NA)))))
    #Select a sample from each trial for D1 and D2
    df$D1 = 0
    df$D2 = 0
    for (i in 1:n_trials) {
      temp.df <- df %>% filter(trial==i)
      df$D1pr[df$trial==i] <- D1pr[i]
      temp <- as.vector(createDataPartition(temp.df$TY, p = D1pr[i], list = FALSE))
      ids <- temp.df$ID[temp]
      df$D1[df$ID %in% ids] <- 1
      
      temp.df2 <- temp.df %>% filter(Y==1)
      temp2 <- as.vector(createDataPartition(temp.df2$TY, p = 0.75, list = FALSE))
      id2 <- temp.df2$ID[temp2]
      df$D2[df$ID %in% id2] <- 1
    }
    df$Delta <- 0
    df$Delta[df$D1==1 | df$D2==1] <- 1
  }
  
  cols <- c("trial", "A", "R", "Y", "S", xvy, "Z", "Delta", "vax")
  dfa <- df %>% dplyr::select(any_of(cols))
  
  return(dfa)
}
