
# setup and data generation -----------------------------------------------
setwd("H:/sim_code/")

## setup -------------------------------------------------------------------
t1 <- Sys.time()

library(tidyverse)
# A custom version of the geex package is used to avoid issues with large
# objects.
# install.packages("geex_1.1.1.tar.gz", repos = NULL)
library(geex)
library(nnet) # multinom() function for multinomial logistic regression
library(Hmisc)
library(caret)

# Extract arguments for analysis. 
# args = commandArgs(trailingOnly=TRUE)
# antibody_type <- args[1]
# target <- args[2]
# vars = args[3]
# mnum <- as.integer(args[4])
# n_boot <- as.integer(args[5])
# n_perm <- as.integer(args[6])
# estimate_weights = as.integer(args[7])
# truncation = as.integer(args[8])
# data_location = args[9]

# the below variables will come from the args
vars <- "X1 + X2"
n_trials <- 5
n_boot <- 3
n_perm <- 10
iteration <- 100
CC_sampling <- FALSE
estimate_weights <- FALSE
setting <- 2
n_t <- 1500

set.seed(iteration)

# load helper functions
source("helper_functions/parameter_values.R") # load parameter settings
source("helper_functions/generate_simulated_data.R")
#source("helper_functions/generate_simulated_data_covariate_dists_equal.R")
source("helper_functions/runDSCAP.R")  
source("helper_functions/estFUN.R")

# Define formulas for the outcome models for (i) clinical outcome Y and (ii)
# surrogate outcome S. The same linear predictors are used for both models. 
Ymod = formula(Y ~ X1 + X2)
formula_Y = formula(Y ~ X1 + X2)
Smod = formula(S ~ X1 + X2)
formula_S = formula(S ~ X1 + X2)
formula_T = formula(trial ~ X1 + X2)

# Define formula for trial participation model
Tmod = as.formula(paste0("trial ~ ",vars))

n_trials = 5
n_t = 5e3
trial_var = "trial"
treatment_var = "A"
CC_weight_var = NULL
target_trial = 1

df <- generate_simulated_data(n_trials = n_trials, 
                              n_t = n_t,
                              Ymod,
                              Smod, 
                              gamma.set = gamma.set, 
                              theta.set = theta.set,
                              CC = FALSE, 
                              ZOPT = 1)

# Add variable that defines the target population. 0: for target populations; 1
# otherwise. trial #1 is defined to be the target trial
df$R <- ifelse(df$trial == 1, 0, 1)

# output directory to store results. note that "results" subdir must exist
wd <- getwd()
outdir = paste0(
  wd,
  "/results/setting",
  setting,
  "_ntrials",
  n_trials#,
  # "/"
)
if (!dir.exists(outdir)) { dir.create(outdir) }

# data analysis -----------------------------------------------------------

## point estimation --------------------------------------------------------------

# Estimate all models.
RESULT <- RunDSCAP(
  data = df,
  formula_S = Smod,
  formula_Y = Ymod,
  formula_T = Tmod,
  target_trial = 1,
  # sim = T,
  estimate_weights = FALSE,
  CC_sampling = CC_sampling,
  alpha_level = 0.05
)

if(CC_sampling){
  # If some some strata have an estimated weight of Inf (i.e., zero probability of
  # being sampled in the case-cohort sampling, they are excluded from the further
  # analyses).
  df = df %>%
    filter(!(CC_stratum %in% RESULT$excluded_CC_strata))
  # Include estimated weights to data.
  df = df %>%
    left_join(
      RESULT$weights_df
    )
}

# write results
outfile = paste0(
  outdir,
  c("/trt_effects_", "/cor_est_"),
  iteration,
  ".csv"
)

RESULT$naive_trt_effects_df$trial <- as.factor(RESULT$naive_trt_effects_df$trial)

write.csv(
  bind_rows(
    RESULT$standardized_trt_effects_df %>%
      mutate(type = "standardized"),
    RESULT$ipw_trt_effects_df %>%
      mutate(type = "ipw"),
    RESULT$naive_trt_effects_df %>%
      mutate(type = "naive")
  ),
  outfile[1],
  row.names = FALSE
)

write.csv(
  bind_rows(
    RESULT$cor_standardized_df %>%
      mutate(type = "standardized"),
    RESULT$cor_ipw_df %>% 
      mutate(type = "ipw"),
    RESULT$cor_naive_df %>%
      mutate(type = "naive")
  ),
  outfile[2],
  row.names = FALSE
)

## conditional independence tests ------------------------------------------


# Fit null model (which does not contain trial as covariate).
glm_null <- glm(as.formula(paste0("Y ~ ", vars)), binomial, df %>%
                  filter(A == 0))
# Fit alternative model (which contains the interaction between trial and all
# covariates which were  already in the null model).
glm_w_trial <-update(glm_null, ~ . * trial)
# Perform likelihood ratio test. 
glm_lr_test = lmtest::lrtest(glm_null, glm_w_trial)

# Extract LRT statistic and corresponding p-value.
lrt_test_statistic <- glm_lr_test[2,4]
lrt_p_value <- glm_lr_test[2,5]

# Initialize data set in which the values of `trial` will be permuted. 
df_placebo_perm = df %>%
  filter(A == 0)
# Extract trial variable which will be permuted.
trial_vec = df_placebo_perm$trial

lrt_test_statistic_permutations = rep(NA, n_perm)
for(i in 1:n_perm){
  df_placebo_perm$trial <- sample(trial_vec, replace = FALSE)
  
  glm_w_trial <-update(glm_null, ~ . * trial, data = df_placebo_perm)
  lrt_test_statistic_permutations[i] <- lmtest::lrtest(glm_null, glm_w_trial)[2,4] #extract LRT stat
}

lrt_permuted_p_value <- mean(lrt_test_statistic_permutations > lrt_test_statistic)
lrt_out <- c("p-value" = lrt_p_value, "p-value (permutations)" = lrt_permuted_p_value, "LRT statistic" = lrt_test_statistic)

outfile_lrt = paste0(
  outdir,
  "/lrt_",
  iteration,
  ".csv"
)
write.csv(lrt_out, outfile_lrt, row.names=FALSE)


## bootstrap variance estimation --------------------------------------------------------------

# Initialize dataframe in which the bootstrap replicates will be saved. Note the
# the results of a single bootstrap replication are saved over many rows. This
# facilities further processing.
bootstrap_replicates_df = data.frame(
  estimand = character(0),
  type = character(0),
  trial = character(0),
  estimate = numeric(0)
)

# Perform bootstrap by resampling within each trial.
for (i in 1:n_boot) {
  # Resample within each trial
  booti <- df %>%
    group_by(trial) %>%
    slice_sample(prop = 1.0, replace = TRUE)
  
  try(expr =  {
    temp <- RunDSCAP(
      data = booti,
      formula_S = Smod,
      formula_Y = Ymod,
      formula_T = Tmod,
      target_trial = 1,
      # sim = T,
      estimate_weights = FALSE,
      CC_sampling = CC_sampling,
      alpha_level = 0.05
    )
    bootstrap_replicates_df = bootstrap_replicates_df %>%
      bind_rows(
        tibble(
          estimand = c("cor_p", "cor_s", "beta"),
          type = "standardized",
          trial = NA,
          estimate = temp$cor_standardized_df[1, ] %>% as.numeric()
        ),
        tibble(
          estimand = c("cor_p", "cor_s", "beta"),
          type = "ipw",
          trial = NA,
          estimate = temp$cor_standardized_df[1, ] %>% as.numeric()
        ),
        tibble(
          estimand = c("cor_p", "cor_s", "beta"),
          type = "naive",
          trial = NA,
          estimate = temp$cor_naive_df[1, ] %>% as.numeric()
        ),
        tibble(
          estimand = "VE",
          type = "standardized",
          trial = temp$standardized_trt_effects_df %>% pull(trial),
          estimate = temp$standardized_trt_effects_df %>% pull(VE_est)
        ),
        tibble(
          estimand = "mean_diff_S",
          type = "standardized",
          trial = temp$standardized_trt_effects_df %>% pull(trial),
          estimate = temp$standardized_trt_effects_df %>% pull(mean_diff_S_est)
        ),
        tibble(
          estimand = "VE",
          type = "ipw",
          trial = temp$ipw_trt_effects_df %>% pull(trial),
          estimate = temp$ipw_trt_effects_df %>% pull(VE_est)
        ),
        tibble(
          estimand = "mean_diff_S",
          type = "ipw",
          trial = temp$ipw_trt_effects_df %>% pull(trial),
          estimate = temp$ipw_trt_effects_df %>% pull(mean_diff_S_est)
        ),
        tibble(
          estimand = "VE",
          type = "naive",
          trial = temp$naive_trt_effects_df %>% pull(trial) %>% as.factor(),
          estimate = temp$naive_trt_effects_df %>% pull(VE_est)
        ),
        tibble(
          estimand = "mean_diff_S",
          type = "naive",
          trial = temp$naive_trt_effects_df %>% pull(trial) %>% as.factor(),
          estimate = temp$naive_trt_effects_df %>% pull(mean_diff_S_est)
        )
      )
  }
  )
}

bootstrap_inferences_df = bootstrap_replicates_df %>%
  group_by(estimand, type, trial) %>%
  summarise(
    CI_lower = quantile(estimate, 0.025, na.rm = TRUE),
    CI_upper = quantile(estimate, 0.975, na.rm = TRUE),
    mean = mean(estimate, na.rm = TRUE),
    se = sd(estimate, na.rm = TRUE)
  )

outfile_bootstrap = paste0(
  outdir,
  "/bootstrap_",
  iteration,
  ".csv"
)

write.csv(bootstrap_inferences_df, outfile_bootstrap, row.names=FALSE)


## sandwich variance estimation -----------------------------------------------------

if(CC_sampling){
  
  # If some of the estimated weights are equal to 1, this will break the sandwich
  # estimator because the corresponding bread matrix will not be invertible. 
  problem_weights = RESULT$weights_df %>%
    filter(CC_stratum != "Placebo", weight == 1) %>%
    summarize(n() >= 1) %>%
    pull() %>% # Returns TRUE if there is an estimated weight of exactly 1.
    any()
  
  if (estimate_weights & problem_weights) {
    warning(
      "At least one weight stratum has an estimated weight equal to one. The original sandwich estimator fails because the bread matrix is not invertible. A modified sandwich estimator is obtained by treating the problematic weight stratum's weight as fixed."
    )
    # The problematic stratum is joined with the placebo stratum. This forces the
    # code to treat the corresponding (estimated) weight as fixed. 
    problematic_strata = RESULT$weights_df %>%
      filter(CC_stratum != "Placebo", weight == 1) %>%
      pull(CC_stratum) %>%
      unique()
    df = df %>%
      mutate(CC_stratum = ifelse(CC_stratum %in% problematic_strata, "Placebo", CC_stratum))
    
    # The models and weights are re-estimated, now with the modified weight
    # strata. Note that all results remain unchanged, except that there is now
    # one weight stratum fewer.
    RESULT <- RunDSCAP(
      data = df,
      formula_S = Smod,
      formula_Y = Ymod,
      target_trial = 1,
      sim = F, 
      estimate_weights = TRUE
    )
  }
}

# Extract estimated parameter vector that corresponds to the set of stacked
# estimating equations.
theta = extract_coefs(RESULT, estimate_weights, target_trial = 1)

if(CC_sampling){
  # Attach estimated weights to original data frame. These weights will be used if
  # the weights are treated as known for the sandwich variance estimator.
  df = df %>%
    select(-weight) %>% left_join(RESULT$weights_df, by = "CC_stratum")
}

# Compute sandwich estimate
if(CC_sampling){
  m_est = m_estimate(
    estFUN = estFUN_taudelta_CC,
    data = df,
    outer_args = list(
      models_tbl = RESULT$outcome_model_fits_df %>%
        dplyr::select(A_modified, outcome_model_fit_Y, outcome_model_fit_S),
      weights_df = RESULT$weights_df,
      target_trial = 1
    ),
    inner_args = list(
      estimate_weights = estimate_weights
    ),
    roots = theta,
    compute_roots = FALSE,
    deriv_control = setup_deriv_control(method = "simple"),
    compute_vcov = TRUE
  )
} else {
  m_est = m_estimate(
    estFUN = estFUN_taudelta,
    data = df,
    outer_args = list(
      models_tbl = RESULT$outcome_model_fits_df %>%
        dplyr::select(A_modified, outcome_model_fit_Y, outcome_model_fit_S),
      target_trial = 1
    ),
    roots = theta,
    compute_roots = FALSE,
    deriv_control = setup_deriv_control(method = "simple"),
    compute_vcov = TRUE
  )
}


vcov_m_est <- data.frame(vcov(m_est))
colnames(vcov_m_est) = names(theta)


# outdir2 = paste0(
#   outdir,
#   "/vcovs"
# )
# if (!dir.exists(outdir2)) { dir.create(outdir2) }

outfile_vcov = paste0(
  outdir,
  "/vcov_",
  iteration,
  ".csv"
)
write.csv(vcov_m_est, outfile_vcov, row.names = FALSE)

t2 <- Sys.time()
print(t2 - t1)

