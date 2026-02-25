# Preliminaries -----------------------------------------------------------
t1 <- Sys.time()

# Load all R packages
library(tidyverse)
# A custom version of the geex package is used to avoid issues with large
# objects.
# install.packages("geex_1.1.1.tar.gz", repos = NULL)
# library(geex, lib.loc = "/home/srosin/sim_code/pkgs/")
library(nnet) # multinom() function for multinomial logistic regression
library(Hmisc)

# Set up parallel computing
if (parallelly::supportsMulticore()) {
  plan("multicore", workers = parallel::detectCores() - 1)
} else {
  plan(multisession, workers = parallel::detectCores() - 1)
}
# Extract arguments for analysis.
args = commandArgs(trailingOnly = TRUE)

setting <- as.numeric(args[1])
# Number of trials in each simulation.
n_trials <- as.numeric(args[2])
# Number of subjects in each trial.
n_t <- as.numeric(args[3])
# Number of bootstrap replications for inference.
n_boot <- as.numeric(args[4])
# Number of replications for permutation test for conditional independence.r
n_perm <- as.numeric(args[5])
vars <- as.character(args[6])
# Case-cohort sampling indicator. If TRUE, the data are generated according to a
# case-cohort sampling design. If FALSE, the data are generated according to a
# full-cohort design.
CC_sampling <- as.logical(args[7])
# Indicator for whether the weights are (re-)estimated for the bootstrap.
estimate_weights <- as.logical(args[8])

# Define formulas for the outcome models for (i) clinical outcome Y and (ii)
# surrogate outcome S. The same linear predictors are used for both models. In
# the simulations, the corresponding models are automatically stratified by
# treatment.
formula_Y = as.formula(paste0("Y ~ ",vars))
formula_S = as.formula(paste0("S ~ ",vars))

# Define formula for trial participation model given covariates X (but not
# treatment).
formula_T = as.formula(paste0("trial ~ ",vars))


df <- generate_simulated_data(n_trials = n_trials, 
                              n_t = n_t,
                              Ymod,
                              Smod, 
                              gamma.set, 
                              theta.set,
                              CC = CC_sampling, 
                              ZOPT = NA)


# Helper Functions --------------------------------------------------------

# Source helper functions.
source("R/helper-functions/DSCAP-estimators.R")
source("R/helper-functions/treatment-effect-estimators.R")

## Simulation Function --------------------------------------------------

# Function to analyze one data set.
analyze <- function() {
  
}

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
  dplyr::summarise(
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

t2 <- Sys.time()
print(t2 - t1)