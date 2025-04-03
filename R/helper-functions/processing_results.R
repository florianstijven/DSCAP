# Helper function to extract trial names from column names. When saving results
# in .csv format, some illegal characters are modified in the column names.
trial_name_convert = function(modified_name) {
  # The last dot should be removed with the new results.
  original_name = switch(
    modified_name,
    J.J..Brazil. = "J&J (Brazil)",
    AstraZeneca. = "AstraZeneca",
    Novavax. = "Novavax",
    Moderna. = "Moderna",
    J.J..Colombia.. = "J&J (Colombia)",
    J.J..S..America.. = "J&J (S. America)",
    J.J..S..Africa.. = "J&J (S. Africa)",
    J.J..USA.. = "J&J (USA)",
    AstraZeneca = "AstraZeneca",
    Novavax = "Novavax",
    Moderna = "Moderna",
    J.J..Colombia. = "J&J (Colombia)",
    J.J..S..America. = "J&J (S. America)",
    J.J..S..Africa. = "J&J (S. Africa)",
    J.J..USA. = "J&J (USA)",
    modified_name
  )
  return(original_name)
}

# Helper function to read in all results for a single analysis. 
read_results = function(target, modn, surr_type, estimate_weights, nX, truncation) {
  # Position of parameters in the estimated covariance matrix ----
  n_trials = ifelse(modn == 1, 8, 5)
  # First position of the treatment effect parameters.
  trt_effect_start = (2 * n_trials * nX - nX) 
  # Position of parameters corresponding to the standardized mean clinical
  # outcomes.
  clin_stand_parm_pos = (trt_effect_start + 1):(trt_effect_start + 1 + n_trials)
  
  # Position of parameters corresponding to the standardized mean surrogate
  # outcomes.
  surr_stand_parm_pos = (trt_effect_start + (1 + n_trials) + 1):(trt_effect_start  + (1 + n_trials) * 2 - 1)
  
  # Position of parameters corresponding to tau 1 to `n_trials`.
  tau_1_to_8_pos = (trt_effect_start + (1 + n_trials) * 2):(trt_effect_start + (1 + n_trials) * 3 - 2)
  
  # Position of correlation parameters (Pearson correlation and beta; in this order).
  corr_pos = (trt_effect_start + (1 + n_trials) * 3 - 1):(trt_effect_start + (1 + n_trials) * 3)
  
  
  target_trial = switch(
    target,
    jjsa = "J&J (S. America)",
    AZ = "AstraZeneca",
    jjbr = "J&J (Brazil)",
    jjusa  = "J&J (USA)",
    jjcol = "J&J (Colombia)",
    nov = "Novavax"
  )

  # Import results ----
  if (estimate_weights) weights_chr = "estwts" else weights_chr = "kwnwts"
  if (truncation) population = "truncated" else population = "full"
  analysis_infile = paste0(surr_type, "_", target, "_", population, "_M", modn, "_", weights_chr, ".csv")
  
  # Import data that contains the results of the analysis with parameters `target`,
  # `modn`, and `v`.
  bootstrap_df <- read.csv(paste0(raw_results_dir, "bootstrap_", analysis_infile))
  results_vcov_df <- read.csv(paste0(raw_results_dir, "vcov_", analysis_infile))
  lrt_results_df <- read.csv(paste0(raw_results_dir, "lrt_", analysis_infile))
  
  trt_effects_df = read.csv(paste0(
    raw_results_dir,
    "trt_effects_",
    analysis_infile
  )) %>%
    select(
      trial,
      VE_est,
      mean_diff_S_est,
      VE_CI_lower,
      VE_CI_upper,
      mean_diff_S_CI_lower,
      mean_diff_S_CI_upper, 
      type
    )
  
  cor_est_df = read.csv(paste0(
    raw_results_dir, "cor_est_", analysis_infile
  ))
  # Data wrangling to simplify plotting ----
  
  # # Get the order in which the trials were saved.
  # trials_chr_ordered = standardized_trt_effects_df %>%
  #   pull(trial)
  
  # Build data set with confidence intervals for all estimates.
  pm_positions = c(clin_stand_parm_pos,
                   surr_stand_parm_pos,
                   tau_1_to_8_pos,
                   corr_pos)
  
  sandwich_se_df = tibble(
    estimand = c(
      rep("mean_Y", n_trials + 1),
      rep("mean_diff_S", n_trials),
      rep("VE", n_trials),
      "rho_p",
      "beta"
    ),
    type = rep("standardized", 3 * n_trials + 3),
    trial_name = colnames(results_vcov_df)[pm_positions] %>% stringr::str_remove(pattern = ".*\\.\\.\\.") %>%
      stringr::str_remove(pattern = "1"), # This last operation should be removed with the newest version of the results.
    se = results_vcov_df[pm_positions, pm_positions] %>%
      as.matrix() %>%
      diag() %>%
      sqrt()
  ) %>%
    mutate(trial = sapply(trial_name, trial_name_convert))
  
  # Add upper and lower limit of confidence intervals, if available, in two
  # separate columns.
  trt_effects_df_long = trt_effects_df %>%
    pivot_longer(cols = !c(trial, type)) %>%
    mutate(
      estimand = ifelse(stringr::str_detect(name, "VE"), "VE", "mean_diff_S"),
      value_type = ifelse(
        stringr::str_detect(name, "CI_lower"),
        "CI_lower",
        ifelse(stringr::str_detect(name, "CI_upper"), "CI_upper", "estimate")
      )
    ) %>% select(-name) %>%
    pivot_wider(names_from = "value_type", values_from = "value")
  
  all_results_df = bind_rows(trt_effects_df_long,
                             cor_est_df %>%
                               pivot_longer(
                                 cols = !c("type"),
                                 values_to = "estimate", 
                                 names_to = "estimand"
                               ) ) 
  
  
  # Attach sandwich SE estimates.
  all_results_df = all_results_df %>%
    left_join(
      sandwich_se_df %>%
        select(-trial_name) %>%
        mutate(trial = ifelse(estimand %in% c("rho_p", "beta"), NA, trial),
               estimand = ifelse(estimand == "rho_p", "cor_p", estimand))
    )
  
  # Compute CIs based on SE estimates.
  all_results_df = bind_rows(
    all_results_df %>%
      filter(type == "naive"),
    all_results_df %>%
      filter(estimand == "VE", type == "standardized") %>%
      mutate(
        log_RR_se = sqrt(((1 - estimate) ** -2) * se ** 2),
        # Compute CI limits on the log(1 - VE) scale.
        CI_lower = log(1 - estimate) + 1.96 * log_RR_se,
        CI_upper = log(1 - estimate) - 1.96 * log_RR_se
      ) %>%
      # Transform CIs back to the VE scale.
      mutate(CI_lower = 1 - exp(CI_lower), CI_upper = 1 - exp(CI_upper)),
    all_results_df %>%
      filter(estimand %in% c("mean_diff_S", "beta"), type == "standardized") %>%
      mutate(
        CI_lower = estimate - 1.96 * se,
        CI_upper = estimate + 1.96 * se
      ),
    all_results_df %>%
      filter(estimand == "cor_p", type == "standardized") %>%
      mutate(
        # Standard error on Fisher's Z scale. 
        fish_z_se = (1 / (1 - estimate ** 2)) * se,
        # Compute CI limits.e.
        CI_lower = tanh(atanh(estimate) - 1.96 * fish_z_se),
        CI_upper = tanh(atanh(estimate) + 1.96 * fish_z_se)
      ),
    all_results_df %>%
      filter(estimand == "cor_s", type == "standardized")
  )

  
  # Add CIs and SEs based on the bootstrap.
  all_results_df = all_results_df %>%
    full_join(
      bootstrap_df %>%
        select(-mean) %>%
        rename(
          CI_lower_bs = "CI_lower",
          CI_upper_bs = "CI_upper",
          se_bs = "se"
        )
    ) %>%
    mutate(estimate_weights = estimate_weights)
  
  return(all_results_df)
}
