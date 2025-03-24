RunDSCAP <- function(data,
                     formula_Y,
                     formula_S,
                     target_trial,
                     sim = F,
                     estimate_weights = FALSE)
{
  df = data
  # Determine which level of the trial factor corresponds to the target trial. 
  target_trial_level_n = which(levels(data$trial) == target_trial)
  

  n_trials <- length(levels(data$trial))
  
  # Re-estimate weights if required.
  if (estimate_weights) {
    df <- df %>% group_by(CC_stratum) %>% mutate(prop.delta = sum(Delta == 1) /
                                                   n(),
                                                 weight = 1 / prop.delta) %>%
      dplyr::select(-prop.delta) %>%
      # PLacebo patients should get a weight of one.
      mutate(weight = ifelse(A == 0, 1, weight))
  }
  
  # Data frame that contains the estimated weight for each category in
  # `CC_stratum`.
  weights_df = df %>%
    group_by(CC_stratum) %>%
    slice_head(n = 1) %>%
    select(CC_stratum, weight)
  
  problematic_strata = NA
  if (any(weights_df$weight == Inf)) {
    problematic_strata = weights_df %>%
      filter(weight == Inf) %>%
      pull(CC_stratum)
    warning(paste0(
      paste(
        c("Estimated weight(s) of Inf for strata", problematic_strata),
        collapse = " "
      ),
      ". Corresponding observations are ignored for this analysis."
    ))
    weights_df = weights_df %>%
      filter(!(CC_stratum %in% problematic_strata))
    df = df  %>%
      filter(!(CC_stratum %in% problematic_strata))
  }

  
  glm_fits_df = df %>%
    mutate(trial_modified = ifelse(A == 0, "Placebo", as.character(trial)),
           trial_modified = factor(trial_modified)) %>%
    filter(trial != target_trial) %>%
    group_by(trial_modified) %>%
    summarize(
      glm_fit_Y = glm(formula_Y, data = pick(everything()), family = binomial, x = FALSE, y = FALSE) %>% list(),
      glm_fit_S = glm(
        formula_S,
        data = pick(everything()),
        family = gaussian,
        subset = (Delta == 1),
        weights = weight, 
        x = FALSE,
        y = FALSE
      ) %>% list()
    ) %>%
    # Extract coefficients vectors.
    ungroup() %>%
    mutate(
      glm_coef_Y = purrr::map(glm_fit_Y, coef),
      glm_coef_S = purrr::map(glm_fit_S, coef)
    )

  # For each subject in the target trial, we predict the outcome using the model
  # estimated in any of the other trials.
  target_trial_df = df %>% 
    filter(trial == target_trial) 
  
  glm_fits_df = glm_fits_df %>%
    mutate(
      predicted_Y = purrr::map(
        .x = glm_fit_Y,
        .f = function(glm_fit) {
          predict(glm_fit, newdata = target_trial_df, type = "response")
        }
      ),
      predicted_S = purrr::map(
        .x = glm_fit_S,
        .f = function(glm_fit) {
          predict(glm_fit, newdata = target_trial_df, type = "response")
        }
      )
    )
  
  
  # Compute the standardized means as the means of the predicted outcomes, based
  # on any of the trial-specific models, of subjects from the target trial. 
  standardized_means_df = glm_fits_df %>%
    mutate(
      mean_Y = purrr::map_dbl(
        .x = predicted_Y,
        .f = mean
      ),
      mean_S = purrr::map_dbl(
        .x = predicted_S,
        .f = mean
      )
    ) %>%
    select(trial_modified, mean_Y, mean_S)

  # Compute the naive trial-specific means. 
  naive_means_df = df %>%
    mutate(S = ifelse(is.na(S), 0, S)) %>%
    group_by(trial, vax) %>%
    summarise(
      mean_Y = mean(Y),
      n = n(),
      # n_obs = sum(!is.na(S)),
      # probs_obs = mean(!is.na(S)),
      # mean_weight = mean(ifelse(!is.na(S), weight, NA), na.rm = TRUE),
      mean_S = mean(weight * S),
      var_S = var(weight * S) 
    )
  
  # Compute treatment effects from the estimated means. For the naive
  # trial-specific estimates, we also compute confidence intervals. 
  naive_trt_effects_df = naive_means_df %>%
    pivot_wider(
      names_from = c("vax"),
      values_from = c("mean_Y", "n", "mean_S", "var_S")
    ) %>%
    mutate(
      VE_est = 1 - mean_Y_1 / mean_Y_0,
      log_RR_est = log(1 - VE_est),
      log_RR_se = sqrt(((1 - mean_Y_1) / (mean_Y_1 * n_1)) + ((1 - mean_Y_0) / (mean_Y_0 * n_0))),
      mean_diff_S_est = mean_S_1 - mean_S_0,
      mean_diff_S_se = sqrt((var_S_1 / n_1)  + (var_S_0 / n_0)),
    ) %>%
    mutate(
      VE_CI_lower = 1 - exp(log_RR_est + 1.96 * log_RR_se),
      VE_CI_upper = 1 - exp(log_RR_est - 1.96 * log_RR_se),
      mean_diff_S_CI_lower = mean_diff_S_est - 1.96 * mean_diff_S_se,
      mean_diff_S_CI_upper = mean_diff_S_est + 1.96 * mean_diff_S_se
    )
  
  standardized_trt_effects_df = standardized_means_df %>%
    bind_rows(
      naive_means_df %>%
        filter(trial == target_trial, vax == 1) %>%
        rename(trial_modified = trial) %>%
        select(trial_modified, mean_Y, mean_S)
    ) %>%
    filter(trial_modified != "Placebo") %>%
    rename(mean_Y_1 = mean_Y, mean_S_1 = mean_S) %>%
    cross_join(
      standardized_means_df %>%
        filter(trial_modified == "Placebo") %>%
        rename(mean_Y_0 = mean_Y, mean_S_0 = mean_S) %>%
        select(-trial_modified)
    ) %>%
    mutate(VE_est = 1 - mean_Y_1 / mean_Y_0,
           mean_diff_S_est = mean_S_1 - mean_S_0,) %>%
    rename(trial = "trial_modified")
  
  # Compute the measures of surrogacy based on the naive trial-specific and
  # standardized treatment effect estimates.
  cor_naive_df = naive_trt_effects_df %>%
    group_by() %>%
    summarize(
      cor_p = cor(VE_est, mean_diff_S_est, method = "pearson"),
      cor_s = cor(VE_est, mean_diff_S_est, method = "spearman"),
      beta = coef(lm(VE_est ~ mean_diff_S_est))[2]
    )

  cor_standardized_df = standardized_trt_effects_df %>%
    group_by() %>%
    summarize(
      cor_p = cor(VE_est, mean_diff_S_est, method = "pearson"),
      cor_s = cor(VE_est, mean_diff_S_est, method = "spearman"),
      beta = coef(lm(VE_est ~ mean_diff_S_est))[2]
    )
  
  
  # Combine everything into a list and return this list. 
  obj = 
    list(
      glm_fits_df = glm_fits_df %>%
        select(-predicted_Y, -predicted_S),
      standardized_means_df = standardized_means_df,
      standardized_trt_effects_df = standardized_trt_effects_df,
      naive_means_df = naive_means_df,
      naive_trt_effects_df = naive_trt_effects_df,
      cor_standardized_df = cor_standardized_df,
      cor_naive_df = cor_naive_df,
      weights_df = weights_df,
      target_trial = target_trial,
      excluded_CC_strata = problematic_strata
    )
  gc()
  
  return(obj)
}

extract_coefs = function(obj, estimate_weights) {
  # The ordering of the trial-specific estimates is determined by the ordering
  # of the rows in `obj$glm_fits_df`.
  trials_modified_ordering = obj$glm_fits_df %>%
    pull(trial_modified) %>%
    as.character()
  # Ordering excluding placebo.
  trials_ordering = obj$glm_fits_df %>%
    filter(trial_modified != "Placebo") %>%
    pull(trial_modified) %>%
    as.character()

  # Extract coefficients of the outcome regression models for Y.
  ncov = length(obj$glm_fits_df$glm_coef_Y[[1]])
  coefs_vec_Y = obj$glm_fits_df %>%
    pull(glm_coef_Y) %>%
    unlist()
  names(coefs_vec_Y) = paste(names(coefs_vec_Y), rep(trials_modified_ordering, each = ncov), sep = " - ")
  # Extract coefficients of the outcome regression models for S. We don't need
  # the model for the placebo group.
  coefs_vec_S = obj$glm_fits_df %>%
    filter(trial_modified != "Placebo") %>%
    pull(glm_coef_S, name = "trial_modified") %>%
    unlist()
  names(coefs_vec_S) = paste(names(coefs_vec_S), rep(trials_ordering, each = ncov), sep = " - ")
  
  # Extract standardized mean estimates for Y and S. 
  standardized_mean_Y = left_join(
    tibble(trial_modified = trials_modified_ordering),
    obj$standardized_means_df) %>%
    pull(mean_Y)
  names(standardized_mean_Y) = paste("standardized_mean_Y", trials_modified_ordering, sep = " - ")
  
  standardized_mean_S = left_join(
    tibble(trial_modified = trials_ordering),
    obj$standardized_means_df) %>%
    pull(mean_S)
  names(standardized_mean_S) = paste("standardized_mean_S", trials_ordering, sep = " - ")
  
  # Extract mean estimates for Y and S for the target trial.
  target_trial_mean_Y = obj$naive_means_df %>%
    filter(trial == obj$target_trial, vax == 1) %>%
    pull(mean_Y)
  names(target_trial_mean_Y) = paste("mean_Y", target_trial, sep = " - ")

  target_trial_mean_S = obj$naive_means_df %>%
    filter(trial == obj$target_trial, vax == 1) %>%
    pull(mean_S)
  names(target_trial_mean_S) = paste("mean_S", target_trial, sep = " - ")
  
  # Extract naive treatment effect estimate for the target trial. 
  target_trial_VE_est = obj$standardized_trt_effects_df %>% filter(trial == obj$target_trial) %>%
    pull(VE_est, name = trial)
  names(target_trial_VE_est) = paste("VE_est", target_trial, sep = " - ")
  
  # Extract standardized VE estimates. 
  standardized_VE_est = left_join(tibble(trial = trials_ordering),
                                  obj$standardized_trt_effects_df) %>%
    pull(VE_est)
  names(standardized_VE_est) = paste("standardized_VE_est", trials_ordering, sep = " - ")
  
  # Extract estimated Pearson correlation and regression slope for the
  # standardized estimates.
  standardized_cor_p = obj$cor_standardized_df %>%
    pull(cor_p)
  names(standardized_cor_p) = "standardized_cor_p"
  standardized_beta = obj$cor_standardized_df %>%
    pull(beta)
  names(standardized_beta) = "standardized_beta"
  
  # Join all estimates.
  estimates_vec = c(
    coefs_vec_Y,
    coefs_vec_S,
    target_trial_mean_Y,
    standardized_mean_Y,
    target_trial_mean_S,
    standardized_mean_S,
    target_trial_VE_est,
    standardized_VE_est,
    standardized_cor_p,
    standardized_beta
  )
  # If the weights are treated as being estimated, they are included in the
  # returned vector of estimates.
  if (estimate_weights) {
    estimated_weights = obj$weights_df %>%
      # The weight for placebo patients is known anyhow, the
      # corresponding weight is thus always treated as known.
      filter(CC_stratum != "Placebo") %>% 
      pull(weight, name = CC_stratum)

    estimates_vec = c(estimates_vec, 1 / estimated_weights)
  }
  return(estimates_vec)
}


