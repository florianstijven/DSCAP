# Preliminaries ----------------------------------------------------
# load packages
library(tidyverse)
library(ggpubr)

# Specify options for saving the plots to files
figures_dir = "results/simulations/figures/"
tables_dir = "results/simulations/tables/"

# Read in simulation results
simulations_results_tbl <- readRDS(
  file = "results/simulations/raw-results/simulations_results_tbl.rds"
)

# Add true values to the results data set. We need the true values to compute
# the bias, MSE, and coverage of the confidence intervals.
true_values_tbl <- readRDS(
  file = "results/simulations/raw-results/true_values_tbl.rds"
) %>%
  select(-type) %>%
  rename(estimand = estimate)

# Join the simulations results with the estimands.
simulations_results_tbl <- simulations_results_tbl %>%
  left_join(true_values_tbl %>%
              filter(target_trial == "standardized") %>%
              select(-target_trial, -different_placebo_model),
            by = c("setting", "n_trials", "measure", "treatment"),
            relationship = "many-to-one")


# Compute summaries from the simulation results
simulations_summary_tbl = simulations_results_tbl %>%
  filter(!(stringr::str_detect(measure, "p_value"))) %>%
  group_by(setting, n_trials, n_t, CC_sampling, measure, type, treatment, model_O_correct, model_T_correct) %>%
  summarise(
    coverage = mean((estimand >= CI_lower_bs) &
                      (estimand <= CI_upper_bs), na.rm = TRUE),
    mean_bias = mean(estimate - estimand, na.rm = TRUE),
    median_bias = median(estimate - estimand, na.rm = TRUE),
    emp_SD = sd(estimate, na.rm = TRUE),
    mean_SE = mean(SE_bs, na.rm = TRUE),
    MSE = mean((estimand - estimate)^2, na.rm = TRUE)
  )

# We can drop the entries for the IPW and standardization estimators when the
# model on which the estimator does not depend is misspecified. For example, the
# IPW estimator does not depend on the model for the outcome, so we can drop the
# entries for the IPW estimator when the outcome model is misspecified. 
simulations_summary_tbl = simulations_summary_tbl %>%
  filter(
    !(type == "ipw" & model_O_correct == FALSE),
    !(type == "standardization" & model_T_correct == FALSE),
    !(type == "naive" & (model_O_correct == FALSE | model_T_correct == FALSE))
  )

# We rename the scenarios. Scenarios 1--3 without CC sampling and correctly
# specified models correspond to the main scenarios. Scenarios 3 with CC is
# renamed to `3 (CC)`. 
simulations_summary_tbl = simulations_summary_tbl %>%
  mutate(setting = case_when(
    setting == "3" & CC_sampling == TRUE ~ "3 (CC)",
    TRUE ~ as.character(setting)
  ))

error_rates_LRT_tbl = simulations_results_tbl %>%
  filter(str_detect(measure, "p_value"),,
         # We only consider the LRT for correctly specified models.
         model_O_correct,
         model_T_correct) %>%
  group_by(setting, n_trials, n_t, CC_sampling, measure) %>%
  summarise(type_I_error_rate = mean(estimate < 0.05, na.rm = TRUE))

# Summary of the Simulation Results ---------------------------------------

## True Trial-Level Effects -----------------------------------------------

# Plot the true trial-level effects (unstandardized and standardized) for each
# setting and treatment arm. 
true_effects_plot = true_values_tbl %>%
  filter(measure %in% c("VE_est", "mean_diff_S_est"),
         setting != "X") %>%
  pivot_wider(names_from = "measure", values_from = "estimand") %>%
  ggplot(aes(
    x = mean_diff_S_est,
    y = VE_est,
    colour = as.factor(treatment)
  )) +
  geom_point() +
  scale_y_continuous(name = "VE") +
  scale_x_continuous(name = "Mean Difference in S") +
  scale_color_discrete(name = "Trial") +
  facet_grid(target_trial~setting, scales = "free", labeller = labeller(setting = as_labeller(c("1" = "Setting 1", "2" = "Setting 2", "3" = "Setting 3")))) +  
  theme(legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(0, "cm"))

ggsave(
  plot = true_effects_plot,
  filename = "ma-plots.pdf",
  path = figures_dir,
  height = double_height,
  width = double_width,
  dpi = res,
  device = "pdf",
  units = "cm"
)


## Estimation Accuracies --------------------------------------------------

### Bias ------------------------------------------------------------------

# Function that plots the bias of the estimators for the DSCAPs with the
# following options: (i) which DSCAP to plot (cor_s_est, cor_p_est, or
# beta_est), (ii) whether to plot the mean or median bias, and (iii) whether to
# plot the results under misspecification (only setting 1).
plot_bias <- function(measure_var, bias_type = "mean", misspecification = FALSE) {
  data <- simulations_summary_tbl %>%
    filter(measure %in% measure_var,
           setting != "X")
  
   if (misspecification) {
    data <- data %>% filter(setting == "1")
  } else {
    data <- data %>% filter(model_O_correct == TRUE, model_T_correct == TRUE)
  }


  bias_col <- if (bias_type == "mean") "mean_bias" else "median_bias"
  # Create measure labels as character strings to be parsed
  measure_labels <- c(
    "cor_s_est" = "rho[s]",
    "cor_p_est" = "rho[p]",
    "beta_est" = "rho[beta]"
  )
  
  # Create setting labeller
  setting_labeller <- as_labeller(c("1" = "Setting 1", "2" = "Setting 2", "3" = "Setting 3", "3 (CC)" = "Setting 3 (CC)"))

  # We need custom labels for the misspecified results.
   if (misspecification) {
     data$misspecification_labels = case_when(
       data$model_O_correct == FALSE & data$model_T_correct == FALSE ~ "O and T misspecified",
       data$model_O_correct == FALSE & data$model_T_correct == TRUE ~ "O misspecified",
       data$model_O_correct == TRUE & data$model_T_correct == FALSE ~ "T misspecified",
       TRUE ~ "Both correct"
     )
   }
  
  if (misspecification) {
    ggplot_temp = data %>%
      mutate(measure = factor(measure, 
                              levels = names(measure_labels),
                              labels = unname(measure_labels))) %>%
      ggplot(aes(
        x = measure,
        y = .data[[bias_col]],
        color = type,
        shape = type
      )) +
      scale_x_discrete(labels = function(x) parse(text = x)) +
      facet_grid(. ~ misspecification_labels, scales = "free_y")
  } else {
    ggplot_temp = data %>%
      mutate(measure = factor(measure, 
                              levels = names(measure_labels),
                              labels = unname(measure_labels))) %>%
      ggplot(aes(
        x = n_t,
        y = .data[[bias_col]],
        color = type,
        shape = type
      )) +
      scale_x_continuous(name = "Number of Subjects per Trial (thousands)", labels = function(x) x/1000) +
      facet_grid(measure ~ setting, labeller = labeller(measure = label_parsed, setting = setting_labeller), scales = "free_y")
  }
  
  ggplot_temp = ggplot_temp +
    geom_point(position = position_dodge(width = 0.1)) +
    geom_abline(intercept = 0, slope = 0) +
    scale_y_continuous(name = paste(str_to_title(bias_type), "Bias")) +
    scale_color_discrete(name = "Estimator") +
    scale_shape_discrete(name = "Estimator") +
    theme(legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(0, "cm"))
  
  return(ggplot_temp)
}

# Plot the bias across different settings.
bias_plots <- tibble(
  measure = list(c("cor_s_est", "cor_p_est", "beta_est")),
  plot_mean_misspecification = map(measure, ~plot_bias(.x, bias_type = "mean", misspecification = TRUE)),
  plot_median_misspecification = map(measure, ~plot_bias(.x, bias_type = "median", misspecification = TRUE)),
  plot_mean = map(measure, ~plot_bias(.x, bias_type = "mean", misspecification = FALSE)),
  plot_median = map(measure, ~plot_bias(.x, bias_type = "median", misspecification = FALSE))
)

walk2(bias_plots$measure, bias_plots$plot_mean_misspecification, ~ggsave(
  plot = .y,
  filename = paste0("bias-mean-misspecification-", paste(.x, collapse = "-"), ".pdf"),
  path = figures_dir,
  height = double_height,
  width = double_width,
  dpi = res,
  device = "pdf",
  units = "cm"
))

walk2(bias_plots$measure, bias_plots$plot_median_misspecification, ~ggsave(
  plot = .y,
  filename = paste0("bias-median-misspecification-", paste(.x, collapse = "-"), ".pdf"),
  path = figures_dir,
  height = double_height,
  width = double_width,
  dpi = res,
  device = "pdf",
  units = "cm"
))

walk2(bias_plots$measure, bias_plots$plot_mean, ~ggsave(
  plot = .y,
  filename = paste0("bias-mean-", paste(.x, collapse = "-"), ".pdf"),
  path = figures_dir,
  height = double_height,
  width = double_width,
  dpi = res,
  device = "pdf",
  units = "cm"
))

walk2(bias_plots$measure, bias_plots$plot_median, ~ggsave(
  plot = .y,
  filename = paste0("bias-median-", paste(.x, collapse = "-"), ".pdf"),
  path = figures_dir,
  height = double_height,
  width = double_width,
  dpi = res,
  device = "pdf",
  units = "cm"
))

### MSE -------------------------------------------------------------------

# Function that plots the MSE of the estimators for the DSCAPs with the
# following options: (i) which DSCAP to plot (cor_s_est, cor_p_est, or
# beta_est), (ii) whether to plot the results under case-cohort sampling.
plot_mse <- function(measure_var, bias_type = "mean", misspecification = FALSE) {
  data <- simulations_summary_tbl %>%
    filter(measure %in% measure_var,
           setting != "X")
  
  if (misspecification) {
    data <- data %>% filter(setting == "1")
  } else {
    data <- data %>% filter(model_O_correct == TRUE, model_T_correct == TRUE)
  }
  
  
  # Create measure labels as character strings to be parsed
  measure_labels <- c(
    "cor_s_est" = "rho[s]",
    "cor_p_est" = "rho[p]",
    "beta_est" = "rho[beta]"
  )
  
  # Create setting labeller
  setting_labeller <- as_labeller(c("1" = "Setting 1", "2" = "Setting 2", "3" = "Setting 3", "3 (CC)" = "Setting 3 (CC)"))
  
  # We need custom labels for the misspecified results.
  if (misspecification) {
    data$misspecification_labels = case_when(
      data$model_O_correct == FALSE & data$model_T_correct == FALSE ~ "O and T misspecified",
      data$model_O_correct == FALSE & data$model_T_correct == TRUE ~ "O misspecified",
      data$model_O_correct == TRUE & data$model_T_correct == FALSE ~ "T misspecified",
      TRUE ~ "Both correct"
    )
  }
  
  if (misspecification) {
    ggplot_temp = data %>%
      mutate(measure = factor(measure, 
                              levels = names(measure_labels),
                              labels = unname(measure_labels))) %>%
      ggplot(aes(
        x = measure,
        y = MSE,
        color = type,
        shape = type
      )) +
      scale_x_discrete(labels = function(x) parse(text = x)) +
      facet_grid(. ~ misspecification_labels, scales = "free_y")
  } else {
    ggplot_temp = data %>%
      mutate(measure = factor(measure, 
                              levels = names(measure_labels),
                              labels = unname(measure_labels))) %>%
      ggplot(aes(
        x = n_t,
        y = MSE,
        color = type,
        shape = type
      )) +
      scale_x_continuous(name = "Number of Subjects per Trial (thousands)", labels = function(x) x/1000) +
      facet_grid(measure ~ setting, labeller = labeller(measure = label_parsed, setting = setting_labeller), scales = "free_y")
  }
  
  ggplot_temp = ggplot_temp +
    geom_point(position = position_dodge(width = 0.1)) +
    geom_abline(intercept = 0, slope = 0) +
    scale_y_continuous(name = "MSE") +
    scale_color_discrete(name = "Estimator") +
    scale_shape_discrete(name = "Estimator") +
    theme(legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(0, "cm"))
  
  return(ggplot_temp)
}

# Plot the MSE across different settings.
mse_plots <- tibble(
  measure = list(c("cor_s_est", "cor_p_est", "beta_est")),
  plot = map(measure, ~plot_mse(.x, misspecification = FALSE)),
  plot_misspecification = map(measure, ~plot_mse(.x, misspecification = TRUE))
)

walk2(mse_plots$measure, mse_plots$plot_misspecification, ~ggsave(
  plot = .y,
  filename = paste0("mse-misspecification-", paste(.x, collapse = "-"), ".pdf"),
  path = figures_dir,
  height = double_height,
  width = double_width,
  dpi = res,
  device = "pdf",
  units = "cm"
))

walk2(mse_plots$measure, mse_plots$plot, ~ggsave(
  plot = .y,
  filename = paste0("mse-", paste(.x, collapse = "-"), ".pdf"),
  path = figures_dir,
  height = double_height,
  width = double_width,
  dpi = res,
  device = "pdf",
  units = "cm"
))


## Performance of Inferences ----------------------------------------------

### Coverage --------------------------------------------------------------

# Function that plots the coverage of the confidence intervals for the DSCAPs 
# with the following options: (i) which DSCAP to plot (cor_s_est, cor_p_est, or
# beta_est), (ii) whether to plot the results under case-cohort sampling.
plot_coverage <- function(measure_var, misspecification = FALSE) {
  coverage_limits = c(0.65, 1)
  
  data <- simulations_summary_tbl %>%
    mutate(type = as.factor(type)) %>%
    filter(measure %in% measure_var,
           setting != "X",
           type != "naive")
  
  if (misspecification) {
    data <- data %>% filter(setting == "1")
  } else {
    data <- data %>% filter(model_O_correct == TRUE, model_T_correct == TRUE)
  }
  
  
  # Create measure labels as character strings to be parsed
  measure_labels <- c(
    "cor_s_est" = "rho[s]",
    "cor_p_est" = "rho[p]",
    "beta_est" = "rho[beta]"
  )
  
  # Create setting labeller
  setting_labeller <- as_labeller(c("1" = "Setting 1", "2" = "Setting 2", "3" = "Setting 3", "3 (CC)" = "Setting 3 (CC)"))
  
  # We need custom labels for the misspecified results.
  if (misspecification) {
    data$misspecification_labels = case_when(
      data$model_O_correct == FALSE & data$model_T_correct == FALSE ~ "O and T misspecified",
      data$model_O_correct == FALSE & data$model_T_correct == TRUE ~ "O misspecified",
      data$model_O_correct == TRUE & data$model_T_correct == FALSE ~ "T misspecified",
      TRUE ~ "Both correct"
    )
  }
  
  if (misspecification) {
    data$misspecification_labels = case_when(
      data$model_O_correct == FALSE & data$model_T_correct == FALSE ~ "O and T misspecified",
      data$model_O_correct == FALSE & data$model_T_correct == TRUE ~ "O misspecified",
      data$model_O_correct == TRUE & data$model_T_correct == FALSE ~ "T misspecified",
      TRUE ~ "Both correct"
    )
  }
  
  if (misspecification) {
    ggplot_temp = data %>%
      mutate(measure = factor(measure, 
                              levels = names(measure_labels),
                              labels = unname(measure_labels))) %>%
      ggplot(aes(
        x = measure,
        y = coverage,
        color = type,
        shape = type
      )) +
      scale_x_discrete(labels = function(x) parse(text = x)) +
      facet_grid(. ~ misspecification_labels, scales = "free_y")
  } else {
    ggplot_temp = data %>%
      mutate(measure = factor(measure, 
                              levels = names(measure_labels),
                              labels = unname(measure_labels))) %>%
      ggplot(aes(
        x = n_t,
        y = coverage,
        color = type,
        shape = type
      )) +
      scale_x_continuous(name = "Number of Subjects per Trial (thousands)", labels = function(x) x/1000) +
      facet_grid(measure ~ setting, labeller = labeller(measure = label_parsed, setting = setting_labeller), scales = "free_y")
  }
  
  ggplot_temp = ggplot_temp +
    geom_point(position = position_dodge(width = 0.1)) +
    geom_abline(intercept = 0.95, slope = 0) +
    scale_y_continuous(name = "Coverage", limits = coverage_limits) +
    scale_color_discrete(name = "Estimator", drop = FALSE, breaks = unique(data$type)) +
    scale_shape_discrete(name = "Estimator", drop = FALSE, breaks = unique(data$type)) +
    theme(legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(0, "cm"))
  
  return(ggplot_temp)
}

# Plot the coverage across different settings.
coverage_plots <- tibble(
  measure = list(c("cor_s_est", "cor_p_est", "beta_est")),
  plot_misspecification = map(measure, ~plot_coverage(.x, misspecification = TRUE)),
  plot = map(measure, ~plot_coverage(.x, misspecification = FALSE))
)

walk2(coverage_plots$measure, coverage_plots$plot_misspecification, ~ggsave(
  plot = .y,
  filename = paste0("coverage-misspecification-", paste(.x, collapse = "-"), ".pdf"),
  path = figures_dir,
  height = double_height,
  width = double_width,
  dpi = res,
  device = "pdf",
  units = "cm"
))

walk2(coverage_plots$measure, coverage_plots$plot, ~ggsave(
  plot = .y,
  filename = paste0("coverage-", paste(.x, collapse = "-"), ".pdf"),
  path = figures_dir,
  height = double_height,
  width = double_width,
  dpi = res,
  device = "pdf",
  units = "cm"
))

# Tables ---------------------------------------------------------------

## True Values and Event Proportions ------------------------------------------------
options(width = 1000) 
# Print the true trial-level standardized effects.
sink(file = paste0(tables_dir, "true-trial-level-effects.txt"))
true_values_tbl %>%
  filter(measure %in% c("VE_est", "mean_diff_S_est")) %>%
  pivot_wider(names_from = "measure", values_from = "estimand") %>%
  print(n = 500)
sink()

# Print the true causal association parameters.
sink(file = paste0(tables_dir, "true-causal-association-parameters.txt"))
true_values_tbl %>%
  filter(measure %in% c("cor_s_est", "cor_p_est", "beta_est")) %>%
  pivot_wider(names_from = "measure", values_from = "estimand") %>%
  select(-treatment, -trial) %>%
  print(n = 500)
sink()

# Print the event proportions.
sink(file = paste0(tables_dir, "event-proportions.txt"))
true_values_tbl %>%
  filter(measure == "prop_events", target_trial == "naive") %>%
  select(setting, treatment, trial, estimand) %>%
  mutate(treatment = ifelse(treatment == 0, "Control (0)", "Active (1)")) %>%
  pivot_wider(names_from = "treatment", values_from = "estimand") %>%
  print(n = 500)
sink()

## Estimators and Inferences --------------------------------------------------------------

# Print the performance measures for the estimators including trial-level
# standardized effects and DSCAPs.

sink(file = paste0(tables_dir, "simulation-results-summary.txt"))
simulations_summary_tbl %>%
  filter(measure %in% c("cor_s_est", "cor_p_est", "beta_est")) %>%
  arrange(setting, measure, n_trials, n_t) %>%
  ungroup() %>%
  print(n = 500)
sink()

# Print the simulations results in a formatted tables used in the main text.
sink(file = paste0(tables_dir, "formatted-simulation-results-summary.txt"))
# For the naive estimator, we only need on column with the bias. For the other
# standardized estimators, we need a column with the bias and a column with the
# coverage.
cat("** Simulations results for scenarios 1--4 with correctly specified models. **\n\n")
simulations_summary_tbl %>%
  filter(measure %in% c("cor_s_est", "cor_p_est", "beta_est")) %>%
  filter(model_O_correct == TRUE, model_T_correct == TRUE, setting != "X") %>%
  select(setting,
         n_trials,
         n_t,
         CC_sampling,
         measure,
         type,
         mean_bias,
         coverage) %>%
  pivot_wider(names_from = "type",
              values_from = c("mean_bias", "coverage")) %>%
  arrange(setting, measure, n_trials, n_t) %>%
  ungroup() %>%
  select(-treatment, -model_O_correct, -n_trials, -CC_sampling) %>%
  select(
    setting,
    n_t,
    measure,
    mean_bias_naive,
    `mean_bias_doubly robust`,
    `coverage_doubly robust`,
    `mean_bias_ipw`,
    `coverage_ipw`,
    `mean_bias_standardization`,
    `coverage_standardization`
  ) %>%
  # Round the numeric columns to 3 decimal places for better formatting.
  mutate(across(where(is.numeric), ~ round(., 3))) %>%
  print(n = 500)
cat("\n\n** Simulations results for scenario 1 with misspecified models. **\n\n")
simulations_summary_tbl %>%
  filter(measure %in% c("cor_s_est", "cor_p_est", "beta_est")) %>%
  filter(setting == "1") %>%
  select(setting, n_trials, n_t, CC_sampling, measure, type, model_O_correct, model_T_correct, mean_bias, coverage) %>%
  pivot_wider(names_from = "type", values_from = c("mean_bias", "coverage")) %>%
  arrange(setting, measure, n_trials, n_t) %>%
  ungroup() %>%
  select(-treatment, -n_trials, -CC_sampling) %>%
  select(
    setting,
    n_t,
    measure,
    model_O_correct,
    model_T_correct,
    mean_bias_naive,
    `mean_bias_doubly robust`,
    `coverage_doubly robust`,
    `mean_bias_ipw`,
    `coverage_ipw`,
    `mean_bias_standardization`,
    `coverage_standardization`
  ) %>%
  # Round the numeric columns to 3 decimal places for better formatting.
  mutate(across(where(is.numeric), ~ round(., 3))) %>%
  arrange(-model_T_correct,-model_O_correct) %>%
  print(n = 500)
sink()



# Print the error rates for the LRT.
sink(file = paste0(tables_dir, "LRT-results.txt"))
error_rates_LRT_tbl %>%
  pivot_wider(names_from = "measure", values_from = "type_I_error_rate") %>%
  print(n = 500)
sink()



