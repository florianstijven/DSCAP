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
              select(-target_trial),
            by = c("setting", "n_trials", "measure", "treatment"),
            relationship = "many-to-one")


# Compute summaries from the simulation results
simulations_summary_tbl = simulations_results_tbl %>%
  filter(measure != "permutation_LRT_p_value") %>%
  group_by(setting, n_trials, n_t, CC_sampling, measure, type, treatment) %>%
  summarise(
    coverage = mean((estimand >= CI_lower_bs) &
                      (estimand <= CI_upper_bs), na.rm = TRUE),
    mean_bias = mean(estimate - estimand, na.rm = TRUE),
    median_bias = median(estimate - estimand, na.rm = TRUE),
    emp_SD = sd(estimate, na.rm = TRUE),
    mean_SE = mean(SE_bs, na.rm = TRUE),
    MSE = mean((estimand - estimate)^2, na.rm = TRUE)
  )

error_rates_permutation_tbl = simulations_results_tbl %>%
  filter(measure %in% c("permutation_LRT_p_value_S", "permutation_LRT_p_value_Y")) %>%
  group_by(setting, n_trials, n_t, CC_sampling, measure) %>%
  summarise(
    type_I_error_rate = mean(estimate < 0.05, na.rm = TRUE)
  )

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
# beta_est), (ii) whether to plot the mean or median bias, (iii) whether to plot
# the results under case-cohort sampling.
plot_bias <- function(measure_var, bias_type = "mean", include_cc = TRUE) {
  data <- simulations_summary_tbl %>%
    filter(measure %in% measure_var,
           setting != "X")


  data <- data %>% filter(CC_sampling == include_cc)


  bias_col <- if (bias_type == "mean") "mean_bias" else "median_bias"
  # Create measure labels as character strings to be parsed
  measure_labels <- c(
    "cor_s_est" = "rho[s]",
    "cor_p_est" = "rho[p]",
    "beta_est" = "rho[beta]"
  )
  
  # Create setting labeller
  setting_labeller <- as_labeller(c("1" = "Setting 1", "2" = "Setting 2", "3" = "Setting 3"))
  
  data %>%
    mutate(measure = factor(measure, 
                            levels = names(measure_labels),
                            labels = unname(measure_labels))) %>%
    ggplot(aes(
      x = n_t,
      y = .data[[bias_col]],
      color = type,
      shape = type
    )) +
    geom_point(position = position_dodge(width = 0.1)) +
    geom_abline(intercept = 0, slope = 0) +
    scale_y_continuous(name = paste(str_to_title(bias_type), "Bias")) +
    scale_x_continuous(name = "Number of Subjects per Trial (thousands)", labels = function(x) x/1000) +
    scale_color_discrete(name = "Estimator") +
    scale_shape_discrete(name = "Estimator") +
    facet_grid(measure ~ setting, labeller = labeller(measure = label_parsed, setting = setting_labeller), scales = "free_y") +
    theme(legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(0, "cm"))
}

# Plot the bias across different settings.
bias_plots <- tibble(
  measure = list(c("cor_s_est", "cor_p_est", "beta_est")),
  plot_mean_cc = map(measure, ~plot_bias(.x, bias_type = "mean", include_cc = TRUE)),
  plot_mean = map(measure, ~plot_bias(.x, bias_type = "mean", include_cc = FALSE)),
  plot_median_cc = map(measure, ~plot_bias(.x, bias_type = "median", include_cc = TRUE)),
  plot_median = map(measure, ~plot_bias(.x, bias_type = "median", include_cc = FALSE))
)

walk2(bias_plots$measure, bias_plots$plot_mean_cc, ~ggsave(
  plot = .y,
  filename = paste0("bias-mean-case-cohort-", paste(.x, collapse = "-"), ".pdf"),
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

walk2(bias_plots$measure, bias_plots$plot_median_cc, ~ggsave(
  plot = .y,
  filename = paste0("bias-median-case-cohort-", paste(.x, collapse = "-"), ".pdf"),
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
plot_mse <- function(measure_var, include_cc = TRUE) {
  data <- simulations_summary_tbl %>%
    filter(measure %in% measure_var,
           setting != "X")
  
  if (include_cc) {
    data <- data %>% filter(CC_sampling == TRUE)
  } else {
    data <- data %>% filter(CC_sampling == FALSE)
  }
  
  # Create measure labels as character strings to be parsed
  measure_labels <- c(
    "cor_s_est" = "rho[s]",
    "cor_p_est" = "rho[p]",
    "beta_est" = "rho[beta]"
  )
  
  # Create setting labeller
  setting_labeller <- as_labeller(c("1" = "Setting 1", "2" = "Setting 2", "3" = "Setting 3"))
  
  data %>%
    mutate(measure = factor(measure,
                            levels = names(measure_labels),
                            labels = unname(measure_labels))) %>%
    ggplot(aes(
      x = n_t,
      y = MSE,
      color = type,
      shape = type
    )) +
    geom_point(position = position_dodge(width = 0.1)) +
    scale_y_continuous(name = "MSE", transform = "log10") +
    scale_x_continuous(name = "Number of Subjects per Trial (thousands)", labels = function(x) x/1000) +
    scale_color_discrete(name = "Estimator") +
    scale_shape_discrete(name = "Estimator") +
    facet_grid(measure ~ setting, scales = "free_y", labeller = labeller(measure = label_parsed, setting = setting_labeller)) +
    theme(legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(0, "cm"))
}

# Plot the MSE across different settings.
mse_plots <- tibble(
  measure = list(c("cor_s_est", "cor_p_est", "beta_est")),
  plot_cc = map(measure, ~plot_mse(.x, include_cc = TRUE)),
  plot_no_cc = map(measure, ~plot_mse(.x, include_cc = FALSE))
)

walk2(mse_plots$measure, mse_plots$plot_cc, ~ggsave(
  plot = .y,
  filename = paste0("mse-case-cohort-", paste(.x, collapse = "-"), ".pdf"),
  path = figures_dir,
  height = double_height,
  width = double_width,
  dpi = res,
  device = "pdf",
  units = "cm"
))

walk2(mse_plots$measure, mse_plots$plot_no_cc, ~ggsave(
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
plot_coverage <- function(measure_var, include_cc = TRUE) {
  coverage_limits = c(0.65, 1)
  
  data <- simulations_summary_tbl %>%
    mutate(type = as.factor(type)) %>%
    filter(measure %in% measure_var,
           setting != "X",
           type != "naive")
  
data <- data %>% filter(CC_sampling == include_cc)

  
  # Create measure labels as character strings to be parsed
  measure_labels <- c(
    "cor_s_est" = "rho[s]",
    "cor_p_est" = "rho[p]",
    "beta_est" = "rho[beta]"
  )
  
  # Create setting labeller
  setting_labeller <- as_labeller(c("1" = "Setting 1", "2" = "Setting 2", "3" = "Setting 3"))
  
  data %>%
    mutate(measure = factor(measure,
                            levels = names(measure_labels),
                            labels = unname(measure_labels))) %>%
    ggplot(aes(
      x = n_t,
      y = coverage,
      color = type,
      shape = type
    )) +
    geom_point(position = position_dodge(width = 0.1)) +
    geom_abline(intercept = 0.95, slope = 0) +
    scale_y_continuous(name = "Coverage", limits = coverage_limits) +
    scale_color_discrete(name = "Estimator", drop = FALSE, breaks = unique(data$type)) +
    scale_x_continuous(name = "Number of Subjects per Trial (thousands)", labels = function(x) x/1000) +
    scale_shape_discrete(name = "Estimator", drop = FALSE, breaks = unique(data$type)) +
    facet_grid(measure ~ setting, scales = "free_y", labeller = labeller(measure = label_parsed, setting = setting_labeller)) +
    theme(legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(0, "cm"))
}

# Plot the coverage across different settings.
coverage_plots <- tibble(
  measure = list(c("cor_s_est", "cor_p_est", "beta_est")),
  plot_cc = map(measure, ~plot_coverage(.x, include_cc = TRUE)),
  plot_no_cc = map(measure, ~plot_coverage(.x, include_cc = FALSE))
)

walk2(coverage_plots$measure, coverage_plots$plot_cc, ~ggsave(
  plot = .y,
  filename = paste0("coverage-case-cohort-", paste(.x, collapse = "-"), ".pdf"),
  path = figures_dir,
  height = double_height,
  width = double_width,
  dpi = res,
  device = "pdf",
  units = "cm"
))

walk2(coverage_plots$measure, coverage_plots$plot_no_cc, ~ggsave(
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

# With case-cohort sampling
sink(file = paste0(tables_dir, "simulation-results-summary-case-cohort.txt"))
simulations_summary_tbl %>%
  filter(measure %in% c("cor_s_est", "cor_p_est", "beta_est"),
         CC_sampling == TRUE) %>%
  select(-treatment) %>%
  arrange(setting, measure, n_trials, n_t) %>%
  print(n = 500)
sink()

# Without case-cohort sampling
sink(file = paste0(tables_dir, "simulation-results-summary.txt"))
simulations_summary_tbl %>%
  filter(measure %in% c("cor_s_est", "cor_p_est", "beta_est"),
         CC_sampling == FALSE) %>%
  select(-treatment) %>%
  arrange(setting, measure, n_trials, n_t) %>%
  print(n = 500)
sink()

# Print the error rates for the permutation LRT. 
sink(file = paste0(tables_dir, "LRT-results.txt"))
error_rates_permutation_tbl %>%
  print(n = 500)
sink()



