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
  select(-n_t, -CC_sampling, -type) %>%
  rename(estimand = estimate)

# Join the simulations results wit the estimands.
simulations_results_tbl <- simulations_results_tbl %>%
  left_join(true_values_tbl,
            by = c("setting", "n_trials", "measure", "treatment"))


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
  filter(measure == "permutation_LRT_p_value") %>%
  group_by(setting, n_trials, n_t, CC_sampling) %>%
  summarise(
    type_I_error_rate = mean(estimate < 0.05, na.rm = TRUE)
  )

# Summary of the Simulation Results ---------------------------------------

## Estimation Accuracies --------------------------------------------------

### Bias ------------------------------------------------------------------

mean_bias = simulations_summary_tbl %>%
  filter(measure %in% c("cor_s_est", "cor_p_est", "beta_est")) %>%
  ggplot(aes(
    x = n_t,
    y = mean_bias,
    color = CC_sampling
  )) +
  geom_point(position = position_dodge(width = 0.1)) +
  geom_line() +
  geom_abline(intercept = 0, slope = 0) +
  scale_y_continuous(name = "Mean Bias") +
  scale_color_discrete(name = "Case Cohort Design") +
  facet_grid(measure ~ setting, scales = "free_y") + 
  theme(legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(0, "cm"))

ggsave(
  plot = mean_bias,
  filename = "mean-bias.pdf",
  path = figures_dir,
  height = double_height,
  width = double_width,
  dpi = res,
  device = "pdf",
  units = "cm"
)

median_bias = simulations_summary_tbl %>%
  filter(measure %in% c("cor_s_est", "cor_p_est", "beta_est")) %>%
  ggplot(aes(
    x = n_t,
    y = median_bias,
    color = CC_sampling
  )) +
  geom_point(position = position_dodge(width = 0.1)) +
  geom_line() +
  geom_abline(intercept = 0, slope = 0) +
  scale_y_continuous(name = "Median Bias") +
  scale_color_discrete(name = "Case Cohort Design") +
  facet_grid(measure ~ setting, scales = "free_y") + 
  theme(legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(0, "cm"))

ggsave(
  plot = median_bias,
  filename = "median-bias.pdf",
  path = figures_dir,
  height = double_height,
  width = double_width,
  dpi = res,
  device = "pdf",
  units = "cm"
)

### MSE -------------------------------------------------------------------

mse_plot = simulations_summary_tbl %>%
  filter(measure %in% c("cor_s_est", "cor_p_est", "beta_est")) %>%
  ggplot(aes(
    x = n_t,
    y = MSE,
    color = CC_sampling
  )) +
  geom_point(position = position_dodge(width = 0.1)) +
  geom_line() +
  scale_y_continuous(name = "MSE", transform = "log10") +
  scale_color_discrete(name = "Case Cohort Design") +
  facet_grid(measure ~ setting, scales = "free_y") + 
  theme(legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(0, "cm"))

ggsave(
  plot = mse_plot,
  filename = "mse.pdf",
  path = figures_dir,
  height = double_height,
  width = double_width,
  dpi = res,
  device = "pdf",
  units = "cm"
)


## Performance of Inferences ----------------------------------------------

### Coverage --------------------------------------------------------------
coverage_limits = c(0.65, 1)

coverage_bs = simulations_summary_tbl %>%
  filter(measure %in% c("cor_s_est", "cor_p_est", "beta_est")) %>%
  ggplot(aes(
    x = n_t,
    y = coverage,
    color = CC_sampling
  )) +
  geom_point(position = position_dodge(width = 0.1)) +
  geom_line() +
  geom_abline(intercept = 0.95, slope = 0) +
  scale_y_continuous(name = "Coverage", limits = coverage_limits) +
  scale_color_discrete(name = "Case Cohort Design") +
  facet_grid(measure ~ setting, scales = "free_y") + 
  theme(legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(0, "cm"))

ggsave(
  plot = coverage_bs,
  filename = "coverage.pdf",
  path = figures_dir,
  height = double_height,
  width = double_width,
  dpi = res,
  device = "pdf",
  units = "cm"
)

# Tables ---------------------------------------------------------------

sink(file = paste0(tables_dir, "simulation-results-summary.txt"))
simulations_summary_tbl %>%
  filter(measure %in% c("cor_s_est", "cor_p_est", "beta_est")) %>%
  arrange(setting, measure, n_trials, n_t) %>%
  print(n = 500)
sink()

sink(file = paste0(tables_dir, "LRT-results.txt"))
error_rates_permutation_tbl
sink()
