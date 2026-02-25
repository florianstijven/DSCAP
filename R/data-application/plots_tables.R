# Setup -------------------------------------------------------------------
# Load required packages
library(tidyverse)
library(patchwork)

# Specify options for saving the plots to files
figures_dir = "results/figures/"
tables_dir = "results/tables/"
raw_results_dir = "results/raw-results/"

# Source functions that help for reading in the raw-result files.
source("R/helper-functions/processing_results.R")

## Analysis Parameters ---------------------------------------------------

# Read in argument that determines the target population.
args = commandArgs(trailingOnly=TRUE)
target <- args[1]

# Models of analysis: 1 corresponds to 8-trial analysis, 2 for the analysis with 
# J&J (Colombia) and J&J (Brazil) left out. 
modn <- 1:2
# Type of surrogate endpoint.
surr_type <- c("spike","neut")
# Indicator whether weights were treated as estimated.
estimate_weights = 1
# Number of covariates in the regression models used for standardization.
nX = 5
# Truncated population
truncation = 0:1

# Combine all parameters for the analyses/plots into one data set.
plot_parameters_tbl = expand_grid(
  target,
  modn, 
  surr_type,
  estimate_weights,
  truncation
)

# We only have truncation in combination with the five-trial analysis.
plot_parameters_tbl = plot_parameters_tbl %>%
  filter(truncation == (modn - 1))


target_trial = switch(
  target,
  jjsa = "J&J (S. America)",
  AZ = "AstraZeneca",
  jjbr = "J&J (Brazil)",
  jjusa  = "J&J (USA)",
  jjcol = "J&J (Colombia)",
  nov = "Novavax"
)

# Reading results ----------------------------------------------------------

# Read in the raw data files and combine all results into a single data set. 
all_results_tbl = plot_parameters_tbl %>%
  rowwise(everything()) %>%
  reframe(read_results(target, modn, surr_type, estimate_weights, nX, truncation))

# Add extra variable indicating the J&J trials.
all_results_tbl = all_results_tbl %>%
  mutate(jnj_trial = ifelse(str_detect(trial, "J&J"), TRUE, FALSE))

# Plotting -------------------------------------------------------------

# Function to make MA plots.
ma_plot = function(modn = 1, bootstrap = FALSE, estimate_weights = TRUE, target, truncation) {
  trials_chr = c(
    "AstraZeneca",
    "J&J (Brazil)",
    "J&J (Colombia)",
    "J&J (S. Africa)",
    "J&J (S. America)",
    "J&J (USA)",
    "Moderna",
    "Novavax"
  )
  
  plotting_data = left_join(
    all_results_tbl %>%
      select(
        c(
          "estimate",
          "estimand",
          "type",
          "trial",
          "jnj_trial",
          "target",
          "modn",
          "surr_type",
          "estimate_weights",
          "truncation"
        )
      ) %>%
      filter(estimand %in% c("mean_diff_S", "VE")) %>%
      pivot_wider(names_from = "estimand", values_from = "estimate"),
    all_results_tbl %>%
      select(
        c(
          "CI_lower",
          "CI_upper",
          "CI_lower_bs",
          "CI_upper_bs",
          "estimand",
          "type",
          "trial",
          "jnj_trial",
          "target",
          "modn",
          "surr_type",
          "estimate_weights",
          "truncation"
        )
      ) %>%
      filter(estimand %in% c("mean_diff_S", "VE")) %>%
      pivot_longer(cols = c(
        "CI_lower", "CI_upper", "CI_lower_bs", "CI_upper_bs"
      )) %>%
      pivot_wider(
        names_from = c("estimand", "name"),
        values_from = "value"
      )
  ) %>%
    filter(
      type %in% c("naive", "standardized"),
      modn == modn,
      !is.na(trial),
      estimate_weights == .env$estimate_weights,
      target == .env$target,
      truncation == .env$truncation
    )  %>%
    mutate()
  
  # Use bootstrap CIs if requested.
  if (bootstrap) {
    plotting_data = plotting_data %>%
      mutate(
        VE_CI_lower = VE_CI_lower_bs,
        VE_CI_upper = VE_CI_upper_bs,
        mean_diff_S_CI_lower = mean_diff_S_CI_lower_bs,
        mean_diff_S_CI_upper = mean_diff_S_CI_upper_bs
      )
  }
  
  # Limits for the y-axis depend on the settings. The lower limit for the y-axis
  # is set to a negative value in the 8-trial analyses because the standardized
  # estimates may be negative there. 
  if (modn == 1) y_lims = c(-0.5, 1) else y_lims = c(0, 1)
  # Plot for binding Ab.
  ma_ggplot = plotting_data %>%
    filter(modn == .env$modn) %>%
    mutate(surr_type = ifelse(surr_type == "spike", "Binding Ab", "Neutralizing Ab"),
           type = ifelse(type == "naive", "Unstandardized", "Standardized"),
           type = factor(type, levels = c("Unstandardized", "Standardized"))) %>%
    ggplot(aes(
      x = mean_diff_S,
      y = VE,
      color = trial,
      shape = trial
    )) +
    geom_errorbar(
      aes(ymin = VE_CI_lower, ymax = VE_CI_upper),
      width = 0.05,
      color = "darkgrey"
    ) +
    geom_errorbar(
      aes(xmin = mean_diff_S_CI_lower, xmax = mean_diff_S_CI_upper),
      width = 0.05,
      color = "darkgrey"
    ) +
    geom_point(size = 2, show.legend = TRUE) +
    geom_hline(yintercept = 1,
               linetype = "dashed",
               color = "lightgrey") +
    # ylim(0, 1) + 
    coord_cartesian(ylim = y_lims,
                    xlim = c(0, 3)) +
    scale_shape_manual(
      name = "Trial",
      breaks = trials_chr,
      values = stringr::str_detect(trials_chr, "J&J") + 16
    ) +
    scale_color_manual(
      name = "Trial",
      breaks = trials_chr,
      values = viridisLite::viridis(8)
    ) +
    xlab(bquote(tilde( ~ delta) ^ a)) +
    ylab(expression(tilde( ~ tau) ^ a ~ " (VE)")) +
    facet_grid(surr_type~type) +
    labs(shape = "Trial") +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(nrow = 2, byrow = TRUE),
           shape = guide_legend(nrow = 2, byrow = TRUE))
  return(ma_ggplot)
}

# Save all plots.
plot_parameters_tbl %>%
  cross_join(tibble(bootstrap = c(TRUE, FALSE))) %>%
  mutate(temp = purrr::pmap(
    .l = list(
      target = target,
      modn = modn,
      estimate_weights = estimate_weights,
      bootstrap = bootstrap,
      truncation = truncation
    ),
    .f = function(target,
                  modn,
                  estimate_weights,
                  bootstrap, 
                  truncation) {
      if (estimate_weights)
        weights_chr = "estwts"
      else
        weights_chr = "kwnwts"
      if (bootstrap)
        bs_chr = "_bs"
      else
        bs_chr = ""
      outfile = paste0("ma_plot_",
                       target,
                       "_",
                       truncation,
                       "_",
                       "M",
                       modn,
                       "_",
                       weights_chr,
                       bs_chr,
                       ".pdf")
      ma_plot(modn, bootstrap, estimate_weights, target, truncation)
      ggsave(
        filename = outfile,
        device = "pdf",
        path =
          figures_dir,
        width = double_width,
        height = double_height, units = unit
      )
    }
  ))

# Tables --------------------------------------------------------------

## Analysis-Specific Tables -------------------------------------------

# Helper function to save tables for trial-level estimates given `surr_type`,
# `target`, and `modn`.
save_table_trial_level = function(surr_type, target, modn, estimate_weights, truncation) {
  if (estimate_weights)
    weights_chr = "estwts"
  else
    weights_chr = "kwnwts"
  
  temp_tbl = all_results_tbl %>%
    filter(
      target == .env$target,
      surr_type == .env$surr_type,
      modn == .env$modn,
      estimate_weights == .env$estimate_weights,
      truncation == .env$truncation
    ) %>%
    filter(!is.na(trial))
  
  outfile_table = paste0(tables_dir, "trial_level_", surr_type, "_", target, "_", truncation, "_M", modn, "_", weights_chr, ".csv")
  write.csv(temp_tbl, outfile_table, row.names=FALSE)
}

# Helper function to save tables for the association parameter estimates
save_table_surr_measures = function(surr_type, target, modn, estimate_weights, truncation) {
  if (estimate_weights)
    weights_chr = "estwts"
  else
    weights_chr = "kwnwts"

  temp_tbl = all_results_tbl %>%
    filter(is.na(trial)) %>%
    filter(
      target == .env$target,
      surr_type == .env$surr_type,
      modn == .env$modn,
      estimate_weights == .env$estimate_weights,
      truncation == .env$truncation
    )
  outfile_table = paste0(tables_dir, "surr_measures_", surr_type, "_", target, "_", truncation, "_M", modn, "_", weights_chr, ".csv")
  write.csv(temp_tbl, outfile_table, row.names=FALSE)
}

# Save tables
plot_parameters_tbl %>%
  rowwise() %>%
  summarise(
    save_table_trial_level(surr_type, target, modn, estimate_weights, truncation),
    save_table_surr_measures(surr_type, target, modn, estimate_weights, truncation)
  )

## Tables for Paper ----------------------------------------------------

# Table with summary of results for the main text. We present the results for
# the full-population 8-trial analysis and truncated-population 6-trial analysis
# in the main text of the paper.
all_results_tbl %>%
  filter(is.na(trial)) %>%
  filter(
    target == .env$target,
  ((truncation == 1) & (modn == 2)) | ((truncation == 0) & (modn == 1))
  ) %>%
  mutate(modn = ifelse(modn == 1, "8-trial analysis", "6-trial analysis")) %>%
  select(modn, type, estimand, estimate, CI_lower, CI_upper, CI_lower_bs, CI_upper_bs, surr_type) %>%
  write.csv(file = paste0(tables_dir, "surrogacy-main-results.csv"), row.names=FALSE)

