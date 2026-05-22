# Setup -------------------------------------------------------------------
# Load required packages
library(tidyverse)
library(patchwork)
library(scales)

# Specify options for saving the plots to files
figures_dir = "results/data-application/figures/"
tables_dir = "results/data-application/tables/"
raw_results_dir = "results/data-application/raw-results/"

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
  reframe(read_rds(
    file = paste0(
      raw_results_dir,
      "DSCAP_estimates_",
      ifelse(truncation, "truncated", "full"),
      "_",
      target,
      "_",
      ifelse(estimate_weights, "estwts", "kwnwts"),
      "_",
      modn,
      "_",
      surr_type
      ,
      ".rds"
    )
  ))

# Add extra variable indicating the J&J trials.
all_results_tbl = all_results_tbl %>%
  mutate(jnj_trial = ifelse(str_detect(trial, "J&J"), TRUE, FALSE))

# Read in p-values for the LRTs.
lrt_pvalues_tbl = plot_parameters_tbl %>%
  rowwise(everything()) %>%
  reframe(read_rds(
    file = paste0(
      raw_results_dir,
      "LRT_",
      ifelse(truncation, "truncated", "full"),
      "_",
      target,
      "_",
      ifelse(estimate_weights, "estwts", "kwnwts"),
      "_",
      modn,
      "_",
      surr_type
      ,
      ".rds"
    )
  ))

# Plotting -------------------------------------------------------------

# Define transformation for VE scale.
transform_VE = new_transform(
  name = "VE",
  transform = function(x) -1 * log(1 - x),
  inverse = function(x) 1 - exp(-x)
)

transform_id = new_transform(
  name = "identity",
  transform = identity,
  inverse = identity
)

# Function to make MA plots.
ma_plot = function(modn = 1,
                   estimate_weights = TRUE,
                   target,
                   truncation,
                   type,
                   range) {
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
  
  plotting_data = all_results_tbl %>%
    filter(
      type %in% c("naive", .env$type),
      modn == modn,!is.na(trial),
      estimate_weights == .env$estimate_weights,
      target == .env$target,
      truncation == .env$truncation
    )   %>%
    mutate()
  
  # Pivot the `plotting_data` to wide format with a seperate column for VE_est and
  # mean_diff_S_est. This makes it easier to plot the estimates and confidence
  # intervals in the MA plot.
  plotting_data = plotting_data %>%
    pivot_wider(
      id_cols = c(surr_type, trial, type),
      names_from = measure,
      values_from = c(estimate, CI_lower_bs, CI_upper_bs)
    ) %>%
    mutate(
      type = ifelse(type == "naive", "Unstandardized", "Standardized"),
      type = factor(type, levels = c("Unstandardized", "Standardized")),
      surr_type = ifelse(surr_type == "spike", "Binding Ab", "Neutralizing Ab"),
    )
  # If `range` is "zoomed", we zoom in on the area around the point estimates,
  # causing the CIs to not be fully contained in the plot. If `range` is "fit",
  # we show the full range of the estimates and CIs.
  transform_f = transform_id
  breaks = waiver()
  if (range == "zoomed1") {
    y_lims = c(-0.5, 1)
    x_lims = c(-0.5, 3)
  } else if (range == "fit") {
    y_lims = c(-10, 1.5)
    x_lims = c(-1, 3)
    transform_f = transform_VE
    breaks = c(-10, -2, 0, 0.5, 0.75, 0.90, 0.95, 0.99)
  } else if (range == "zoomed2") {
    y_lims = c(-0, 1)
    x_lims = c(-0.5, 3)
  }
   else {
    stop("Invalid value for `range`. Must be one of 'zoomed1', 'zoomed2', or 'fit'.")
  }
  # Plot for binding Ab.
  ma_ggplot = plotting_data %>%
    ggplot(aes(
      x = estimate_mean_diff_S_est,
      y = estimate_VE_est,
      color = trial,
      shape = trial
    )) +
    geom_errorbar(
      aes(ymin = CI_lower_bs_VE_est, ymax = CI_upper_bs_VE_est),
      width = 0.05,
      color = "darkgrey"
    ) +
    geom_errorbar(
      aes(xmin = CI_lower_bs_mean_diff_S_est, xmax = CI_upper_bs_mean_diff_S_est),
      width = 0.05,
      color = "darkgrey"
    ) +
    geom_point(size = 2, show.legend = TRUE) +
    geom_hline(yintercept = 1,
               linetype = "dashed",
               color = "lightgrey") +
    # ylim(0, 1) +
    coord_cartesian(ylim = y_lims, xlim = x_lims) +
    scale_shape_manual(
      name = "Trial",
      breaks = trials_chr,
      values = stringr::str_detect(trials_chr, "J&J") + 16
    ) +
    scale_y_continuous(transform = transform_f, breaks = breaks) +
    scale_color_manual(name = "Trial",
                       breaks = trials_chr,
                       values = viridisLite::viridis(8)) +
    xlab(bquote(tilde( ~ delta)^a)) +
    ylab(expression(tilde( ~ tau)^a ~ " (VE)")) +
    facet_grid(surr_type ~ type) +
    labs(shape = "Trial") +
    theme(legend.position = "bottom",
          legend.text = element_text(margin = margin(r = 0, unit = "pt"))) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE),
           shape = guide_legend(nrow = 2, byrow = TRUE))
  return(ma_ggplot)
}

# Save all plots.
plot_parameters_tbl %>%
  expand_grid(type = c("standardized", "doubly robust", "ipw"),
                     range = c("zoomed1", "fit", "zoomed2")) %>%
  mutate(temp = purrr::pmap(
    .l = list(
      target = target,
      modn = modn,
      estimate_weights = estimate_weights,
      truncation = truncation,
      type = type,
      range = range
    ),
    .f = function(target,
                  modn,
                  estimate_weights,
                  truncation,
                  type,
                  range) {
      if (estimate_weights)
        weights_chr = "estwts"
      else
        weights_chr = "kwnwts"
      outfile = paste0("ma_plot_",
                       target,
                       "_",
                       truncation,
                       "_",
                       "M",
                       modn,
                       "_",
                       weights_chr,
                       "_",
                       type,
                       "_",
                       range,
                       ".pdf")
      ma_plot(modn, estimate_weights, target, truncation, type, range)
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
save_table_results = function(surr_type, target, modn, estimate_weights, truncation) {
  if (estimate_weights)
    weights_chr = "estwts"
  else
    weights_chr = "kwnwts"

  all_results_tbl_filtered <- all_results_tbl %>%
    filter(
      target == .env$target,
      surr_type == .env$surr_type,
      modn == .env$modn,
      estimate_weights == .env$estimate_weights,
      truncation == .env$truncation
    ) %>%
    select(
      -target, -modn, -surr_type, -truncation, -estimate_weights
    )
  
  # Print and save to file the results for each type of estimator separately.
  sink(paste0(tables_dir, "trial_level_", surr_type, "_", target, "_", truncation, "_M", modn, "_", weights_chr, ".txt"))
  cat("**Doubly Robust Estimator**\n\n")
  all_results_tbl_filtered %>%
    filter(type == "doubly robust") %>%
    select(-type) %>%
    print()
  cat("\n\n")
  cat("**IPW Estimator**\n\n")
  all_results_tbl_filtered %>%
    filter(type == "ipw") %>%
    select(-type) %>%
    print()
  cat("\n\n")
  cat("**Standardized Estimator**\n\n")
  all_results_tbl_filtered %>%
    filter(type == "standardized") %>%
    select(-type) %>%
    print()
  cat("\n\n")
  cat("**Naive Estimator**\n\n")
  all_results_tbl_filtered %>%
    filter(type == "naive") %>%
    select(-type) %>%
    print()
  sink()
}

# Save tables
plot_parameters_tbl %>%
  rowwise() %>%
  summarise(
    save_table_results(surr_type, target, modn, estimate_weights, truncation)
  )

## Tables for Paper ----------------------------------------------------

# Table with summary of results for the main text. We present the results for
# the full-population 8-trial analysis and truncated-population 6-trial analysis
# in the main text of the paper.
sink(paste0(tables_dir, "surrogacy-main-results.txt"))
cat("**Direct Standardization Estimator --- 8-trial-analysis**\n\n")
all_results_tbl %>%
  filter(is.na(trial)) %>%
  filter(modn == 1, type == "direct standardization") %>%
  mutate(standardized = paste0(round(estimate, 2), " (", round(CI_lower_bs, 2), ", ", round(CI_upper_bs, 2), ")")) %>%
  select(measure, standardized, surr_type) %>%
  left_join(
    all_results_tbl %>%
      filter(is.na(trial)) %>%
      filter(modn == 1, type == "naive") %>%
      mutate(naive = paste0(round(estimate, 2), " (", round(CI_lower_bs, 2), ", ", round(CI_upper_bs, 2), ")")) %>%
      select(measure, naive, surr_type),
    by = c("measure", "surr_type")
  ) %>%
  select(surr_type, measure, naive, standardized)
cat("\n\n")
cat("**Direct Standardization Estimator --- 6-trial-analysis**\n\n")
all_results_tbl %>%
  filter(is.na(trial)) %>%
  filter(modn == 2, type == "direct standardization") %>%
  mutate(standardized = paste0(round(estimate, 2), " (", round(CI_lower_bs, 2), ", ", round(CI_upper_bs, 2), ")")) %>%
  select(measure, standardized, surr_type) %>%
  left_join(
    all_results_tbl %>%
      filter(is.na(trial)) %>%
      filter(modn == 2, type == "naive") %>%
      mutate(naive = paste0(round(estimate, 2), " (", round(CI_lower_bs, 2), ", ", round(CI_upper_bs, 2), ")")) %>%
      select(measure, naive, surr_type),
    by = c("measure", "surr_type")
  ) %>%
  select(surr_type, measure, naive, standardized)
sink()

# Table with summary of results for the appendix, where we present the results
# for the IPW and standardized estimators for all analyses.
sink(paste0(tables_dir, "surrogacy-appendix-results.txt"))
cat("**IPW Estimator --- 8-trial-analysis**\n\n")
all_results_tbl %>%
  filter(is.na(trial)) %>%
  filter(modn == 1, type == "ipw") %>%
  mutate(ipw = paste0(round(estimate, 2), " (", round(CI_lower_bs, 2), ", ", round(CI_upper_bs, 2), ")")) %>%
  select(measure, ipw, surr_type) %>%
  left_join(
    all_results_tbl %>%
      filter(is.na(trial)) %>%
      filter(modn == 1, type == "naive") %>%
      mutate(naive = paste0(round(estimate, 2), " (", round(CI_lower_bs, 2), ", ", round(CI_upper_bs, 2), ")")) %>%
      select(measure, naive, surr_type),
    by = c("measure", "surr_type")
  ) %>%
  select(surr_type, measure, naive, ipw)
cat("\n\n")
cat("**IPW Estimator --- 6-trial-analysis**\n\n")
all_results_tbl %>%
  filter(is.na(trial)) %>%
  filter(modn == 2, type == "ipw") %>%
  mutate(ipw = paste0(round(estimate, 2), " (", round(CI_lower_bs, 2), ", ", round(CI_upper_bs, 2), ")")) %>%
  select(measure, ipw, surr_type) %>%
  left_join(
    all_results_tbl %>%
      filter(is.na(trial)) %>%
      filter(modn == 2, type == "naive") %>%
      mutate(naive = paste0(round(estimate, 2), " (", round(CI_lower_bs, 2), ", ", round(CI_upper_bs, 2), ")")) %>%
      select(measure, naive, surr_type),
    by = c("measure", "surr_type")
  ) %>%
  select(surr_type, measure, naive, ipw)
cat("\n\n")
cat("**Doubly Robust Estimator --- 8-trial-analysis**\n\n")
all_results_tbl %>%
  filter(is.na(trial)) %>%
  filter(modn == 1, type == "doubly robust") %>%
  mutate(standardized = paste0(round(estimate, 2), " (", round(CI_lower_bs, 2), ", ", round(CI_upper_bs, 2), ")")) %>%
  select(measure, standardized, surr_type) %>%
  left_join(
    all_results_tbl %>%
      filter(is.na(trial)) %>%
      filter(modn == 1, type == "naive") %>%
      mutate(naive = paste0(round(estimate, 2), " (", round(CI_lower_bs, 2), ", ", round(CI_upper_bs, 2), ")")) %>%
      select(measure, naive, surr_type),
    by = c("measure", "surr_type")
  ) %>%
  select(surr_type, measure, naive, standardized)
cat("\n\n")
cat("**Doubly Robust Estimator --- 6-trial-analysis**\n\n")
all_results_tbl %>%
  filter(is.na(trial)) %>%
  filter(modn == 2, type == "doubly robust") %>%
  mutate(standardized = paste0(round(estimate, 2), " (", round(CI_lower_bs, 2), ", ", round(CI_upper_bs, 2), ")")) %>%
  select(measure, standardized, surr_type) %>%
  left_join(
    all_results_tbl %>%
      filter(is.na(trial)) %>%
      filter(modn == 2, type == "naive") %>%
      mutate(naive = paste0(round(estimate, 2), " (", round(CI_lower_bs, 2), ", ", round(CI_upper_bs, 2), ")")) %>%
      select(measure, naive, surr_type),
    by = c("measure", "surr_type")
  ) %>%
  select(surr_type, measure, naive, standardized)
sink()

# Table with p-values for the LRT p-values based on the asymptotic chi-squared
# distribution.
sink(paste0(tables_dir, "lrt-pvalues.txt"))
lrt_pvalues_tbl
sink()
