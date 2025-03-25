# Rscripts

This directory contains the following Rscripts:

* `data-exploration.R`. This Rscript generates figures and tables that are saved
into subdirectories of `results/`. This Rscript assumes that the original data
`processed_data.csv` are present in the `data/` directory. Hence, this script
cannot be run without the original data.
* `estimate_dscap.R`. This Rscript performs the main analysis. This Rscript saves 
  the results in the `results/raw-results/` directory. It takes the 
  following command-line arguments (in the presented order), which parameterize the type of analsysis that is performed.
  - Antibody type: `"spike"` for binding Ab and `"neut"` for neutralizing Ab. 
  - Target trial: This should be `"AZ"` for AstraZeneca as target trial.
  - GLM formula: The right-hand side of the formula for the logistic and linear regression models that 
    are used for the standardization estimators.
  - Selected trials: `1` for the 8-trial analysis, `2` for the 6-trial analysis 
    that excludes J&J (Brazil) and J&J (Colombia).
  - Number of bootstrap replications.
  - Number of replications in the permutation test.
  - Treat the case-cohort sampling weights as known or estimated. `1` treats them as estimated, 
    `2` treats them as known.
  - Truncated target population. `1` for a truncated target population where only subjects
  from the target trial are retained whose risk score lies in the intersection of the supports for
  the risk score across the selected trials. `0` for no truncation (i.e., using the full population).
  - Location of the data. Should be `"data/processed_data_synthetic.csv` for the 
  analysis with the synthetic data and `"data/processed_data.csv"` for the analysis 
  with the original data (which cannot be made available).
* `plots_tables.R`. This Rscript takes the "raw results" saved in `results/raw-results/`
after running `estimate_dscap.R` and summarizes the results in figures and
tables.

# Helper Functions

The `helper-functions/` directory contains helper functions which are used in
the above Rscripts. 

* `helper-functions/estFUN.R` contains the `estFUN_taudelta()` function which implements
the estimating function underlying the DSCAP methodology.
* `helper-functions/processing_results.R` contains helper functions for further processing
the results saved in `results/raw-results/`.
* `helper-functions/runDSCAP.R` contains the `RunDSCAP()` function that implements
the standardized estimators for the treatment effects and surrogacy measures. 




