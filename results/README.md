This directory contains three subdirectories which contain all results of the 
analyses performed whith Rscripts located elsewhere in this repo. 

`figures/` and `tables/` contain summaries of the results in the form of figures 
and tables. 

`raw-results/` contains the so-called raw results of the standardized analyses.
The following types of files are included, where an asterisk refers to unspecified
text.

* `bootstrap_*`. Summary of the percentile bootstrap for standardized and naive
trial-level treatment effects and standardized and naive surrogacy measures.
* `cor_est_*`. Estimates of the standardized and naive surrogacy measures.
* `lrt_`. Results of the likelihood ratio test for conditional exchangeability. 
The first value is the p-value based on the asymptotic chi-squared distribution. 
The second value is the permutation p-value. The third value is the likelihood-ratio test
statistic.
* `trt_effects_*`. Standardized and naive treatment effect estimates with standard errors and
confidence intervals based on the sandwich estimator.
* `vcov_*`. The sandwich estimator of the variance-covariance matrix for all parameter estimates. Note that, because
we're using many stacked estimating equation to capture the uncertainty in all
parameters, this is a very large matrix containing superfluous information. For
instance, the estimated variance for the case-cohort sampling proportions are
included; these variance are not important.

The `*` term above is built as follows `a_b_Mc_d` where 

* `a` denotes the target trial. This is `AZ` (AstraZeneca) in all analyses.
* `b` denotes whether the `full` or `truncated` population is used for the target trial.
* `c` denotes the 8-trial analysis (`1`) or the 6-trial analysis (`2`).
* `d` denotes the handling of the case-cohort sampling weights. 




