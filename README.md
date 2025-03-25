# DSCAP - Doubly Standardized Causal Association Parameters

This repository contains the code used for the analyses in *add ref*.

## Project Structure

The project is organized into the following directories:
* R/: Contains the R scripts and functions used for the analyses.
* data/: Directory where the data should be stored. This directory also contains an Rscript to process the original data and to generate a synthetic data set from the original data.
  Note that the original data are not included in this repository.
* results/: When running the code, the results (including tables and figures) will be saved here.

The Makefile allows one to run all Rscripts in the correct order by running `make all` in the console. 

## Reproducibility

The results with the synthetic data can be reproduced by executing the Makefile
(by running `make all` in the command line). This will run all statistical
analyses (which produces results saved to `results/raw-results/`) and will
produce all plots and tables presented in the paper (which are saved to
`results/figures/` and `results/tables/`. Note that this requires the presence
of `data/processed_data_synthetic.csv`. Reproducibility is also facilitated
through the use of the `renv` R package.

Since the analyses are computationally intensive, we ran the analysis on a
computing cluster, for which we used the `run-analysis.sh` batch script. If
these analyses are run on a personal computer, there should be enough RAM
available because the sandwich estimator based on the `geex` R package is greedy
in its use of RAM. Running the Makefile is expected to take 10 hours on a recent
personal computer.

## renv

The project uses the `renv` R package for dependency management. This helps to
ensure that the code can be run across devices using similar environments (i.e.,
using the same versions of the required R packages).

# License

*to add*




