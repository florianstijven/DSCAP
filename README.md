# DSCAP - Doubly Standardized Causal Association Parameters

This repository contains the code associated in *add ref*.

## Project Structure

The project is organized into the following directories:
* R/: Contains the R scripts and functions used for the analyses.
* data/: Directory where the data should be stored. This directory also contains an Rscript to process the original data and to generate a synthetic data set from the original data.
  Note that the original data are not included in this repository.
* results/: When running the code, the results (including tables and figures) will be saved here.

The Makefile allows one to run all Rscripts in the correct order by running `make all` in the console. 

## renv

The project uses the `renv` R package for dependency management. 
This helps to ensure that the code can be run across devices using similar environments (i.e., using the same versions of the required R packages).

# License

*to add*




