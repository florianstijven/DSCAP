#!/bin/bash
#
#SBATCH --job-name=restore-renv
#SBATCH --ntasks=1 --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --cluster=wice
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -A lp_doctoralresearch

export OMP_NUM_THREADS=1

module load cluster/wice/batch_sapphirerapids
module load R/4.4.2-gfbf-2024a
module load R-bundle-CRAN/2024.11-foss-2024a

Rscript -e "renv::status()" -e "renv::restore()" - "renv::status()"