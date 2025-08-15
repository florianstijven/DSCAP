#!/bin/bash
#SBATCH --job-name=dscap-restore
#SBATCH --ntasks=1 --cpus-per-task=2
#SBATCH --time=1:00:00
#SBATCH --cluster=genius
#SBATCH --partition=batch
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -A lp_doctoralresearch

module load R

export OMP_NUM_THREADS=1

Rscript -e 'renv::rebuild("mvtnorm", type = "source")'
# Rscript -e 'options(pkgType = "source"); renv::restore()'



