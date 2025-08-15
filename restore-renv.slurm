#!/bin/bash
#SBATCH --job-name=dscap-restore
#SBATCH --ntasks=1 --cpus-per-task=2
#SBATCH --time=1:00:00
#SBATCH --cluster=wice
#SBATCH --partition=batch
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -A lp_doctoralresearch

module load cluster/wice/batch
module load R

export OMP_NUM_THREADS=1

Rscript -e 'renv::install("lme4@1.1-36", type = "source")'
Rscript -e 'options(pkgType = "source"); renv::restore()'



