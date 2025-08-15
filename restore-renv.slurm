#!/bin/bash
#SBATCH --job-name=dscap
#SBATCH --ntasks=1 --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --cluster=wice
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -A lp_doctoralresearch

module use /apps/leuven/rocky8/icelake/2024a/modules/all/

module load R

export OMP_NUM_THREADS=1

Rscript -e "renv::restore()"




