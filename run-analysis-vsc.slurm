#!/bin/bash
#SBATCH --job-name=dscap
#SBATCH --ntasks=1 --cpus-per-task=36
#SBATCH --time=10:00:00
#SBATCH --cluster=wice
#SBATCH --partition=batch
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -A lp_doctoralresearch

module load cluster/wice/batch
module load NLopt
module load R/4.4.2-gfbf-2024a

export OMP_NUM_THREADS=1

Rscript -e "renv::restore()"

make all -j 8



