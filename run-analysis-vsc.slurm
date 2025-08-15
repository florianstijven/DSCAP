#!/bin/bash
#SBATCH --job-name=dscap
#SBATCH --ntasks=1 --cpus-per-task=36
#SBATCH --time=6:00:00
#SBATCH --cluster=wice
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -A lp_doctoralresearch

export OMP_NUM_THREADS=1

Rscript -e "renv::restore()"

make all -j 8



