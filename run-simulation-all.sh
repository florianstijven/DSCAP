#!/bin/bash
#SBATCH --output=par-%J.out
#SBATCH --ntasks=1 --cpus-per-task=36 --nodes=1
#SBATCH --time=20:00:00
#SBATCH --cluster=wice
#SBATCH --mail-type=END,FAIL,REQUEUE,STAGE_OUT
#SBATCH -A lp_doctoralresearch

module load cluster/wice/batch
module load R/4.4.2-gfbf-2024a
module load R-bundle-CRAN/2024.11-foss-2024a

export OMP_NUM_THREADS=1

Rscript -e "renv::restore()"

make simulation