#!/bin/bash
#SBATCH --job-name=dscap
#SBATCH --ntasks=1 --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --cluster=wice
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -A lp_doctoralresearch

module use /apps/leuven/rocky8/icelake/2023b/modules/all/
module load GSL
module load CMake
module load R-bundle-CRAN/2024.06-foss-2023b

export OMP_NUM_THREADS=1

Rscript -e 'renv::rebuild("KernSmooth", type = "source")'
Rscript -e 'renv::rebuild(type = "source")'



