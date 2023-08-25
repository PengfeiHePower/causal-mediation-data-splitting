#!/bin/bash
#SBATCH --job-name=hima_sparse
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH --time=4:00:00
#SBATCH --output=log/hima_sparse.SLURMout

module purge
module load GCC/8.3.0 OpenMPI/3.1.4 R/4.0.2

Rscript --vanilla hima.R --Filename "data/data5001000.RData"
Rscript --vanilla hima.R --Filename "data/data5002500.RData"
Rscript --vanilla hima.R --Filename "data/data5005000.RData"