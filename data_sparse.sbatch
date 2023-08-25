#!/bin/bash
#SBATCH --job-name=data_sparse
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH --time=4:00:00
#SBATCH --output=data_sparse.SLURMout

module purge
module load GCC/8.3.0 OpenMPI/3.1.4 R/4.0.2

Rscript --vanilla dataGenerationSparse.R --Sample 500 --Mediator 1000
Rscript --vanilla dataGenerationSparse.R --Sample 500 --Mediator 2500
Rscript --vanilla dataGenerationSparse.R --Sample 500 --Mediator 5000