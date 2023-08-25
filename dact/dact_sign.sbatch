#!/bin/bash
#SBATCH --job-name=dact_sign
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH --time=4:00:00
#SBATCH --output=log/dact_sign.SLURMout

module purge
module load GCC/8.3.0 OpenMPI/3.1.4 R/4.0.2

Rscript --vanilla dact.R --Filename "data/dataSign1.RData"
Rscript --vanilla dact.R --Filename "data/dataSign3.RData"
Rscript --vanilla dact.R --Filename "data/dataSign5.RData"
Rscript --vanilla dact.R --Filename "data/dataSign7.RData"
Rscript --vanilla dact.R --Filename "data/dataSign9.RData"
Rscript --vanilla dact.R --Filename "data/dataSign13.RData"