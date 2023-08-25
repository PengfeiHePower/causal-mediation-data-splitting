#!/bin/bash
#SBATCH --job-name=hima_comp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH --time=4:00:00
#SBATCH --output=log/hima_comp.SLURMout

module purge
module load GCC/8.3.0 OpenMPI/3.1.4 R/4.0.2

Rscript --vanilla hima.R --Filename "data/dataComp_rho1.RData"
Rscript --vanilla hima.R --Filename "data/dataComp_rho2.RData"
Rscript --vanilla hima.R --Filename "data/dataComp_rho3.RData"
Rscript --vanilla hima.R --Filename "data/dataComp_rho4.RData"
Rscript --vanilla hima.R --Filename "data/dataComp_rho5.RData"
Rscript --vanilla hima.R --Filename "data/dataComp_rho6.RData"