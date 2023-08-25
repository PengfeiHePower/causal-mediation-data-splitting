#!/bin/bash
#SBATCH --job-name=dact_toep20
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH --time=4:00:00
#SBATCH --output=log/dact_toep20.SLURMout

module purge
module load GCC/8.3.0 OpenMPI/3.1.4 R/4.0.2

Rscript --vanilla dact.R --Filename "data/dataToep_r1_b20.RData"
Rscript --vanilla dact.R --Filename "data/dataToep_r2_b20.RData"
Rscript --vanilla dact.R --Filename "data/dataToep_r3_b20.RData"
Rscript --vanilla dact.R --Filename "data/dataToep_r4_b20.RData"
Rscript --vanilla dact.R --Filename "data/dataToep_r5_b20.RData"
Rscript --vanilla dact.R --Filename "data/dataToep_r6_b20.RData"