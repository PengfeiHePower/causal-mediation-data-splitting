#!/bin/bash
#SBATCH --job-name=hima_toep10
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH --time=4:00:00
#SBATCH --output=log/hima_toep10.SLURMout

module purge
module load GCC/8.3.0 OpenMPI/3.1.4 R/4.0.2

Rscript --vanilla hima.R --Filename "data/dataToep_r1_b10.RData"
Rscript --vanilla hima.R --Filename "data/dataToep_r2_b10.RData"
Rscript --vanilla hima.R --Filename "data/dataToep_r3_b10.RData"
Rscript --vanilla hima.R --Filename "data/dataToep_r4_b10.RData"
Rscript --vanilla hima.R --Filename "data/dataToep_r5_b10.RData"
Rscript --vanilla hima.R --Filename "data/dataToep_r6_b10.RData"