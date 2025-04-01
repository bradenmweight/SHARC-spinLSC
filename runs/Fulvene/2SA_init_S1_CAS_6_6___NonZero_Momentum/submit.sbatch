#!/bin/bash
#SBATCH -p standard
#SBATCH -J FUL_sLSC
#SBATCH -o output.slurm
#SBATCH -t 5-00:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
##SBATCH --mem 8GB

# module load gcc
# module load openblas
# module load blas
# module load openmpi
# module load impi/2017.5

sh run.sh
