#!/bin/bash
##SBATCH -p action -A action
#SBATCH -p debug
#SBATCH -J GET_INITCONDS
#SBATCH -o output
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
##SBATCH --mem 10GB

module load gaussian
module load openbabel

g16 < geometry.com > geometry.out
babel -i out geometry.out -o com geometry.xyz

