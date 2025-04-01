#!/bin/bash
##SBATCH -p action -A action
#SBATCH -p action
#SBATCH -J GET_INITCONDS
#SBATCH -o output
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem 40GB

python3 make_InitConds.py Fulvene.xyz


