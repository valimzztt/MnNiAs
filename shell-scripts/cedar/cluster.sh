#!/bin/bash
#SBATCH --account=def-rottler
#SBATCH --partition=default
#SBATCH --time=23:00:00 
#SBATCH --mail-user=valimzztt@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --ntasks=70             # number of MPI processes
#SBATCH --mem-per-cpu=2G      # memory; default unit is megabytes
cd MnNiAs-scf-conc
python cluster-exp-conc1.py