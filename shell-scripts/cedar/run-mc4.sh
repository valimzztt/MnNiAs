#!/bin/bash
#SBATCH --account=def-rottler
#SBATCH --partition=default
#SBATCH --time=23:00:00 
#SBATCH --mail-user=valimzztt@gmail.com
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=2
#SBATCH --mem-per-cpu=10G      # memory; default unit is megabytes

cd MnNiAs-smol
module load mpi
mpiexec -n 8 python run-mc-smol2.py