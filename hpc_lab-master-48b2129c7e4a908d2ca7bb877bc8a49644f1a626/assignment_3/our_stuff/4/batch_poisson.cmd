#!/bin/bash
#SBATCH -D /home/hpc/t1221/di34mip/
#SBATCH -o ./logs/a3_4_cg/%j_%N.txt
#SBATCH --clusters=mpp3
#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=manuel.rothenberg@tum.de
mpiexec -n 6 ./code/a3_4_cg/poisson 0.0625 1000 0.0001