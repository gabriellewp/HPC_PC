#!/bin/bash
#SBATCH -o /home/hpc/t1221/di34mip/logs/a1_3_gauss/%j_%N.txt
#SBATCH -D /home/hpc/t1221/di34mip/code/a1_3_gauss
#SBATCH --clusters=mpp3
#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=manuel.rothenberg@tum.de
mpiexec ./gauss