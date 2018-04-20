#!/bin/bash
#SBATCH -D /home/hpc/t1221/di34mip/
#SBATCH -o ./logs/a2_4_dgemm/%j_%N.txt
#SBATCH --clusters=mpp3
#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=manuel.rothenberg@tum.de
export OMP_NUM_THREADS=64
mpiexec ./code/a2_4_dgemm/dgemm 1024 1 8 1
