#!/bin/bash
#SBATCH -o /home/hpc/t1221/di34mip/logs/a1_4_dgemm/%j_%N.txt
#SBATCH -D /home/hpc/t1221/di34mip/code/a1_4_dgemm
#SBATCH --clusters=mpp3
#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=manuel.rothenberg@tum.de
export OMP_NUM_THREADS=4
mpiexec ./dgemm
