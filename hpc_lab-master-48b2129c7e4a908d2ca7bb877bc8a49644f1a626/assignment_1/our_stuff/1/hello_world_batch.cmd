#!/bin/bash
#SBATCH -o /home/hpc/t1221/di34mip/logs/a1_1_hello_world/%j_%N.txt
#SBATCH -D /home/hpc/t1221/di34mip/code/a1_1_hello_world
#SBATCH --clusters=mpp3
#SBATCH --nodes=2-2
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=manuel.rothenberg@tum.de
mpiexec ./hello_world.out