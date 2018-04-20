#!/bin/bash
#SBATCH -o /home/hpc/t1221/lu26xum/a2_2/%j_%N.txt
#SBATCH -D /home/hpc/t1221/lu26xum/a2_2
#SBATCH --clusters=mpp3
#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=gw.poerwawinata@tum.de
export OMP_NUM_THREADS=45
mpiexec ./phicalculation_reduction
