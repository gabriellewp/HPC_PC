#!/bin/bash
#SBATCH -o /home/hpc/t1221/lu26xum/%j_%N.txt
#SBATCH -D /home/hpc/t1221/lu26xum
#SBATCH --nodes=3
#SBATCH --mail-user=gw.poerwawinata@tum.de
srun mpirun -np 8 IMB-MPI1 SendRecv Reduce
