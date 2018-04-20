#!/bin/bash
#SBATCH -o /home/hpc/t1221/lu26xum/%j_%N.txt
#SBATCH -D /home/hpc/t1221/lu26xum
#SBATCH --nodes=5
#SBATCH --mail-user=gw.poerwawinata@tum.de
srun mpiexec -n 4 ./mpibcast2 4
