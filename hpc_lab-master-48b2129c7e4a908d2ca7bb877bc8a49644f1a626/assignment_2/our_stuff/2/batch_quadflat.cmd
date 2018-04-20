#!/bin/bash
#SBATCH -o /home/hpc/t1221/lu26xum/a2_2/%j_%N.txt
#SBATCH -D /home/hpc/t1221/lu26xum/a2_2
#SBATCH --constraint=snc4,flat
#SBATCH --clusters=mpp3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=gw.poerwawinata@tum.de
srun numactl --physcpubind=1 ./stream.omp.AVX2.72M.10x.icc
