#!/bin/bash
#SBATCH --qos=shared
#SBATCH --time=24:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --constraint=haswell

#export OMP_PROC_BIND=true
#export OMP_PLACES=threads
export OMP_NUM_THREADS=16

srun ./xilrsd < params_rsd.inp
