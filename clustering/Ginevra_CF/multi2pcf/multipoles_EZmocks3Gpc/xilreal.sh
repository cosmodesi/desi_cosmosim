#!/bin/bash
#SBATCH --qos=shared
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --constraint=haswell

#export OMP_PROC_BIND=true
#export OMP_PLACES=threads
export OMP_NUM_THREADS=32

ifort -fopenmp xilreal.f90 -o xilreal
jobfork/jobfork_omp jobfork_real_list
