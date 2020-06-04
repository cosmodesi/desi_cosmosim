#! /bin/bash

#SBATCH -J mocks_for_covariances
#SBATCH -N 59
#SBATCH -C haswell
#SBATCH -A desi
#SBATCH -o %x.o%j
#SBATCH -q debug
#SBATCH -t 30:00
#SBATCH --mail-user=krolewski@berkeley.edu
#SBATCH --mail-type=ALL

module unload python/2.7-anaconda-4.4

source /global/common/software/m3035/conda-activate.sh 3.7

bcast-pip https://github.com/bccp/nbodykit/archive/master.zip

#srun -N 16 -n 128 python fourier_space_for_covariances.py 0 64
srun -N 59 -n 466 python fourier_space_for_covariances.py 768 1000
#srun -N 64 -n 512 python fourier_space_for_covariances.py 256 512
#srun -N 64 -n 512 python fourier_space_for_covariances.py 512 768
