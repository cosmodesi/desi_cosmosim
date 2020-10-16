#!/bin/bash
#SBATCH -p debug
#SBATCH -o nbkit-example
#SBATCH --time=5
#SBATCH --nodes=2
#SBATCH --tasks-per-node=1
#SBATCH --constraint=haswell
# load nbodykit

module purge
module load python h5py-parallel

source /global/common/software/m3035/conda-activate.sh 3.7

# run the main nbodykit example
srun python power_spectrum_cosmosim_folding.py /project/projectdirs/desi/cosmosim/proto_sim/IC/Panphasia/Examples/192/DESI_IC_192.%d.hdf5 1 256 1 test.dat
