#!/bin/bash
#SBATCH -p regular
#SBATCH -o logs/test_HMF_CF
#SBATCH --time=120
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --constraint=haswell

module load python

# run the main nbodykit example
srun python HMF_CF.py /project/projectdirs/desi/cosmosim/proto_sim/DESI_ROCKSTAR_HALO_CATALOGS/SWIFT/1296/z1/halos_swift_Np1296_z1p000.list  test_run
