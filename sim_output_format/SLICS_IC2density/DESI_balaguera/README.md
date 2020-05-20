# DESI


This c++ code reads binary files with SLICS IC and intepolate to a mesh.
balaguera@iac.es 2020

Compile with make file:
make

Execute with flag -h for help.
Executable get_dens_field.exe


Execute with flag -r with parameters as explained below:

Input parameter are:

1 -> <Path>: Path to the IC folder

2 -> <nIC>: Number x of the SLAC realization.

3 -> <MAS>: Interpolation scheme: 0 (NGP), 1 (CIC), 2 (TSC), 3 (PCS)

4 -> <Nres>: Number of cells per dimention (e.g. 256 for a 256³ mesh)

5 -> <Outputdir>: Name of output directory


The code reads the 64 files located in the path <Path>/LOS<nIC>/xv##.ic
where ## is from 0 to 63 (64 sub-volumes used in the generation of IC)
The output is written in a single-precision binary file (C order)
at Outputdir/SLICS_IC_LOS<nIC>_Nres<Nres>_MAS<MAS>.dat

Example:

./get_dens_field.exe -r SLICS_IC 1025 1 256 Output

This reads the 64 files SLICS_IC/LOS1025/xv*.ic 
and geneates a file Outputdir/SLICS_IC_LOS1025_Nres256_MAS0.dat

