# multi2pcf
Code to compute the two-point correlation function multipoles (monopole, quadrupole and hexadecapole) of different DESI mocks. 
The code is parallelised in OMP and optimised to be run on NERSC using the MPI wrapper [jobfork](https://github.com/cheng-zhao/jobfork). 

# UNIT 2Gpc mock
A real- and a redshift-space (rsd) code version are provided, with corresponding slurm scripts and input parameter files.

1) Create an outputs folder in the same directory.

2) Compile the codes as:

>real: *ifort -fopenmp xilreal.f90 -o xilreal*

>rsd: *ifort -fopenmp xilrsd.f90 -o xilrsd*

3) Launch the slurm scripts using:

>real: *sbatch xilrsd.sh*

>rsd: *sbatch xilreal.sh*

The params_rsd.inp and params_real.inp files have the following entries:
- input mock catalogue
- output DD counts
- output xi(s, μ)
- output xi_l (multipoles 2PCF)
- Lbox
- redshift
- Omega_m
- Omega_l
- scale factor = 1/(1+redshift)

Questions, bugs, queries ---> ginevra.favole@port.ac.uk
