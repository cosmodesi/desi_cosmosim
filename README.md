# multi2pcf
Code to compute the two-point correlation function multipoles (monopole, quadrupole and hexadecapole) of different DESI mock realisations. 
The code is parallelised in OMP and optimised to be run on NERSC using the MPI wrapper [jobfork](https://github.com/cheng-zhao/jobfork). 

# UNIT 1, 2, 3 Gpc mocks - 1 box
A real- and a redshift-space (rsd) code versions are provided, with corresponding slurm scripts and input parameter files. The multipoles are computed with 200 linear bins in *0 < s < 200 Mpc/h* and 120 linear bins in *0 < μ < 1*. 

1) Create an *outputs* folder in the same directory. Here all the outputs will be stored.

2) Compile the codes as:

>real: *ifort -fopenmp xilreal.f90 -o xilreal*

>rsd: *ifort -fopenmp xilrsd.f90 -o xilrsd*

3) Launch the slurm scripts using:

>real: *sbatch xilrsd.sh*

>rsd: *sbatch xilreal.sh*

The params_rsd.inp and params_real.inp files have the following entries:
- input mock catalogue
- output DD counts: *#s(Mpc/h), μ, DD*
- output ξ(s, μ): *#s(Mpc/h), μ, ξ(s,μ)*
- output ξ<sub>l</sub>(s) (multipoles 2PCF): *#s(Mpc/h), ξ<sub>0</sub>(s), ξ<sub>2</sub>(s), ξ<sub>4</sub>(s)*
- Lbox
- redshift
- Omega_m
- Omega_l
- scale factor = 1/(1+redshift)


# EZmocks 1, 2, 3 Gpc mocks - 1000 boxes 
A real- and a redshift-space (rsd) code versions are provided, with corresponding slurm scripts and input parameter files optimised using the  [jobfork](https://github.com/cheng-zhao/jobfork) wrapper.

1) Create a *multipoles* folder in the same directory. Here all the outputs will be stored.

2) Compile the codes as:

>real: *ifort -fopenmp xilreal.f90 -o xilreal*

>rsd: *ifort -fopenmp xilrsd.f90 -o xilrsd*

3) Install jobfork in the same directory so that the slurm scripts can execute the command:

>real: *jobfork/jobfork_omp jobfork_real_list*

>rsd: *jobfork/jobfork_omp jobfork_rsd_list*


4) Launch the slurm scripts using:

>real: *sbatch xilrsd.sh*

>rsd: *sbatch xilreal.sh*


The params_xxx.inp files are named after the 1000 mock realisations. Each one of them has the following entries:
- input mock catalogue
- output DD counts: *#s(Mpc/h), μ, DD*
- output ξ(s, μ): *#s(Mpc/h), μ, ξ(s,μ)*
- output ξ<sub>l</sub>(s) (multipoles 2PCF): *#s(Mpc/h), ξ<sub>0</sub>(s), ξ<sub>2</sub>(s), ξ<sub>4</sub>(s)*
- Lbox
- redshift
- Omega_m
- Omega_l
- scale factor = 1/(1+redshift)

Scripts are provided to re-bin the pair counts in 40 linear bins in *0 < s < 200 Mpc/h* (i.e. 5 Mpc/h wide) and 120 linear bins in *0 < μ < 1*. These script read the *xilrsd_xxx.txt* and *xilreal_xxx.txt* files in the *multipoles* repository above and output the corresponding re-binned measurements in the *multipoles_5Mpc*. To run the re-binning script:

1) Create a *multipoles_5Mpc* folder in the same directory.

2) Run:

> real: python rebinned_xilreal.py
> rsd: python rebinned_xilrsd.py

Scripts are provided to compute the mean values of the multipoles and store them in the *multipoles_5Mpc* folder. Run them as:

>real: *ifort mean_rebinned_xilrsd -o meanrsd*
>        *./meanrsd*

>rsd: *ifort mean_rebinned_xilreal -o meanreal*
>       *./meanreal*

# EZmocks covariances 
The covariances o the EZmocks are computed from the multipoles in 40 linear bins in *0 < s < 200 Mpc/h* and 120 linear bins in *0 < μ < 1*, which are located in the *multipoles_EZmocks3Gpc/multipoles_5Mpc* repository. To compute the covariances, you will first need to calculate the mean values of the multipoles (see above). Then the auto- and cross-covariances between ξ<sub>0</sub>(s), ξ<sub>2</sub>(s) and ξ<sub>4</sub>(s) are obtained in 40 x 40 matrix format by running:


>real: *ifort covariances_real.f90 -o cov*
>        *./cov*

>rsd: *ifort covariances_rsd.f90 -o cov*
>       *./cov*


To assemble the 9 individual  covariances into a unique 3 x 3 matrix and output it in eBOSS format run:

> real: python covreal_eBOSSformat.py

> rsd: python covrsd_eBOSSformat.py

these scripts will output, in the *covariances* repository, the assembled (3x3) covariance and normalised covariance matrix in eBOSS format in the 

Questions, bugs, queries ---> ginevra.favole@port.ac.uk
