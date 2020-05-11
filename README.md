# Requirements

* [Python](http://www.python.org) 3.5 or 3.6

* [Astropy](http://www.astropy.org)

* [configparser](https://docs.python.org/3/library/configparser.html)

* [CAMB](https://camb.readthedocs.io/en/latest/)

# gallightcone

This is a simple code to generate a lightcone from Shadab's HOD output.

To run this code use:

```python 
python build_gallightcone_multibox.py configfile X
```
configfile will look like (see test.config for an example):
```
[dir]
dir_out     Output directory
dir_gcat    Directory where the input out_*p.list.gcat files are located
file_alist  File containing list of simulation scale factors
file_camb   CAMB configuration file

[sim]
boxL        Length of the simulation box [Mpc/h]
shellwidth  Width of shell [Mpc/h]

[sample]
galtype     Type of galaxy used - 1=LRG,2=QSO,3=ELG 
```

X is the shell number (0 being z=0)
e.g. if you define shellwidth = 25 in the config file and choose X=1 then you are selecting the shell between 25-50 Mpc/h 

The code will output galaxy catalogs in each shell, which you will have to combine.

If you have questions or comments please email yomori@stanford.edu
