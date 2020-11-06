"""
HMF_CF.py

Calculate the Halo mass function and Correlation function from Halo catalogs
created by Rockstar. Do this for several different low mass cuts

Takes the file path as a command line argument
"""
import numpy as np
from Corrfunc.theory.xi import xi
import sys

import multiprocessing


def read_Rockstar(path):
    Mh, x, y, z, PID = np.loadtxt('%s'%path,usecols=(2,8,9,10,41),unpack=True)
    
    mask = np.where((PID == -1))
    
    Mh = Mh[mask]
    x = x[mask]
    y = y[mask]
    z = z[mask]
    
    return x,y,z,Mh

def HMF(Mhalo):
    M_min = 8.0
    M_max = 15.0
    dlogM = 0.2
    M_b = np.arange(M_min,M_max,dlogM)
    M_b1 = M_b + dlogM*0.5
    nbin = len(M_b)
    V_s = (500)**3 # (Mpc/h)**3
    
    Mlh = np.log10(Mhalo[Mhalo>0])
    
    Nh, bins_edges = np.histogram(Mlh,bins=np.append(M_b,M_max))
    dHMF=Nh/(V_s*dlogM)
    
    return M_b1, dHMF

def CF(Mhalo,x,y,z):
    
    Mlh = np.log10(Mhalo[Mhalo>0])
    
    nthreads = multiprocessing.cpu_count()
    box = 500.
    nbins = 20
    bins = 10**np.linspace(np.log10(0.5),np.log10(50),nbins+1)
    
    #-------------------------------------------------------------
    # Correlation function samples: logM > 10.5, 11.0, 11.5
    #-------------------------------------------------------------

    mh1 = np.where(Mlh > 10.5)
    mh2 = np.where(Mlh > 11.0)
    mh3 = np.where(Mlh > 11.5)

    print('Number of haloes Mlh > 10.5', len(x[mh1]))
    print('Number of haloes Mlh > 11.0', len(x[mh2]))
    print('Number of haloes Mlh > 11.5', len(x[mh3]))

    xi_h1 = xi(box, nthreads, bins, x[mh1], y[mh1], z[mh1], output_ravg=True)
    print("xi1 done")
    xi_h2 = xi(box, nthreads, bins, x[mh2], y[mh2], z[mh2], output_ravg=True)
    print("xi2 done")
    xi_h3 = xi(box, nthreads, bins, x[mh3], y[mh3], z[mh3], output_ravg=True)
    print("xi3 done")


    xir1 = xi_h1['xi']
    xir2 = xi_h2['xi']
    xir3 = xi_h3['xi']

    rr1 = xi_h1['ravg']
    rr2 = xi_h2['ravg']
    rr3 = xi_h3['ravg']
    
    return rr1, xir1, rr2, xir2, rr3, xir3


path = sys.argv[1]

pathout = sys.argv[2]

print("Reading in data")
x, y, z, Mh = read_Rockstar(path)
print("Measuring HMF")
M_b1, dHMF = HMF(Mh)
print("Saving HMF")
HMF_out = np.zeros((len(M_b1),2))
HMF_out[:,0] = M_b1
HMF_out[:,1] = dHMF

np.savetxt(pathout + "_HMF.txt",HMF_out)
print("Measuring CF")
rr1, xir1, rr2, xir2, rr3, xir3 = CF(Mh,x,y,z)

CF_out = np.zeros((len(rr1),2))
CF_out[:,0] = rr1
CF_out[:,1] = xir1

np.savetxt(pathout + "_CF1.txt",CF_out)

CF_out[:,0] = rr2
CF_out[:,1] = xir2

np.savetxt(pathout + "_CF2.txt",CF_out)


CF_out[:,0] = rr3
CF_out[:,1] = xir3

np.savetxt(pathout + "_CF3.txt",CF_out)

print("Done!")


