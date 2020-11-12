import numpy as np
import time
from nbodykit.lab import *
from argparse import ArgumentParser
import pandas as pd
ap = ArgumentParser(description='PropagatorCalc')
ap.add_argument('--data',type=str)
ap.add_argument('--ran',type=str)
ap.add_argument('--out',type=str)
#ap.add_argument('--pk',type=str)
ap.add_argument('--Gf', type=float,default=0.61207)
ap.add_argument('--b', type=float, default=1.4)
ap.add_argument('--Nmu',type=int,default=10)
ap.add_argument('--Nmesh',type=int,default=512)
ns = ap.parse_args()
Nmu = ns.Nmu
Nmesh=ns.Nmesh
N=Nmesh
path = '/global/project/projectdirs/desi/users/UNIT-BAO-RSD-challenge/Reconstruction/Stage1/'
initname = '/global/cscratch1/sd/yuyu22//unitic/den%04d.bin'%Nmesh
box_size=1000
t0=time.time()
print(f"Loading initial density field {initname}", flush=True)
init_dat = np.fromfile(initname,dtype=np.float32,sep='')
print("Done")
arr_init = init_dat.reshape((N,N,N),order='F')
print("Generating mesh")
mesh_init = ArrayMesh(arr_init,BoxSize=box_size)
print("Done")
#datfile = path+ns.data
#ranfile = path+ns.ran
datfile = ns.data
ranfile = ns.ran

#May have to add a 'z_rsd' column to names depending on how teams format their files.
names = ['x','y','z']
print(f"Loading {datfile}", flush=True)
try:
    catp = CSVCatalog(datfile,names=names)
except:
    catdf = np.genfromtxt(datfile, names=names, usecols=(0, 1, 2))
    print("Loaded with numpy successfully")

    catp = ArrayCatalog(catdf)

print(f"Loading {ranfile}", flush=True)
try:
    catran = CSVCatalog(ranfile,names=names)
except:
    randf = np.genfromtxt(ranfile, names=names, usecols=(0, 1, 2))
    print("Loaded with numpy successfully")

    catran = ArrayCatalog(randf)

print("Done")

xp = catp['x'].compute()
yp = catp['y'].compute()
zp = catp['z'].compute()
#replace zp with z_rsd when relevant
positionp = np.stack([xp,yp,zp],axis=1)

n_elem = xp.shape[0]


xran = catran['x'].compute()
yran = catran['y'].compute()
zran = catran['z'].compute()
positionran = np.stack([xran,yran,zran],axis=1)

n_rand = xran.shape[0]

catp['Position'] = positionp
catran['Position'] = positionran
compensated=True
interlaced=True
meshp = catp.to_mesh(Nmesh=Nmesh,BoxSize=box_size,compensated=compensated,interlaced=interlaced)
meshran = catran.to_mesh(Nmesh=Nmesh,BoxSize=box_size,compensated=compensated,interlaced=interlaced)

array_d = meshp.compute()
array_s = meshran.compute()

arrayrecon = array_d-array_s

mesh_recon = ArrayMesh(arrayrecon,BoxSize=box_size)

LOS = [0,0,1]
kmin=0.005 #leftmost edge of first bin
dk=0.005
kmax=1 #rightmost edge of last bin
kbins = np.arange(kmin, kmax, dk)
kbincenters = 0.5 * (kbins[1:] + kbins[:-1])

r_cross = FFTPower(mesh_init, mode='2d', Nmesh=Nmesh, Nmu=2*Nmu, dk=dk, second=mesh_recon,los=LOS, kmin=kmin, kmax=kmax)
r_cross_1d = FFTPower(mesh_init, mode='1d', Nmesh=Nmesh, dk=dk, second=mesh_recon, kmin=kmin, kmax=kmax)

r_auto_recon = FFTPower(mesh_recon,mode='2d', Nmesh=Nmesh, Nmu=2*Nmu, dk=dk, los=LOS, poles=[0,2,4], kmin=kmin, kmax=kmax)
r_auto_recon_1d = FFTPower(mesh_recon, mode='1d', Nmesh=Nmesh, dk=dk, kmin=kmin, kmax=kmax)

r_auto_init = FFTPower(mesh_init,mode='2d', Nmesh=Nmesh, Nmu=2*Nmu, dk=dk, los=LOS, kmin=kmin, kmax=kmax)
r_auto_init_1d = FFTPower(mesh_recon, mode='1d', Nmesh=Nmesh, dk=dk, kmin=kmin, kmax=kmax)

Gf = ns.Gf
b = ns.b
out = ns.out
pg2d=[]

pg1d=[]
pg1d.append(r_cross_1d.power['k'])
pg1d.append(r_cross_1d.power['power'].real/r_auto_init_1d.power['power'].real/(Gf*b))
np.savetxt(out+'1Dpropagator.txt',np.column_stack([pg1d[0],pg1d[1]]),header='Gf= %lf \nb + %lf \n dk=dk#\nkmean \t \t C(k)' % (Gf,b))


pg2d.append(r_cross.power[:,Nmu]['k'])
pg2d.append(r_cross.power[:,Nmu]['power'].real/r_auto_init.power[:,Nmu]['power'].real/(Gf*b))
pg2d.append(r_cross.power[:,2*Nmu-1]['power'].real/r_auto_init.power[:,2*Nmu-1]['power'].real/(Gf*b))
np.savetxt(out+'2Dpropagator.txt',np.column_stack([pg2d[0],pg2d[1],pg2d[2]]),header='Gf = %lf \nb = %lf \ndk=dk \nkmean \t \t C(k,mu=0.05) \t \t C(k,mu=0.95' %(Gf,b))

reconell=[]
poles = r_auto_recon.poles
#reconell.append(poles['k'])
reconell.append(kbincenters)
for ell in [0,2,4]:
    P = poles['power_%d' %ell].real
    if ell==0:
        if poles.attrs['shotnoise']==0:
            print(f"Nbodykit computed shot noise 0. Correcting with 1/nbar")
            poles.attrs['shotnoise'] = (box_size**3 / n_elem) + (box_size**3 / n_rand)
            print(f"New shot noise = {poles.attrs['shotnoise']}")
        P = P-poles.attrs['shotnoise']
    reconell.append(P)

np.savetxt(out+'reconstructed_multipoles.txt',np.column_stack([reconell[0],reconell[1],reconell[2],reconell[3]]),header='dk=%lf \nshotnoise= %lf \nkcenter \t \t \t P_0-shot \t \t \t P_2 \t \t \t P_4' %(dk, poles.attrs['shotnoise']))


print('%lf seconds' %(time.time()-t0))
