import numpy as np
import time
from nbodykit.lab import *
from argparse import ArgumentParser
ap = ArgumentParser(description='PropagatorCalc')
ap.add_argument('--data',type=str)
ap.add_argument('--ran',type=str)
ap.add_argument('--out',type=str)
#ap.add_argument('--pk',type=str)
ap.add_argument('--Gf', type=float,default=0.61207)
ap.add_argument('--b', type=float, default=1.4)
ap.add_argument('--Nmu',type=int,default=10)
ns = ap.parse_args()
N=512
Nmu = ns.Nmu

path = '/global/project/projectdirs/desi/users/UNIT-BAO-RSD-challenge/Reconstruction/Stage1/'
initname = '/global/cscratch1/sd/yuyu22//unitic/den0512.bin'

t0=time.time()

init_dat = np.fromfile(initname,dtype=np.float32,sep='')
arr_init = init_dat.reshape((N,N,N),order='F')
mesh_init = ArrayMesh(arr_init,BoxSize=1000)

datfile = ns.data
ranfile = ns.ran

#May have to add a 'z_rsd' column to names depending on how teams format their files.
names = ['x','y','z']
catp = CSVCatalog(datfile,names)
catran = CSVCatalog(ranfile,names)


xp = catp['x'].compute()
yp = catp['y'].compute()
zp = catp['z'].compute()
#replace zp with z_rsd when relevant
positionp = np.stack([xp,yp,zp],axis=1)

xran = catran['x'].compute()
yran = catran['y'].compute()
zran = catran['z'].compute()
positionran = np.stack([xran,yran,zran],axis=1)

catp['Position'] = positionp
catran['Position'] = positionran

meshp = catp.to_mesh(Nmesh=512,BoxSize=1000,compensated=True,interlaced=True)
meshran = catran.to_mesh(Nmesh=512,BoxSize=1000,compensated=True,interlaced=True)

array_d = meshp.compute()
array_s = meshran.compute()

arrayrecon = array_d-array_s

mesh_recon = ArrayMesh(arrayrecon,BoxSize=1000)

LOS = [0,0,1]

r_cross = FFTPower(mesh_init, mode='2d', Nmesh=512, Nmu=2*Nmu, dk=0.05, second=mesh_recon,los=LOS)
r_cross_1d = FFTPower(mesh_init, mode='1d', Nmesh=512, dk=0.05, second=mesh_recon)

r_auto_recon = FFTPower(mesh_recon,mode='2d', Nmesh=512, Nmu=2*Nmu, dk=0.05, los=LOS, poles=[0,2,4])
#r_auto_recon_1d = FFTPower(mesh_recon, mode='1d', Nmesh=512, dk=0.05)

r_auto_init = FFTPower(mesh_init,mode='2d', Nmesh=512, Nmu=2*Nmu, dk=0.05, los=LOS)
#r_auto_init_1d = FFTPower(mesh_recon, mode='1d', Nmesh=512, dk=0.05)


Gf = ns.Gf
b = ns.b
out = ns.out

pg2d=[]
pg1d=[]

pg1d.append(r_cross_1d.power['k'])
pg1d.append(r_cross_1d.power['power'].real/r_auto_init_1d.power['power'].real/(Gf*b))
np.savetxt(out+'1Dpropagator.txt',np.column_stack([pg1d[0],pg1d[1]]),header='Gf= %lf \nb + %lf \n dk=0.05 \nkmean \t \t C(k)' % (Gf,b))


pg2d.append(r_cross.power[:,Nmu]['k'])
pg2d.append(r_cross.power[:,Nmu]['power'].real/r_auto_init.power[:,Nmu]['power'].real/(Gf*b))
pg2d.append(r_cross.power[:,2*Nmu-1]['power'].real/r_auto_init.power[:,2*Nmu-1]['power'].real/(Gf*b))
np.savetxt(out+'2Dpropagator.txt',np.column_stack([pg2d[0],pg2d[1],pg2d[2]]),header='Gf = %lf \nb = %lf \ndk=0.05 \nkmean \t \t C(k,mu=%lf) \t \t C(k,mu=%lf)' %(Gf,b,r_cross.power.coords['mu'][Nmu],r_cross.power.coords['mu'][2*Nmu-1]))


reconell=[]
poles = r_auto_recon.poles
reconell.append(poles['k'])
for ell in [0,2,4]:
    P = poles['power_%d' %ell].real
    if ell==0:
        P = P-poles.attrs['shotnoise']
    reconell.append(P)

np.savetxt(out+'reconstructed_multipoles.txt',np.column_stack([reconell[0],reconell[1],reconell[2],reconell[3]]),header='dk=0.05 \nshotnoise= %lf \nkmean \t \t \t P_0-shot \t \t \t P_2 \t \t \t P_4' %(poles.attrs['shotnoise']))


print('%lf seconds' %(time.time()-t0))
