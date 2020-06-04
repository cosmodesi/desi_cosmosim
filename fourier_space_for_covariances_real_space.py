from nbodykit.lab import *
import numpy as np
import time
import sys
import glob
from mpi4py import MPI

# Input arguments: file paths and names, maximum R, number of bins, box size
# Number of bins is actually the number of *bins*, and not edges, i.e.
# if you want 40 bins spaced at 5 h^-1 Mpc between 0 and 200, specify Nbins = 40
#input_path = sys.argv[1]

#Nmesh = int(sys.argv[3])
with TaskManager(cpus_per_task=2, use_all_cpus=True) as tm:
	begin_ind = int(sys.argv[1])
	finish_ind = int(sys.argv[2])

	Nmesh = 512
	BoxSize = 3000
	#BoxSize = int(sys.argv[4])

	binning = 'linear'

	t0 = time.time()

	# Output name
	output_path = '/project/projectdirs/desi/users/akrolew/mock_challenge/UNIT/'

	#Ncores_per_task = 4

	#for i in range(1000):
	#	if rank//Ncores_per_task == i:
	input_path = '/global/project/projectdirs/desi/cosmosim/UNIT-BAO-RSD-challenge/3Gpc/ELG-single-tracer/EZmock/'
	input_files = glob.glob(input_path + '*.dat')[begin_ind:finish_ind]
	
	for input_file in tm.iterate(input_files):
		rank = input_file.split('/')[-1].split('-')[-1].split('.dat')[0]
		print('input file',input_file)
		#output_subdir = 'FS_Krolewski_MockChallenge_1/'
		output_name_pk = 'EZmock_Nmesh512_pk2D_kmu_real_space/%i.txt' % int(rank)
		output_name_pk_l = 'EZmock_Nmesh512_pkl_real_space/%i.txt' % int(rank)

		# Read file
		#data = np.loadtxt(input_path + input_file)

		# Create catalog
		#cat = ArrayCatalog({'RSDPosition': np.concatenate((data[:,:2],data[:,3].reshape(len(data),1)),axis=1)})
		#cat.attrs['BoxSize'] = BoxSize
		names = ['x','y','z','vx','vy','vz']
		cat = CSVCatalog(input_file, names)
		redshift = 0.9873
		Omega_m = 0.3089
		H = 100 * np.sqrt(Omega_m * (1+redshift)**3 + 1-Omega_m)
		zrsd = cat['z'][:,None]
		x = cat['x'][:,None]
		y = cat['y'][:,None]

		x = np.array(x)
		y = np.array(y)
		zrsd = np.array(zrsd)

		cond = ~np.isnan(zrsd)
		zrsd = zrsd[cond]

		x = x[cond]
		y = y[cond]

		x = x.reshape((len(x),1))
		y = y.reshape((len(y),1))
		zrsd = zrsd.reshape((len(zrsd),1))

		rsd_position =x * [1, 0, 0] + y * [0, 1, 0] + zrsd * [0, 0, 1]

		cat = ArrayCatalog({'RSDPosition': rsd_position}, comm=cat.comm)

		cat.attrs['BoxSize'] = BoxSize

		# Create mesh (real space)
		mesh = cat.to_mesh(window='tsc',Nmesh=Nmesh,compensated=True,interlaced=True,position='RSDPosition')
		# Need mesh to be big enough so that nyquist frequency is > kmax
		# pi * Nmesh/Lbox = 3.22 h Mpc^{-1} for Nmesh = 1024 and Lbox = 1000 h^{-1} Mpc, 
		# vs kmax = 2.262 h Mpc^{-1}
		# Interlacing + TSC gives very unbiased results up to the Nyquist frequency:
		# https://nbodykit.readthedocs.io/en/latest/cookbook/interlacing.html
		k_Nyquist = np.pi * Nmesh/BoxSize
		
		'''# Save 2D Power spectrum
		kx = np.real(np.fromfile('/global/cfs/cdirs/desi/users/akrolew/mock_challenge/UNIT/kx.bin',dtype=np.complex64))
		ky = np.real(np.fromfile('/global/cfs/cdirs/desi/users/akrolew/mock_challenge/UNIT/ky.bin',dtype=np.complex64))
		kz = np.real(np.fromfile('/global/cfs/cdirs/desi/users/akrolew/mock_challenge/UNIT/kz.bin',dtype=np.complex64))
		
		c1 = mesh.compute(Nmesh=Nmesh,mode='complex')
		r1 = c1.preview(mesh.attrs['Nmesh'])
		c1 = np.fft.rfftn(r1)/ mesh.attrs['Nmesh'].prod()
		power = c1 * c1.conj() * BoxSize**3.
		power.flat[0] = 0
		
		W = numpy.empty(power.shape, dtype='f4') 
		W[...] = 2.0 
		W[..., 0] = 1.0 
		W[..., -1] = 1.0       
		
		power = power.flatten()
		W = W.flatten()
		
		kmax_for_3d_grid = 0.1
		cond = np.where((np.abs(kx) <= kmax_for_3d_grid) & (np.abs(ky) <= kmax_for_3d_grid)& (np.abs(kz) <= kmax_for_3d_grid))
		
		power = np.real(power)
		kx_cond = kx[cond]
		ky_cond = ky[cond]
		kz_cond = kz[cond]
		power_cond = power[cond]
		
		kx_cond = kx_cond.astype('float32')
		ky_cond = ky_cond.astype('float32')
		kz_cond = kz_cond.astype('float32')
		power_cond = power_cond.astype('float32')
		
		arr = np.array([kx_cond, ky_cond, kz_cond, power_cond]).T
		#arr = np.array([kx, ky, kz, power]).T
		arr.tofile('3d_power/%i.bin' % int(rank))		'''

		dk_actual = 0.01
		r = FFTPower(mesh,mode='2d',dk=dk_actual,kmin=0.0,kmax=k_Nyquist,Nmu=120,los=[0,0,1],poles=[0,2,4])
		Pkmu = r.power
		poles = r.poles

		f = open(output_path + output_name_pk,'w')
		f.write('# Krolewski mock challenge 1, input from ' + input_path + input_file + '\n')
		f.write('# Estimated shot noise subtracted from power spectra\n')
		f.write('# Estimated shot noise: %.5f\n' % (Pkmu[:,0].attrs['shotnoise']))
		f.write('# Code to generate this measurement in ' + __file__ + '\n')
		f.write('# Boxsize = %.1f\n'  % BoxSize)
		f.write('# Nmesh =  %i\n' % Nmesh)
		f.write('# Binning = ' + binning + '\n')
		f.write('# k mu pk Nmodes\n')
		for i in range(Pkmu.shape[1]):
			Pk = Pkmu[:,i]
			for j in range(len(Pk['k'])):
				f.write('%20.8e %20.8e %20.8e %i\n' % (Pk['k'][j], Pkmu.coords['mu'][i], Pk['power'][j].real-Pk.attrs['shotnoise'], 
					Pkmu.data['modes'][:,i][j]))
			Nk_actual = len(Pk['k'])

		f.close()

		f = open(output_path + output_name_pk_l,'w')
		f.write('# Krolewski mock challenge 1, input from ' + input_path + input_file + '\n')
		f.write('# Estimated shot noise subtracted from power spectra\n')
		f.write('# Estimated shot noise: %.5f\n' % (Pkmu[:,0].attrs['shotnoise']))
		f.write('# Code to generate this measurement in ' + __file__ + '\n')
		f.write('# Boxsize = %.1f\n'  % BoxSize)
		f.write('# Nmesh =  %i\n' % Nmesh)
		f.write('# Binning = ' + binning + '\n')
		f.write('# k P0 P2 P4 Nmodes\n')
		P0 = poles['power_0']
		P2 = poles['power_2']
		P4 = poles['power_4']

		for i in range(len(P0)):
			f.write('%20.8e %20.8e %20.8e %20.8e %i\n' % (poles['k'][i], P0.real[i]-poles.attrs['shotnoise'],
				P2.real[i], P4.real[i], poles.data['modes'][i]))

		f.close()
		
		'''# Write out P(kpar,kperp)
		KPAR = kz
		KPERP = np.sqrt(kx**2 + ky**2)
		
		# Define bins
		kpar = np.linspace(0,Nk_actual * dk_actual + dk_actual/2,Nk_actual+1) 
		kperp = np.linspace(0,Nk_actual * dk_actual + dk_actual/2,Nk_actual+1)
		
		# Sort by kpar and bin
		KPAR_argsort = np.argsort(KPAR)
		KPAR = KPAR[KPAR_argsort]
		KPERP = KPERP[KPAR_argsort]
		power = power[KPAR_argsort]
		W = W[KPAR_argsort]
		kpar_inds = np.searchsorted(KPAR, kpar, side='left')
		print('kpar_inds',kpar_inds)
		
		power_out = np.zeros((len(kpar)-1,len(kperp)-1))
		Nmodes_out = np.zeros((len(kpar)-1,len(kperp)-1))

		
		for i in range(len(kpar)-1):
			min = kpar_inds[i]
			max = kpar_inds[i+1]
			KPERP_sel = KPERP[min:max]
			power_sel = power[min:max]
			W_sel = W[min:max]
			
			KPERP_argsort = np.argsort(KPERP_sel)
			KPERP_sel = KPERP_sel[KPERP_argsort]
			power_sel = power_sel[KPERP_argsort]
			W_sel = W_sel[KPERP_argsort]
			kperp_inds = np.searchsorted(KPERP_sel, kperp, side='left')
			for j in range(len(kperp)-1):
				power_out[i,j] = np.sum(power_sel[kperp_inds[j]:kperp_inds[j+1]] * W_sel[kperp_inds[j]:kperp_inds[j+1]])/np.sum(W_sel[kperp_inds[j]:kperp_inds[j+1]])
				Nmodes_out[i,j] = kperp_inds[j+1]-kperp_inds[j]
				
			
		kpar_mesh,kperp_mesh = np.meshgrid(0.5 * (kpar[1:]+kpar[:-1]),0.5 * (kperp[1:] + kperp[:-1])) 
		#np.savetxt('test_rotate.txt',np.array([kpar_mesh.flatten(),kperp_mesh.flatten(),power_out.flatten(),Nmodes_out.flatten()]).T)
		
		kpar_mesh = kpar_mesh.flatten()
		kperp_mesh = kperp_mesh.flatten()
		power_out = power_out.flatten()
		Nmodes_out = Nmodes_out.flatten()
		
		f = open(output_path + output_name_pk_kpar_kperp,'w')
		f.write('# Krolewski mock challenge 1, input from ' + input_path + input_file + '\n')
		f.write('# Estimated shot noise subtracted from power spectra\n')
		f.write('# Estimated shot noise: %.5f\n' % (Pkmu[:,0].attrs['shotnoise']))
		f.write('# Code to generate this measurement in ' + __file__ + '\n')
		f.write('# Boxsize = %.1f\n'  % BoxSize)
		f.write('# Nmesh =  %i\n' % Nmesh)
		f.write('# Binning = ' + binning + '\n')
		f.write('# kpar kperp P(k) Nmodes\n')

		for i in range(len(kpar_mesh)):
			f.write('%20.8e %20.8e %20.8e %i\n' % (kpar_mesh[i], kperp_mesh[i], power_out[i], Nmodes_out[i]))
		f.close()'''

		print('Total time (s): ', time.time()-t0)
		