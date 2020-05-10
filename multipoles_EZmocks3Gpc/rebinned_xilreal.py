import os
import numpy as np
import pandas as pd


for i in range(0,10,1):
	df = pd.read_csv('multipoles/xilreal_00'+str(i)+'.txt',delim_whitespace=True, header=None)
	df.columns =['s','mono', 'quad', 'hexa']  
	bins = np.linspace(0., 200., 41)
	groups = df.groupby(np.digitize(df.s, bins))
	snew=groups.median().s
	mononew=groups.mean().mono
	quadnew=groups.mean().quad
	hexanew=groups.mean().hexa
	df_out = pd.DataFrame(snew)
	df_out['s'] = snew.values
	df_out['mono'] = mononew.values
	df_out['quad'] = quadnew.values
	df_out['hexa'] = hexanew.values

	np.savetxt('multipoles_5Mpc/xilreal_'+str(i)+'.txt', df_out.values,
               fmt='%.5f %.5f %.5f %.5f')


for i in range(10,100,1):
    df = pd.read_csv('multipoles/xilreal_0'+str(i)+'.txt',delim_whitespace=True, header=None)
    df.columns =['s','mono', 'quad', 'hexa']
    bins = np.linspace(0., 200., 41)
    groups = df.groupby(np.digitize(df.s, bins))
    snew=groups.median().s
    mononew=groups.mean().mono
    quadnew=groups.mean().quad
    hexanew=groups.mean().hexa
    df_out = pd.DataFrame(snew)
    df_out['s'] = snew.values
    df_out['mono'] = mononew.values
    df_out['quad'] = quadnew.values
    df_out['hexa'] = hexanew.values
    
    np.savetxt('multipoles_5Mpc/xilreal_'+str(i)+'.txt', df_out.values,
               fmt='%.5f %.5f %.5f %.5f')


for i in range(100,1000,1):
    df = pd.read_csv('multipoles/xilreal_'+str(i)+'.txt',delim_whitespace=True, header=None)
    df.columns =['s','mono', 'quad', 'hexa']
    bins = np.linspace(0., 200., 41)
    groups = df.groupby(np.digitize(df.s, bins))
    snew=groups.median().s
    mononew=groups.mean().mono
    quadnew=groups.mean().quad
    hexanew=groups.mean().hexa
    df_out = pd.DataFrame(snew)
    df_out['s'] = snew.values
    df_out['mono'] = mononew.values
    df_out['quad'] = quadnew.values
    df_out['hexa'] = hexanew.values
    
    np.savetxt('multipoles_5Mpc/xilreal_'+str(i)+'.txt', df_out.values,
               fmt='%.5f %.5f %.5f %.5f')


