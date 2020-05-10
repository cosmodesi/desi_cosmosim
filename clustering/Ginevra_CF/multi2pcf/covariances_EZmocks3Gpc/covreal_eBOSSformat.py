#Assembled (3x3) full covariance matrix and normalised one for the EZmocks
# output in eBOSS format

import os
import numpy as np
import pylab
from pylab import *
import csv
import pandas as pd
import pyfits




def split(array, nrows, ncols):
    """Split a matrix into sub-matrices."""
    r, h = array.shape
    return (array.reshape(h//nrows, nrows, -1, ncols)
                 .swapaxes(1, 2)
                 .reshape(-1, nrows, ncols))


#------------------------------------- read the covariances in (40 x 40) matrix format:


CMASS50JKx0 = np.loadtxt('covariances/covariance_xi0xi0_real.txt')
CMASS50JKx2 = np.loadtxt('covariances/covariance_xi2xi2_real.txt')
CMASS50JKx4 = np.loadtxt('covariances/covariance_xi4xi4_real.txt')
CMASS50JKx0x2 = np.loadtxt('covariances/covariance_xi0xi2_real.txt')
CMASS50JKx2x0=np.transpose(CMASS50JKx0x2)
CMASS50JKx0x4 = np.loadtxt('covariances/covariance_xi0xi4_real.txt')
CMASS50JKx4x0=np.transpose(CMASS50JKx0x4)
CMASS50JKx2x4 = np.loadtxt('covariances/covariance_xi2xi4_real.txt')
CMASS50JKx4x2=np.transpose(CMASS50JKx2x4)


CMASS50JK_full = np.squeeze(np.asarray(np.bmat([[CMASS50JKx0,CMASS50JKx0x2,CMASS50JKx0x4], [CMASS50JKx2x0,CMASS50JKx2,CMASS50JKx2x4], [CMASS50JKx4x0,CMASS50JKx4x2,CMASS50JKx4]])))

norm_CMASS50JK_full=np.zeros([len(CMASS50JK_full), len(CMASS50JK_full)])
column_CMASS50JK_full=np.zeros([len(CMASS50JK_full), len(CMASS50JK_full)])


outf=[]
outnormf=[]
for i in range(len(CMASS50JK_full)):
	for j in range(len(CMASS50JK_full)):
		norm_CMASS50JK_full[i][j]=CMASS50JK_full[i][j]/(np.sqrt(CMASS50JK_full[i][i]*CMASS50JK_full[j][j]))
		column_CMASS50JK_full[i][j]=CMASS50JK_full[i][j]
		outf.append([i,j,column_CMASS50JK_full[i][j]])
		outnormf.append([i,j,norm_CMASS50JK_full[i][j]])
        

df = pd.DataFrame(outf, columns=["row", "columns", "covariance"])      
dfnorm = pd.DataFrame(outnormf, columns=["row", "columns", "covariance"])      
np.savetxt('covariances/covariance_ful_real.txt', df.values,fmt='%ld %ld %.12f')
np.savetxt('covariances/normalised_covariance_full_real.txt', dfnorm.values,fmt='%ld %ld %.12f')
