import numpy as np
import sys
# Flux
Q=np.zeros((72,),dtype='float64')
Q[0::2] = 0.0003
Q[1::2] = 0.03

if sys.byteorder == 'little': Q.byteswap(True)

with open('salt.flux','wb') as fid: Q.tofile(fid)


