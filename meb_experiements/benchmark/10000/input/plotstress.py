#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from scipy.io.netcdf import netcdf_file

cohesion=25e3
intFrictCoeff    = 0.7
sqrtq = np.sqrt(intFrictCoeff**2 + 1) + intFrictCoeff
q = sqrtq**2

ElastMod = 650000.
sigc=2*cohesion/(np.sqrt(intFrictCoeff**2 + 1) - intFrictCoeff)/ElastMod
sigt=-sigc/q

nc=netcdf_file('snapshot.0000000000.t001.nc','r')
sig1 = nc.variables['SIsigI'][:].flatten()
sig2 = nc.variables['SIsigII'][:].flatten()
# sig1 = nc.variables['SIsigI'][-1,:,:,:].flatten()
# sig2 = nc.variables['SIsigII'][-1,:,:,:].flatten()
x = np.linspace(0,sig1.max(),2)
y = np.linspace(sig1.min(),0,2)
z = np.linspace(max(sig1.min(),sig2.min()),min(sig1.max(),sig2.max()),2)
# coordinate axes
plt.figure();
plt.clf();
plt.arrow(sig1.min(), 0., sig1.max()-sig1.min(), 0., color = 'k', width=1.e-5)
plt.arrow(0., min(sig2.min(),sigt), 0., sig2.max()-min(sig2.min(),sigt),
          color = 'k', width = 1.e-5)
plt.plot(sig1, sig2, '.')
plt.plot(x,(x-sigc)/q,'k-')
plt.plot(y,sigt*np.ones(y.shape),'k--')
plt.plot(z,z,'--',color='gray')

plt.grid()

nc.close()

#plt.savefig('stress')
plt.show()
