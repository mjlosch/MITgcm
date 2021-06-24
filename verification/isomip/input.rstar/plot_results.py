# simple python script to generate additional input data
import numpy as np
import matplotlib.pyplot as plt
from MITgcmutils import rdmds
from myutils import *
#
# fwflx=rdmds('SHICE_fwFlux',NaN)
# htflx=rdmds('SHICE_heatFlux',NaN)
# e=rdmds('Eta',NaN)
# t=rdmds('T',NaN)
# v=rdmds('V',NaN)
# u=rdmds('U',NaN)

myIter = np.Inf
myIter = 10
[d3,myIter,info] = rdmds('dynDiag',myIter,returnmeta=True)
s = d3[-2,:,:,:]
t = d3[-3,:,:,:]
u = d3[5,:,:,:]
v = d3[4,:,:,:]
w = d3[2,:,:,:]
rstar = d3[-1,:,:,:]

d2 = rdmds('surfDiag',myIter)
e = d2[0,:,:]
fwflx=d2[6,:,:]
htflx=d2[7,:,:]

yg=rdmds('YG')
drF=rdmds('drF')
rF=rdmds('rF')
# rsurf=rdmds('rSurfC')
# rlow=rdmds('rLowC')
# h0 = rsurf-rlow
# rfac = (e+h0)/np.where(h0==0,np.Inf,h0)
# drstar = rfac*drF
# rstar0 = cumsum(drstar[::-1,:,:],axis=0)[::-1,:,:] + rlow

ix = 20

fig, ax = plt.subplots(4,1,sharex=True,sharey=True)
csf = ax[0].pcolormesh(np.tile(yg[:,ix],(30,1)),rstar[:,:,ix],sq(t[:,:,ix]))
#                       vmin=-2.1)
plt.colorbar(csf, ax = ax[0], orientation='vertical')
csf = ax[1].pcolormesh(np.tile(yg[:,ix],(30,1)),rstar[:,:,ix],sq(s[:,:,ix]));
plt.colorbar(csf, ax = ax[1], orientation='vertical')
csf = ax[2].pcolormesh(np.tile(yg[:,ix],(30,1)),rstar[:,:,ix],sq(u[:,:,ix]));
plt.colorbar(csf, ax = ax[2],orientation='vertical')
csf = ax[3].pcolormesh(np.tile(yg[:,ix],(30,1)),rstar[:,:,ix],sq(v[:,:,ix]));
plt.colorbar(csf, ax = ax[3],orientation='vertical')

day = myIter[0]*1800/86400
ax[0].set_title('theta, day %i'%day)
ax[1].set_title('salinity')
ax[2].set_title('uvel')
ax[3].set_title('vvel')

ax[0].set_xlim([-80,-74])

plt.show()
#fig.savefig('isomip_rstar')
