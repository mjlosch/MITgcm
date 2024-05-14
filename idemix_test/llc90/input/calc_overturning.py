import numpy as np
import sys
sys.path.append('/home/ollie/mlosch/MITgcm/MITgcm/utils/python/MITgcmutils')
sys.path.append('/home/ollie/mlosch/python')
import matplotlib.pyplot as plt
from MITgcmutils import rdmds
from MITgcmutils import llc
from myutils import sq, pcol
import os

myrun = 'test'
#bdir = '/work/ollie/mlosch/egerwing/llc270'
bdir = '/work/ollie/mlosch/idemix_test/llc90'
gdir = '/work/ollie/mlosch/idemix_test/llc90/grid'
fdir = '/home/ollie/mlosch/MITgcm/MITgcm/idemix_test/llc90/grid'
xg,yg=rdmds(os.path.join(gdir,'XG')),rdmds(os.path.join(gdir,'YG'))
xc,yc=rdmds(os.path.join(gdir,'XC')),rdmds(os.path.join(gdir,'YC'))

dxg=rdmds(os.path.join(gdir,'DXG'))
dyg=rdmds(os.path.join(gdir,'DYG'))
drf=rdmds(os.path.join(gdir,'DRF'))
rac=rdmds(os.path.join(gdir,'RAC'))
drf=rdmds(os.path.join(gdir,'DRF'))

#rf=rdmds(os.path.join(gdir,'RF'))[:,0,0]
rf=rdmds(os.path.join(bdir,myrun,'RF'))[:,0,0]
if np.sum(rf-np.abs(rf)) < 0:
    usingPCoords = False
else:
    usingPCoords = True

myiter = np.Inf
myiter = 35040
d3z,myiter,meta = rdmds(os.path.join(bdir,myrun,'diags3D'),myiter,
                        returnmeta=True)
dl = rdmds(os.path.join(bdir,myrun,'diagsLAYERS'),myiter)
nv,nl,ny,nx = dl.shape
nz = d3z.shape[-3]

#dlt=rdmds('layers1TH')
dlr=rdmds(os.path.join(bdir,myrun,'layers1RHO'))

#uflx=d3z[2,...]*np.tile(dyg.reshape((1,ny,nx)),(nz,1,1))*np.tile(drf,(1,ny,nx))
#vflx=d3z[3,...]*np.tile(dxg.reshape((1,ny,nx)),(nz,1,1))*np.tile(drf,(1,ny,nx))

def overturningw(w,rac):

    nz,ny,nx=w.shape
    dxdy=np.tile(rac.reshape((1,ny,nx)),(nz,1,1))
    wfx=llc.faces(w*dxdy)
    # integration over longitude
    wf = ( np.sum(wfx[0]+wfx[1],axis=-1)
         + np.sum(wfx[3]+wfx[4],axis=-2)[:,::-1] )
    # integration over latitudes
    psiw = np.cumsum(wf,axis=-1)

    return psiw

def overturning(udr,vdr,zlu,zlv,dxg,dyg):

    nl,ny,nx=vdr.shape
    ufx=llc.faces(udr*np.tile(dyg.reshape((1,ny,nx)),(nl,1,1)))
    vfx=llc.faces(vdr*np.tile(dxg.reshape((1,ny,nx)),(nl,1,1)))
    zzu = llc.faces(zlu)
    zzv = llc.faces(zlv)
    # get rid of Mediterrean Sea (does not help)
    vfx[0][:,600:670,110:] = 0.
    vfx[0][:,618:627,98:110] = 0.
    zzv[0][:,600:670,110:] = 0.
    zzv[0][:,618:627,98:110] = 0.
    zzu[0][:,600:670,110:] = 0.
    zzu[0][:,618:627,98:110] = 0.
    # integration over longitude
    vf =   np.sum(vfx[0]+vfx[1],axis=-1)[:,1:] \
         - np.sum(ufx[3]+ufx[4],axis=-2)[:,:0:-1]
    # average layer thickness over lontitude
    zl = - ( 0.5*(zzu[0]+zzu[1]).mean(axis=-1) \
           + 0.5*(zzv[3]+zzv[4]).mean(axis=-2)[:,::-1] ).cumsum(axis=0)
    # add a column of zeros at the southern end
    vf = np.hstack((np.zeros((nl,1)),vf))
    # integration over depth
    # psi = -np.cumsum(vf[::-1,:],axis=0)[::-1,:]
    psi = np.cumsum(vf,axis=0)

    return psi, zl

# overturning in z from wvelmass
if usingPCoords:
    cfac = 1. / (1035*9.81)
    zdir = -1.
else:
    cfac = 1.
    zdir = 1.

pfac = 1e-6 * cfac

psiw = overturningw(d3z[4,...],rac)
#psit = overturning(dl[1,...],dyg)
psir,zl = overturning(dl[0,...],dl[1,...],dl[2,...],dl[3,...],dxg,dyg)
yl = np.tile(yg[:3*nx,0].reshape((1,3*nx)),(nl,1))
zlmsk = np.where(zl,1,0)

levs = np.linspace(-30,20,26)
fig, ax = plt.subplots(3,1,sharex=True,sharey=False,squeeze=True,figsize=(8,10))
csf=np.copy(ax)
csf[0]=ax[0].contourf(yg[:3*nx,0],rf[:-1]*cfac*zdir,sq(psiw)*pfac*zdir,
                      levels=levs,extend='both')
csf[1]=ax[1].contourf(yg[:3*nx,0],sq(dlr[:-1,0,0]),sq(psir)*pfac,
                      levels=levs,extend='both')
#csf[2]=ax[2].pcolormesh(yg[:,0],dlt,sq(psit)*pfac)
csf[2]=ax[2].contourf(yl[:,:3*nx],np.array(zl[:,:3*nx])*cfac,
                      sq(psir*zlmsk)*pfac,
                      levels=levs, extend='both')
ax[0].set_title('%s: MOC (Sv) at timestep %i'%(myrun,myiter[0]))
ax[-1].set_xlabel('latitude (degN)')
ax[1].set_ylim(dlr.max(), dlr.min())
ax[1].set_ylim(dlr.max(), dlr[1:].min())
ylab=('depth (m)','$\sigma_{2}$ in kg/m$^3$',
      '$\sigma_{2}$ remapped to depth (m)')
for k,b in enumerate(ax):
    b.set_ylabel(ylab[k])
#    csf[k].set_clim([-30,20])
    plt.colorbar(csf[k],ax=ax[k],orientation='vertical',extend='both')


fig.show()
