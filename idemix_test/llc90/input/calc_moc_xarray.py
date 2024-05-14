import sys
sys.path.append('/home/ollie/mlosch/MITgcm/MITgcm/utils/python/MITgcmutils')
sys.path.append('/home/ollie/mlosch/python')

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import os
from xmitgcm import open_mdsdataset
from xmitgcm.utils import get_grid_from_input, get_extra_metadata

from myutils import *
import xarray as xr
import cmocean.cm as cmo

def mosaic_llc(field):
    return np.vstack([np.hstack([np.vstack([np.rot90(field[i]) for i in [9,8,7]]),
                                 np.vstack([np.rot90(field[i]) for i in [12,11,10]]),
                                 np.vstack([field[i] for i in [0,1,2]]),
                                 np.vstack([field[i] for i in [3,4,5]])]),
                      np.hstack([np.rot90(field[6])*np.tri(90)[::-1,:],np.triu(np.rot90(field[6],k=2))*np.tri(90)[::-1,:],
                                 np.triu(np.rot90(field[6],k=-1)),np.zeros(field[6].shape)])])[30:315,:]

def llc13to5faces(field):
    """
    fld = llc13to5faces(field) returns a list of 5 faces constructed from
    the input of a 13-faces field as returned from
    xmitgcm.open_mdsdataset(...,geometry='llc')
    """
    return [np.vstack((field[0,...],field[1,...],field[2,...])),
            np.vstack((field[3,...],field[4,...],field[5,...])),
            field[6,...],
            np.hstack((field[7,...],field[8,...],field[9,...])),
            np.hstack((field[10,...],field[11,...],field[12,...]))]

deltat=1800.
ny, nx = 1170, 90
prefix=['diags2D','diags3D','diagsGGL90','diagsKrN2']
# prefix=['diagsKrN2']

prefixmonthly = ['diags2DMonthly']
myiters=[3243360]

myrun = 'run03'
#bdir = '/work/ollie/mlosch/egerwing/llc270'
bdir = '/work/ollie/mlosch/idemix_test/llc90'
gdir = '/work/ollie/mlosch/idemix_test/llc90/grid'
fdir = '/home/ollie/mlosch/MITgcm/MITgcm/idemix_test/llc90/grid'
rdir0 = os.path.join(bdir,myrun)

llc90_extra_metadata = get_extra_metadata(domain='llc', nx=90)
grid = get_grid_from_input(os.path.join('/home/ollie/mlosch/MITgcm/nils/llc90/input','tile<NFACET>.mitgrid'),
                           geometry='llc', extra_metadata=llc90_extra_metadata)

ds0 = open_mdsdataset(rdir0,prefix=prefix,iters=myiters,
                     delta_t=deltat,ref_date="1958-1-1 0:0:0",geometry='llc')
ds0.coords['XC'] = grid.XC
ds0.coords['YC'] = grid.YC
ds0.coords['XG'] = grid.XG
ds0.coords['YG'] = grid.YG
coords = ds0.coords.to_dataset().reset_coords()
#ds = ds0.reset_coords(drop=True)
# fix coordinates

dl0 = open_mdsdataset(rdir0,prefix='diagsLAYERS',
                     delta_t=deltat,ref_date="1958-1-1 0:0:0",geometry='llc')


def compute_moc(wflux):
    u0 = xr.concat([xr.concat([wflux.isel(face=0),
                               wflux.isel(face=1),
                               wflux.isel(face=2)], dim='j' ),
                    xr.concat([wflux.isel(face=3),
                               wflux.isel(face=4),
                               wflux.isel(face=5)], dim='j' )], dim='i')
    v0 = xr.concat([xr.concat([wflux.isel(face=7),
                               wflux.isel(face=8),
                               wflux.isel(face=9)], dim='i' ),
                    xr.concat([wflux.isel(face=10)
                               ,wflux.isel(face=11),
                               wflux.isel(face=12)],dim='i')],dim='j')
    wflx = (np.concatenate((u0,np.rot90(v0,k=1,axes=(-2,-1))),
                           axis=-1)).sum(axis=-1)
    # order of integration: from north to south because of Atlantic
    # MOC, requires sign change
    mocstrf = -np.flip(np.flip(wflx,axis=-1).cumsum(axis=-1),axis=-1)
    mocstrf[wflx==0]=0.
    return mocstrf


def compute_moc_layers(dl,msk):

    locmsk = msk.values
    ufx = dl.LaUH1RHO*dl.dyG*locmsk
    vfx = dl.LaVH1RHO*dl.dxG*locmsk
    zzu = dl.LaHw1RHO*locmsk
    zzv = dl.LaHs1RHO*locmsk
    zzu = zzu.where(zzu>0.)
    zzv = zzv.where(zzv>0.)
    # integration over longitude
    vf =  xr.concat([vfx.isel(face= 0),vfx.isel(face= 1),vfx.isel(face= 2)],
                    dim='j_g').sum(dim='i') \
        + xr.concat([vfx.isel(face= 3),vfx.isel(face= 4),vfx.isel(face= 5)],
                    dim='j_g').sum(dim='i')
    uf =  xr.concat([ufx.isel(face= 7),ufx.isel(face= 8),ufx.isel(face= 9)],
                    dim='i_g').sum(dim='j') \
        + xr.concat([ufx.isel(face=10),ufx.isel(face=11),ufx.isel(face=12)],
                    dim='i_g').sum(dim='j')
    uvf = vf.data[:,1:]-uf.data[:,:0:-1]
    # average layer thickness over lontitude
    zv = xr.concat([
        xr.concat([zzv.isel(face= 0),zzv.isel(face= 1),zzv.isel(face= 2)],
                  dim='j_g'),
        xr.concat([zzv.isel(face= 3),zzv.isel(face= 4),zzv.isel(face= 5)],
                  dim='j_g')],dim='i').mean(dim='i', skipna=True)
    zu = xr.concat([
        xr.concat([zzu.isel(face= 7),zzu.isel(face= 8),zzu.isel(face= 9)],
                  dim='i_g'),
        xr.concat([zzu.isel(face=10),zzu.isel(face=11),zzu.isel(face=12)],
                  dim='i_g')],dim='j').mean(dim='j', skipna=True)
    # at this point we have averaged over 2 times 2 faces, now we need
    # to average over these two averages after flipping the
    # coordinates and shifiting them.
    # complicated reordering because of xarray dataArrays
    zur = zu.rename({'i_g':'j_g'}).assign_coords(j_g=range(zu.shape[1]))
    zvr = zv.assign_coords(j_g=range(zv.shape[1]))
    zl = xr.concat([zvr.roll(j_g=1),
                    zur.reindex(j_g=zur.j_g[::-1])],
                   dim='ii').mean(dim='ii',skipna=True).fillna(0.)
    zl = - zl.cumsum(axis=0)
    # add a column of zeros at the southern end
    uvf = np.hstack((np.zeros((uvf.shape[0],1)),uvf))
    # zl  = np.hstack((np.zeros((zl.shape[0],1)),zl))
    # integration over depth
    # psi = -np.cumsum(vf[::-1,:],axis=0)[::-1,:]
    psi = np.cumsum(uvf,axis=0)

    return psi, zl

mytime = -1
pfac = 1e-6
global_mask = coords.hFacC.isel(k=0).where(coords.YC<70).fillna(0)
atlantic_mask = coords.hFacC.isel(k=0).where(
    np.logical_and(coords.YC>-35,coords.YC<70)).where( # Southern Ocean and Arctic
    np.logical_and(coords.XC<20,coords.XC>-98)).where( # most of the non-Atlantic Ocean
    np.logical_or(coords.XC<0,np.logical_or(coords.YC<30,coords.YC>47))).where(
    np.logical_or(coords.XC<-9,np.logical_or(coords.YC<34,coords.YC>38))).where( # Strait of Gibraltar
    np.logical_or(coords.XC>-70,coords.YC>9)).where( # East Pacific
    np.logical_or(coords.XC>-84,coords.YC>14)).where( # Isthmus of Panama etc.
    np.logical_or(coords.XC>-90,coords.YC>18)).where(
    np.logical_or(coords.XC>-70,coords.YC<50)).fillna(0)
indopacific_mask = (global_mask-atlantic_mask).where(
    np.logical_and(coords.YC>-35,coords.YC<70)).fillna(0)
# remove Hudson
indopacific_mask[10,10:,:39] = 0.
# remove Med and parts of Arctic
indopacific_mask[ 2,20:,29:84] = 0.
# remove Bering strait and Chukchy Sea
indopacific_mask[ 7,:,:14] = 0.

msk=atlantic_mask
msk=global_mask
# psiw = compute_moc(ds0.WVEL.mean(dim='time')*(ds0.rA*msk))
# psir, zl = compute_moc_layers(dl0.mean(dim='time'),msk)
psiw = compute_moc(ds0.WVEL.isel(time=-1)*(ds0.rA*msk))
psir, zl = compute_moc_layers(dl0.isel(time=-1),msk)

zlmsk = np.ones(zl.shape)
#zlmsk = np.where(zl,1,0)

yg = xr.concat([coords.YG.isel(face= 0),coords.YG.isel(face= 1),
                coords.YG.isel(face= 2)], dim='j_g').isel(i_g=0)

nl = dl0.layer_1RHO_center.shape[0]
yl = np.tile(yg.data.reshape((1,yg.shape[0])),(nl,1))

vm=30
levs = np.linspace(-vm,vm,61)

fig, ax = plt.subplots(3,1,sharex=True,sharey=False,squeeze=True,figsize=(8,10))
for b in ax: b.cla()
csf=np.copy(ax)
# csf[0]=ax[0].contourf(yg,dl0.Z,(psiw)*pfac,
#                       levels=levs,extend='both',cmap=cmo.balance)
# csf[1]=ax[1].contourf(yg,sq(dl0.layer_1RHO_interface),sq(psir[:-1,:])*pfac,
#                       levels=levs,extend='both',cmap=cmo.balance)
# csf[2]=ax[2].contourf(yl,zl, sq(psir*zlmsk)*pfac, levels=levs, extend='both',
#                       cmap=cmo.balance)
csf[0]=ax[0].pcolormesh(yg,dl0.Z,(psiw)*pfac,
                        vmin=-vm,vmax=vm,cmap=cmo.balance)
csf[1]=ax[1].pcolormesh(yg,sq(dl0.layer_1RHO_interface),sq(psir[:-1,:])*pfac,
                        vmin=-vm,vmax=vm,cmap=cmo.balance)
csf[2]=ax[2].pcolormesh(yl,zl, sq(psir*zlmsk)*pfac,
                        vmin=-vm,vmax=vm,cmap=cmo.balance)

ax[1].invert_yaxis()
ax[0].set_title('%s: MOC (Sv)'%(myrun))
ax[-1].set_xlabel('latitude (degN)')
#ax[1].set_ylim(40,30)
ylab=('depth (m)','$\sigma_{2}$ in kg/m$^3$',
      '$\sigma_{2}$ remapped to depth (m)')
for k,b in enumerate(ax):
    b.set_ylabel(ylab[k])
#    csf[k].set_clim([-30,20])
    plt.colorbar(csf[k],ax=ax[k],orientation='vertical',extend='both')


fig.show()
