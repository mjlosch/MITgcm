import sys
sys.path.append('/albedo/home/mlosch/MITgcm/MITgcm/utils/python/MITgcmutils')
sys.path.append('/albedo/home/mlosch/python')

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
prefix=['diags3D']
# prefix=['diagsKrN2']

prefixmonthly = ['diags2DMonthly']
myiters=[3243360]
myiters = [17520]
myiters = [17520,35040,52608,70128,87648]
myiters = [5014032,
           5031600,
           5049120,
           5066640,
           5084160,
           5101728,
           5119248,
           5136768,
           5154288,
           5171856,
           5189376,
           5206896,
           5224416,
           5241984,
           5259504,
           5277024,
           5294544,
           5312112,
           5329632,
           5347152,
           5364672,
           5382240,
           5399760,
           5417280,
           5434800,
           5434800]
myiters = 'all'

myrun = 'run_kpp/tmp'
myrun = 'run05'
#bdir = '/albedo/work/users/mlosch/egerwing/llc270'
bdir = '/albedo/work/projects/p_idemix_test/idemix_test/llc90'
gdir = '/albedo/work/projects/p_idemix_test/idemix_test/llc90/grid'
fdir = '/albedo/home/mlosch/MITgcm/MITgcm/idemix_test/llc90/grid'
rdir0 = os.path.join(bdir,myrun)

llc90_extra_metadata = get_extra_metadata(domain='llc', nx=90)
grid = get_grid_from_input(os.path.join('/albedo/home/mlosch/MITgcm/nils/llc90/input','tile<NFACET>.mitgrid'),
                           geometry='llc', extra_metadata=llc90_extra_metadata)

ds = open_mdsdataset(rdir0,prefix=prefix,iters=myiters,
                     delta_t=deltat,ref_date="1710-1-1 0:0:0",geometry='llc')
ds.coords['XC'] = grid.XC
ds.coords['YC'] = grid.YC
ds.coords['XG'] = grid.XG
ds.coords['YG'] = grid.YG
coords = ds.coords.to_dataset().reset_coords()
#ds = ds.reset_coords(drop=True)
# fix coordinates

def make_masks(coords):
    latmax = 70
    global_mask = coords.hFacC.isel(k=0)
    global_mask[6,:,:]=0. # delete Arctic face
    global_mask[2,80:,60:]=0.
    global_mask[7,:,:13]=0.
    global_mask[10,:43,:11]=0.
    # remove Hudson
    global_mask[10,30:54,5:39] = 0.
    global_mask[10,30:62,10:39] = 0.
    #
    atlantic_mask = global_mask.where(coords.YC>-35).where( # Southern Ocean
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
    return global_mask, atlantic_mask, indopacific_mask


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


mytime = -1
ds0 = ds.isel(time=-1)
ds0 = ds.mean(dim='time')
pfac = 1e-6
global_mask, atlantic_mask, indopacific_mask = make_masks(coords)

psi_g = compute_moc(ds0.WVEL*(ds0.rA*global_mask))
psi_a = compute_moc(ds0.WVEL*(ds0.rA*atlantic_mask))
psi_p = compute_moc(ds0.WVEL*(ds0.rA*indopacific_mask))

yg = xr.concat([coords.YG.isel(face= 0),coords.YG.isel(face= 1),
                coords.YG.isel(face= 2)], dim='j_g').isel(i_g=0)

vm=15
levs = np.linspace(-vm,vm,61)

fig, ax = plt.subplots(3,1,sharex=True,sharey=True,squeeze=True,figsize=(8,10))
for b in ax: b.cla()
csf=np.copy(ax)
csf[0]=ax[0].pcolormesh(yg,ds.Z,(psi_g)*pfac,
                        vmin=-vm,vmax=vm,cmap=cmo.balance)
csf[1]=ax[1].pcolormesh(yg,ds.Z,(psi_a)*pfac,
                        vmin=-vm,vmax=vm,cmap=cmo.balance)
csf[2]=ax[2].pcolormesh(yg,ds.Z,(psi_p)*pfac,
                        vmin=-vm,vmax=vm,cmap=cmo.balance)

ax[0].set_title('%s: MOC (Sv)'%(myrun))
ax[-1].set_xlabel('latitude (degN)')
#ax[1].set_ylim(40,30)
ylab=('depth (m)','depth (m)','depth (m)')
for k,b in enumerate(ax):
    b.set_ylabel(ylab[k])
#    csf[k].set_clim([-30,20])
    plt.colorbar(csf[k],ax=ax[k],orientation='vertical',extend='both')


fig.show()
