import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import os, sys
from scipy.interpolate import interp1d

sys.path.append('/albedo/home/mlosch/python')
sys.path.append('/albedo/home/mlosch/MITgcm/MITgcm/utils/python/MITgcmutils')
from myutils import *
from MITgcmutils import llc
import cmocean.cm as cmo
import xarray as xr
from xmitgcm import open_mdsdataset
from xmitgcm.utils import get_grid_from_input, get_extra_metadata
import datetime
import glob

bdir = "/albedo/work/projects/p_idemix_tr181/llc90"
runs = ["run16","run16_stormtide"]

def make_masks(coords, withoutArctic=True):
    global_mask = coords.hFacC.isel(k=0)
    # global_mask[6,:,:]=0. # delete Arctic face
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
    if withoutArctic:
        global_mask[6,:,:]=0. # delete Arctic face
        atlantic_mask[6,:,:]=0. # delete Arctic face
        indopacific_mask[6,:,:]=0. # delete Arctic face
    return global_mask, atlantic_mask, indopacific_mask

def compute_moc(wflux):
    # zonal integral
    wflx = zonal_sum(wflux)
    # order of integration: from north to south because of Atlantic MOC, requires sign change
    mocstrf = -np.flip(np.flip(wflx,axis=-1).cumsum(axis=-1),axis=-1)
    mocstrf[wflx==0]=0.
    return mocstrf

def flat2d(x):
    if type(x) is np.ndarray:
        x0 = np.concatenate( [np.concatenate([x[:,0,:,:],x[:,1,:,:],x[:,2,:,:]], axis=-2),
                              np.concatenate([x[:,3,:,:],x[:,4,:,:],x[:,5,:,:]], axis=-2)], axis=-1 )
        y0 = np.concatenate( [np.concatenate([x[:,7,:,:],x[:,8,:,:],x[:,9,:,:]], axis=-1),
                              np.concatenate([x[:,10,:,:],x[:,11,:,:],x[:,12,:,:]], axis=-1)], axis=-2 )
    else:
        x0 = xr.concat( [xr.concat( [x.isel(face=0),x.isel(face=1),x.isel(face=2)], dim = 'j' ),
                         xr.concat( [x.isel(face=3),x.isel(face=4),x.isel(face=5)], dim = 'j' )], dim='i' )
        y0 = xr.concat( [xr.concat( [x.isel(face=7),x.isel(face=8),x.isel(face=9)], dim = 'i' ),
                         xr.concat( [x.isel(face=10),x.isel(face=11),x.isel(face=12)], dim = 'i' )], dim='j' )
    return np.concatenate((x0,np.rot90(y0,k=1,axes=(-2,-1))), axis=-1)

def zonal_sum(fld):
    # zonal integral of scalar field
    return flat2d(fld).sum(axis=-1)


rdir = os.path.join(bdir,runs[0])
rdir1 = os.path.join(bdir,runs[1])

allfiles = glob.glob(os.path.join(rdir,'diags3D.*.meta'))

myiters = []
for f in allfiles:
    myiters.append(int(f.split('.')[-2]))

myiters.sort()


deltat=3600.
ds = open_mdsdataset(rdir, prefix=['diags3D'], iters=myiters[0],
                     delta_t=deltat, geometry='llc')
coords = ds.coords.to_dataset()
global_mask, atlantic_mask, indopacific_mask = make_masks(coords)

gdir = os.path.join(bdir,'grid')
rac=coords.rA.compute() #rdmds(os.path.join(gdir,'RAC'))
# drf=rdmds(os.path.join(gdir,'DRF'))
# rf=rdmds(os.path.join(gdir,'RF'))
# yc=rdmds(os.path.join(gdir,'YC'))
# yg=rdmds(os.path.join(gdir,'YG'))

refdate = datetime.datetime(780,1,1)
amoc_timeseries_0 = []
amoc_timeseries_1 = []
bmoc_timeseries_0 = []
bmoc_timeseries_1 = []
pmoc_timeseries_0 = []
pmoc_timeseries_1 = []
xdays = []
jmoc=197 # 26N
kmoc=33  # Z~-1500m
jy0=128
jy1=165
for it in myiters:
    d0 = open_mdsdataset(rdir, prefix=['diags3D'], iters=it,
                         delta_t=deltat, geometry='llc')
    d1 = open_mdsdataset(rdir1, prefix=['diags3D'], iters=it,
                         delta_t=deltat, geometry='llc')
    amoc0 = compute_moc(d0.WVEL*rac*atlantic_mask)
    pmoc0 = compute_moc(d0.WVEL*rac*indopacific_mask)
    amoc1 = compute_moc(d1.WVEL*rac*atlantic_mask)
    pmoc1 = compute_moc(d1.WVEL*rac*indopacific_mask)
    # now decide what is "amoc strength"
    # amoc at 26N
    # max between 500 and 1500
    amoc_timeseries_0.append(max(np.squeeze(amoc0)[:kmoc,jmoc]))
    amoc_timeseries_1.append(max(np.squeeze(amoc1)[:kmoc,jmoc]))
    # bottom cell min between 1500 to 6000
    bmoc_timeseries_0.append(min(np.squeeze(amoc0)[kmoc:,jmoc]))
    bmoc_timeseries_1.append(min(np.squeeze(amoc1)[kmoc:,jmoc]))
    # pmoc at 26N
    # min between 1500 to 6000, 30S and Eq
    pmoc_timeseries_0.append(min(np.squeeze(pmoc0)[kmoc:,jy0:jy1].ravel()))
    pmoc_timeseries_1.append(min(np.squeeze(pmoc1)[kmoc:,jy0:jy1].ravel()))
    xdays.append(refdate+datetime.timedelta(int(it/24)))
    print(xdays[-1],amoc_timeseries_0[-1]*1e-6,amoc_timeseries_1[-1]*1e-6,
          bmoc_timeseries_0[-1]*1e-6,bmoc_timeseries_1[-1]*1e-6,
          pmoc_timeseries_0[-1]*1e-6,pmoc_timeseries_1[-1]*1e-6)

import csv

header = ['days','CTRL AMOC','STORMTIDE AMOC','CTRL BMOC','STORMTIDE BMOC','CTRL PMOC','STORMTIDE PMOC']
data = [xdays,
        amoc_timeseries_0,amoc_timeseries_1,
        bmoc_timeseries_0,bmoc_timeseries_1,
        pmoc_timeseries_0,pmoc_timeseries_1]
fname = os.path.join(bdir,'postprocessing','timeseries_ABPMOC.csv')
with open(fname,'w') as f:
    write = csv.writer(f)
    write.writerow(header)
    write.writerows(np.asarray(data).transpose())

fig=plt.figure()
plt.plot(xdays,np.asarray(amoc_timeseries_0)*1e-6,label='CTRL')
plt.plot(xdays,np.asarray(amoc_timeseries_1)*1e-6,label='STORMTIDE')
plt.grid()
plt.legend()
plt.title('max(AMOC) at 26N / Sv')
plt.xlabel('calendar years')
fig.savefig(os.path.join(bdir,'timeseries_AMOC'))

plt.show()
