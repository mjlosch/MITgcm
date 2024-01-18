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
dvol = (coords.rA*coords.hFacC*coords.drF).compute()
total_vol = dvol.sum()
Cp = 3994.
rho0=1035.

refdate = datetime.datetime(780,1,1)
thc_timeseries_0 = []
thc_timeseries_1 = []
xdays = []

for it in myiters:
    d0 = open_mdsdataset(rdir, prefix=['diags3D'], iters=it,
                         delta_t=deltat, geometry='llc')
    d1 = open_mdsdataset(rdir1, prefix=['diags3D'], iters=it,
                         delta_t=deltat, geometry='llc')
    thc_timeseries_0.append(float(rho0*Cp*(d0.THETA*dvol).sum().compute()))
    thc_timeseries_1.append(float(rho0*Cp*(d1.THETA*dvol).sum().compute()))
    xdays.append(refdate+datetime.timedelta(int(it/24)))
    print(xdays[-1],thc_timeseries_0[-1],thc_timeseries_1[-1])

import csv

header = ['days','CTRL THC','STORMTIDE THC']
data = [xdays,
        thc_timeseries_0,thc_timeseries_1]
fname = os.path.join(bdir,'postprocessing','timeseries_THC.csv')
with open(fname,'w') as f:
    write = csv.writer(f)
    write.writerow(header)
    write.writerows(np.asarray(data).transpose())

# fig=plt.figure()
# plt.plot(xdays,np.asarray(thc_timeseries_0),label='CTRL')
# plt.plot(xdays,np.asarray(thc_timeseries_1),label='STORMTIDE')
# plt.grid()
# plt.legend()
# plt.title('max(AMOC) at 26N / Sv')
# plt.xlabel('calendar years')
# fig.savefig(os.path.join(bdir,'timeseries_AMOC'))

# plt.show()
