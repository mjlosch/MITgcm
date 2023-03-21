# run "python postproc.py" to process current working directory
import numpy as np
import os
from xmitgcm import open_mdsdataset
from xmitgcm.utils import get_grid_from_input, get_extra_metadata

import xarray as xr

cwd = os.getcwd()
thisrun = cwd.split('/')[-1]
prfx = thisrun.split('_')[0]

if "_" in thisrun:
    rname = 'dummy'
    if "jayne" in thisrun:       rname = 'JAYNE'
    elif "stormtide" in thisrun: rname = 'STORMTIDE'
    elif "nycander" in thisrun:  rname = 'NYCANDER'
    else: print("unknown extention in %s"%thisrun)
else:
    rname = 'CTRL'

rname = '%s_%s'%(prfx,rname)
print('processsing %s as %s'%(thisrun,rname))

# thisrun, rname = 'run17', 'CTRL'
bdir='/albedo/work/projects/p_idemix_tr181/llc90'
gdir='/albedo/work/projects/p_idemix_tr181/llc90/grid'
rdir=os.path.join(bdir,thisrun)
postprocdir = os.path.join(bdir,'postprocessing')

deltat=3600.
ny, nx = 1170, 90
year0,year1 = '1980','2020'
yearstr = year1

cycle = 5
startyear = 1958
refdate = "%i-1-1 0:0:0"%(startyear-(cycle-1)*62)
mycycle = "cycle%i"%(cycle)

print("preocessing %s (%s) for years %s to %s"%(thisrun,rname,year0,year1))

print("read diagsLAYERS: %s"%rdir)
dls = open_mdsdataset(rdir,prefix=['diagsLAYERS'],
                      delta_t=deltat,ref_date=refdate,geometry='llc')
myiters = dls.iter.values

dlm = dls.sel(time=slice(year0,year1)).mean(dim='time').compute()
dlm.to_netcdf(os.path.join(postprocdir,"%s_layers_%s.nc"%(rname,mycycle)))

print("read diags2D, 3D, GGL90, KrN2: %s"%rdir)
ds = open_mdsdataset(rdir,
                     prefix=['diags2D','diags3D','diagsGGL90','diagsKrN2'],
                     iters=myiters,
                     delta_t=deltat,ref_date=refdate,geometry='llc')

dsm = ds.sel(time=slice(year0,year1)).mean(dim='time').compute()
dsm.to_netcdf(os.path.join(postprocdir,"%s_%s.nc"%(rname,mycycle)))

print("read diags2DMonthly: %s"%rdir)
dm = open_mdsdataset(rdir,prefix=['diags2DMonthly'],
                     delta_t=deltat,ref_date=refdate,geometry='llc')

monthlyClim = dm.sel(time=slice(year0,year1)
                     ).groupby('time.month').mean(dim='time')
monthlyMaxMLD = dm.MXLDEPTH.sel(time=slice(year0,year1)
                                ).groupby('time.month').max(dim='time')
monthlyClim.to_netcdf(os.path.join(postprocdir,
                                   "%s_monthly_%s.nc"%(rname,mycycle)))
monthlyMaxMLD.to_netcdf(os.path.join(postprocdir,
                                     "%s_monthly_maxmld_%s.nc"%(rname,mycycle)))

print("read diagsTHflx: %s"%rdir)
ds = open_mdsdataset(rdir,
                     prefix=['diagsTHflx'],
#                     iters=myiters[1:],
                     delta_t=deltat,ref_date=refdate,geometry='llc')

dshf = ds.sel(time=slice(year0,year1)).mean(dim='time').compute()
dshf.to_netcdf(os.path.join(postprocdir,"%s_heatflux_%s.nc"%(rname,mycycle)))
