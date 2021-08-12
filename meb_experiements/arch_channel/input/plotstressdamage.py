# convert to animation:
# ffmpeg -framerate 10 -pattern_type glob -i '*.png' -c:v libx264 -pix_fmt yuv420p stressdamage.mp4
# need this for batch script, but does not work
#import matplotlib as mpl
#mpl.use('Agg') 
#
from MITgcmutils import rdmds, wrmds
from myutils import *
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import os, sys
import datetime
from scipy.io.netcdf import netcdf_file
#import pdb
from cmocean import cm
import matplotlib.colors as colors
import matplotlib.dates as dates
import datetime

ncg= netcdf_file('grid.t001.nc','r')
xg = np.copy(ncg.variables['XG'][:])/1000.
yg = np.copy(ncg.variables['YG'][:])/1000.
ncg.close()
nc = netcdf_file('snapshot.0000000000.t001.nc','r')

time=np.copy(nc.variables['T'][:])
# s1=np.copy(nc.variables['SIsig1'][:])
# s2=np.copy(nc.variables['SIsig2'][:])
# d =np.copy(nc.variables['SIdamage'][:])
# hf=np.copy(nc.variables['SIheff'][:])

refdate=datetime.datetime(1,1,1,0,0)
fdir = "../figs/"+os.path.split(os.getcwd())[-1]

shareaxes=True
stressmax = 2e-4;
#stressmax = 3
horvert = 'vertical'
fig, ax =plt.subplots(2,2,sharex=shareaxes,sharey=shareaxes,squeeze=True)
bx = ax.flatten()
fig.set_size_inches( (9,4), forward=True )

first=True
for t, mytime in enumerate(time):
#    if t<0: continue
    mydate=dates.num2date(mytime/86400. + 1.)
    for k in range(len(bx)): bx[k].cla()
    s1=np.copy(nc.variables['SIsig1'][t,0,:,:])
    s2=np.copy(nc.variables['SIsig2'][t,0,:,:])
    d =np.copy(nc.variables['SIdamage'][t,0,:,:])
    hf=np.copy(nc.variables['SIheff'][t,0,:,:])
    csf1 = bx[0].pcolormesh(xg, yg, sq(s1+s2),
                            vmin=-stressmax,vmax=stressmax,cmap=cm.delta)
    csf2 = bx[1].pcolormesh(xg, yg, sq(s1-s2),
                            norm=colors.LogNorm(vmin=1e-6, vmax=stressmax))
    csf3 = bx[2].pcolormesh(xg, yg, sq(d), vmin=0., vmax=1.)
    csf4 = bx[3].pcolormesh(xg, yg, sq(hf),vmax=1.5,vmin=0.0)
    # csf1 = bx[0].pcolormesh(xg, yg, sq(s1+s2),[t,:,:],
    #                         vmin=-stressmax,vmax=stressmax,cmap=cm.delta)
    # csf2 = bx[1].pcolormesh(xg, yg, sq(s1-s2),[t,:,:],
    #                         norm=colors.LogNorm(vmin=1e-6, vmax=stressmax))
    # csf3 = bx[2].pcolormesh(xg, yg, sq(d)[t,:,:], vmin=0., vmax=1.)
    # csf4 = bx[3].pcolormesh(xg, yg, sq(hf)[t,:,:],vmax=1.5,vmin=.2)
    if first:
        cbh1 = plt.colorbar(csf1,ax=bx[0],orientation = horvert, extend='both')
        cbh2 = plt.colorbar(csf2,ax=bx[1],orientation = horvert, extend='max')
        cbh3 = plt.colorbar(csf3,ax=bx[2],orientation = horvert)
        cbh4 = plt.colorbar(csf4,ax=bx[3],orientation = horvert, extend='max')
        first=False

    #for k in range(len(ax)): bx[k].set_aspect('equal')
    bx[0].title.set_text('$\sigma_{I}$ (normalized)')
    bx[1].title.set_text('$\sigma_{II}$ (normalized)')
    bx[2].title.set_text('damage')
    bx[3].title.set_text('HEFF (m)')
    fig.suptitle(mydate.strftime('%Y/%m/%d %H:%M'))
    fname = 'stressdamage_' + mydate.strftime('%Y%m%d%H%M')
    print(fname)
    plt.draw()
    fig.savefig(os.path.join(fdir,fname))
#    fig.show(); plt.pause(.1)




