import numpy as np
from myutils import *
from matplotlib.mlab import find
import matplotlib.pyplot as plt
ieee='b'
accuracy='float64'

f0=1.4e-4
gravity=9.81
Ho=1000
Lx=444e3
Ly=204e3
dx = 2e3;
dy = dx;
nx=int(Lx/dx)
ny=int(Ly/dy)
nz=3

x = (np.arange(nx,dtype = accuracy)+0.5)*dx;
y = (np.arange(ny,dtype = accuracy)+0.5)*dy;
xx,yy = np.meshgrid(x,y);

dxstr = "%i"%dx

# Flat bottom at z=-Ho
h=-Ho*np.ones((ny,nx),dtype = accuracy);
# channnel walls
h = np.where( ((yy<72e3) | (yy>Ly-72e3)),0,h)
writefield('bathy_channel_'+dxstr+'.bin',h)

h=-Ho*np.ones((ny,nx),dtype = accuracy);
# no channnel walls
# but islands
h = np.where((xx>117e3) & (xx<Lx-117e3) & ((yy<72e3) | (yy>Ly-72e3)),0,h)
writefield('bathy_3c_'+dxstr+'.bin',h)

h[0,:] = 0.
h[-1,:] =0.
writefield('bathyWithWall.bin',h)

h=-Ho*np.ones((ny,nx),dtype = accuracy);
h[:,:-2]=0.
writefield('bathyWall.bin',h)
# constant wind field
uwind = np.ones((xx.shape[0],xx.shape[1]),dtype=accuracy)*20.
writefield('Uwindfield_'+dxstr+'.bin',uwind)

# ocean is at rest

# initial thickness:
hice = 1.*np.ones((xx.shape[0],xx.shape[1]),dtype=accuracy)
writefield('thickness_'+dxstr+'.bin',hice)

# constant zeros
writefield('const_00_'+dxstr+'.bin',np.zeros(hice.shape))
