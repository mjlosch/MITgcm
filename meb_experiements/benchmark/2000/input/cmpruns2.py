# demonstrate fast generation of animation with FuncAnimation
#
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
from scipy.io.netcdf import netcdf_file
from cmocean import cm
import matplotlib.colors as colors
import matplotlib.dates as dates

# load grid data
ncg= netcdf_file('../run04/grid.t001.nc','r')
xg = np.copy(ncg.variables['XG'][:])/1000.
yg = np.copy(ncg.variables['YG'][:])/1000.
ncg.close()
# point to and open the netcdf files with the data
nc = netcdf_file('../run04/snapshot.0000000000.t001.nc','r')
nc2= netcdf_file('../run05/snapshot.tmp.cdf','r')

time=np.copy(nc.variables['T'][:])

# a small helper routine
def sq(a):
    a = np.squeeze(a)
    masked_array=np.ma.masked_where(a==0., a)
    return masked_array

shareaxes=True
stressmax = 3e-4;
horvert = 'horizontal'
fig, ax =plt.subplots(2,3,sharex=shareaxes,sharey=shareaxes,squeeze=True)
bx = ax.flatten()
fig.set_size_inches( (9,8), forward=True )
t = 0
# draw the plot first for t = 0
mydate=dates.num2date(time[t]/86400. + 1.)
sII=np.copy(nc.variables['SIsig1'][t,0,:,:]-
            nc.variables['SIsig2'][t,0,:,:])
d  =np.copy(nc.variables['SIdamage'][t,0,:,:])
hf =np.copy(nc.variables['SIheff'][t,0,:,:])
csf1 = bx[0].pcolormesh(xg, yg, sq(sII),
                        norm=colors.LogNorm(vmin=1e-6, vmax=stressmax))
csf2 = bx[1].pcolormesh(xg, yg, sq(hf),vmax=.35,vmin=0.0,cmap=cm.ice)
csf3 = bx[2].pcolormesh(xg, yg, sq(d), vmin=0., vmax=1.)
sII=np.copy(nc2.variables['SIsig1'][t,0,:,:]-
            nc2.variables['SIsig2'][t,0,:,:])
hf=np.copy(nc2.variables['SIheff'][t,0,:,:])
csf4 = bx[3].pcolormesh(xg, yg, sq(sII),
                        norm=colors.LogNorm(vmin=1e-3, vmax=1))
csf5 = bx[4].pcolormesh(xg, yg, sq(hf),vmax=.35,vmin=0.0,cmap=cm.ice)

# colorbars for each graph
cbh1 = plt.colorbar(csf1,ax=bx[0],orientation = horvert, extend='max')
cbh2 = plt.colorbar(csf2,ax=bx[1],orientation = horvert, extend='max')
cbh3 = plt.colorbar(csf3,ax=bx[2],orientation = horvert)
cbh4 = plt.colorbar(csf4,ax=bx[3],orientation = horvert, extend='max')
cbh5 = plt.colorbar(csf5,ax=bx[4],orientation = horvert, extend='max')

# some labels and titles
bx[0].set_ylabel('MEB')
bx[0].title.set_text('$\sigma_{II}$ (normalized)')
bx[1].title.set_text('HEFF (m)')
bx[2].title.set_text('damage')
bx[3].set_ylabel('VP')
bx[3].title.set_text('$\sigma_{II}$ (normalized)')
bx[4].title.set_text('HEFF (m)')
bx[5].set_axis_off()
th = fig.suptitle(mydate.strftime('%Y/%m/%d %H:%M'))

def animate(t):
    # In this routine you do whatever is necessary to update the data
    # in your plot. In my case it is loading new data and replacing it
    # in the pcolormesh plots with "set_array", and replacing the
    # string in the title by "set_text". Each plotted object has one of
    # methods to replace the data (e.g., "set_line" for plt.plot, etc.)
    sII=np.copy(nc.variables['SIsig1'][t,0,:,:]-
                nc.variables['SIsig2'][t,0,:,:])
    d  =np.copy(nc.variables['SIdamage'][t,0,:,:])
    hf =np.copy(nc.variables['SIheff'][t,0,:,:])
    csf1.set_array(sq(sII).ravel())
    csf2.set_array(sq(hf).ravel())
    csf3.set_array(sq(d).ravel())
    sII=np.copy(nc2.variables['SIsig1'][t,0,:,:]-
                nc2.variables['SIsig2'][t,0,:,:])
    hf =np.copy(nc2.variables['SIheff'][t,0,:,:])
    csf4.set_array(sq(sII).ravel())
    csf5.set_array(sq(hf).ravel())
    mydate=dates.num2date(time[t]/86400. + 1.)
    th.set_text(mydate.strftime('%Y/%m/%d %H:%M'))

# this creates the animation, which can then be save by anim.save('movie.mp4')
anim = FuncAnimation(
    fig, animate, interval=100, frames=len(time)-1)

plt.show()

# for t, mytime in enumerate(time):
#     animate(t)
#     plt.draw()
#     plt.show()

