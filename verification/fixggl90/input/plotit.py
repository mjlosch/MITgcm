import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cmocean.cm as cmo

from xmitgcm import open_mdsdataset
import os

bdir="/Users/mlosch/MITgcm/MITgcm/verification/fixggl90"

runs = { 'zCoord': os.path.join(bdir,'z00'),
         'pCoord': os.path.join(bdir,'p00') }

def sq(a):
    import numpy as np
    a = np.squeeze(a)
    masked_array=np.ma.masked_where(a==0., a)
    return masked_array

dss = []
for name, path in runs.items():
    print(name, path)
    dss.append(open_mdsdataset(path,
                               prefix=['dynDiag','DiagMXL_3d','DiagMXL_2d'],
                               delta_t = 1200., ref_date="2001-01-01",
                               geometry='cartesian'))


nr = dss[0].Z.shape[0]

vname, vm = "GGL90ShP", 1e-9
pckwargs   = dict(norm = colors.LogNorm(vmin=1e-11,vmax=1e-4),
                  cmap=cmo.dense)

vname, vm = "GGL90TKE", 1e-7
pckwargs   = dict(norm = colors.LogNorm(vmin=1e-7,vmax=1e-1),
                  cmap=cmo.dense)

kt = -1

fig,ax = plt.subplots(1,2,sharex=True,sharey=True)

x = dss[0].XC.values/1000
z = dss[0].Z[1:]
#z = np.arange(nr)[1:]
cm=ax[0].pcolormesh(x,z,sq(dss[0][vname].isel(time=kt))[1:,:], **pckwargs)
cm=ax[1].pcolormesh(x,z,sq(dss[1][vname].isel(time=kt))[::-1,:][:-1,:],
                    **pckwargs)

fig.suptitle("%s at %s"%(vname,dss[0].time[kt].data.astype('datetime64[m]')))

# vm = 0.01
# pckwargs   = dict(norm = colors.Normalize(vmin=-vm,vmax=vm),
#                   cmap=cmo.balance)
# z = dss[0].Zp1[:-1]
# cm = ax[0].pcolormesh(x,z,sq(dss[0].UVEL.isel(time=kt)),**pckwargs)
# cm = ax[1].pcolormesh(x,z,sq(dss[1].UVEL.isel(time=kt))[::-1,:],**pckwargs)


for b in ax:
    b.grid()
    b.set_xlabel('X / km')
    if z.min()>=0: b.set_ylim([nr,0])

if z.min()>=0:
    ax[0].set_ylabel('vertical index')
else:
    ax[0].set_ylabel('Z / m')

ax[0].set_title('zCoord')
ax[1].set_title('pCoord')
cbarprops  = dict(orientation = 'horizontal')
fig.colorbar(cm,ax=ax,**cbarprops)

fname = "%s_masked"%vname
fig.savefig(fname)
