# This is a matlab script that generates the input data
# variable x resolution

import numpy as np
import MITgcmutils as mit

prec='float64'
ieee='b'

rhoConst = 999.8
gravity = 9.81
grho = gravity*rhoConst

def writefield(fname,data):
    """Call signatures::

    writefield(filename, numpy.ndarray)

    Write unblocked binary data.
    """

    import sys

    if sys.byteorder == 'little': data.byteswap(True)

    fid = open(fname,"wb")
    data.tofile(fid)
    fid.close()

    # switch back to machine format
    if sys.byteorder == 'little': data.byteswap(True)

# Dimensions of grid
nx=3
ny=1
nz=26
nr = nz
# Vertical grid (meters)
dz = np.array([10., 10., 10., 10., 10., 11., 12., 14., 16., 18.,
     	       21., 24., 27., 31., 35.,
       	       40., 40., 40., 40., 40., 40., 40., 40., 40., 40., 40.]);

print('delR = ',nz,*dz, sep=", ")

print('delR = ',nz,*(dz[::-1]*grho), sep=", ")

zi = np.hstack([0,np.cumsum(dz)])
z = 0.5*(zi[:-1]+zi[1:])

# Initial temperature
gravity=9.81
talpha=2.0e-4
N2=2.e-5
Tz=N2/(gravity*talpha)

Sref=np.ones(nz,)*35
Tref=-Tz*(z - zi[-1])+2
print('Tref =', Tref)
Tref2d=np.tile(Tref,(nx,1)).transpose()
writefield('T_26.init',Tref2d)

writefield('T_26_p.init',Tref2d[::-1,...])

# Flux
Q=np.zeros((72,ny,nx));
Q[0:18,...] = 350;
Q[25:66,...] = -350*3/7;
writefield('Qnet_72.forcing',Q)

taux = 0.1*np.ones((72,ny,nx))
writefield('taux_72.forcing',taux)

# Uini:
uini = 0.5*np.ones((nz,ny,nx))
writefield('Uini.bin',uini)

# load anomaly
def sqinf(a):
    """ replace zeros by Inf
    """
    b = np.copy(np.squeeze(a))
    b[b==0] = np.inf
    return b

def calc_hydrostatic_pressure(s,t,p0,dz,gravity=9.81,rhoConst=1035.):
    from MITgcmutils import dens
    mskz = np.copy(t)
    mskz[mskz!=0]=1.
    dims=np.asarray(t.shape)
    dims[0]=dims[0]+1
    pf = np.zeros(dims)
    grho = gravity*rhoConst
    # initialize integration of non-linear hydrostatic equation
    dp = np.copy(dz)*grho
    p0 = np.copy(dz)*grho
    nr = Sref.shape[0]
    resid = 1
    while resid>1e-15:
        rhoInSitu = dens.mdjwf(s,t,p0/grho)*mskz
        # save old pressure
        dp = np.copy(p0)
        # compute new pressure
        pf[0,...] = 0.
        for k in range(nr):
            dpk =  dz[k,...]*gravity*rhoInSitu[k,...]
            p0[k,...]   = (pf[k,...] + 0.5*dpk)*mskz[k,...]
            pf[k+1,...] = (p0[k,...] + 0.5*dpk)*mskz[k,...]

        # check convergence
        dp = dp-p0
        resid = np.sqrt((dp**2).sum())
        print('hydrostatic pressure: resdiual = %e, '%np.sqrt((dp**2).sum()))

    print()
    return p0, pf, rhoInSitu

p0,pf,rhoInSitu = calc_hydrostatic_pressure(Sref,Tref,dz*grho,dz,
                                            gravity=gravity,
                                            rhoConst=rhoConst)

# layer thickness
print('delP =',nr,*np.diff(pf)[::-1], sep=", ")

dpc = np.zeros((nr,))
dpc[0]=p0[0] # assuming zero surface pressure
dpc[1:] = p0[1:]-p0[:-1]

# topography
dz3 = np.tile(dz.reshape((nz,1,1)),(1,ny,nx))
dp3 = np.tile(np.diff(pf).reshape((nz,1,1)),(1,ny,nx))

hfz = np.ones((nr,ny,nx))
# hfz[-1,0,1] = 0
hfz[-1,0,:] = 0
hfz[-2,0,1] = 0
#hfz[-3,0,1] = 0

H = -(hfz*dz3).sum(axis=0)

Hp = (hfz*dp3).sum(axis=0)

writefield('bathy.bin',H)
writefield('bathy_p.bin',Hp)

# load anomaly
hfp=hfz*np.tile(dpc.reshape((nz,1,1)),(1,ny,nx))
recip_rho = np.tile(1./sqinf(rhoInSitu).reshape((nz,1,1)),(1,ny,nx))
geopotanom = -((recip_rho - 1/rhoConst)*hfp).sum(axis=0)

writefield('geopotanom_mdjwf.bin',geopotanom)

# %%%%%
# tke = 1e-6*ones(1,nz);
# tke(1,1) = 1e-3;
# fid=fopen('TKE.init','w',ieee); fwrite(fid,tke,prec); fclose(fid);
