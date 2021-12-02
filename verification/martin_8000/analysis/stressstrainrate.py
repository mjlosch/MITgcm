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
#import pdb
from cmocean import cm
from stressroutines import *

kT = 0.0
DeltaMin = 1e-10
bdir = "/ace/user/mlosch/arctic/"
# exps = ["arctic_4km","arctic_9km","arctic_18km","arctic_36km"]
# rdir=["spinup03","run00/all","run00/all","run00/all"]
# exps = ["arctic_9km","arctic_18km","arctic_36km"]
# rdir=["run00/all","run00/all","run00/all"]
exps = ["arctic_9km"]
rdir=["run00/all"]

m = Basemap(projection='npstere',lon_0 = -45., boundinglat=68, resolution = 'i')

def myplotbase(m,ax,x,y,v,vmin,vmax,cmap):
    m.drawcoastlines(color = 'k',ax=ax)
    m.drawmapboundary(color='k',ax=ax)
    m.drawparallels(np.arange(50.,80,10.),ax=ax)
    m.drawmeridians(np.arange(0.,360.,30.),ax=ax)
    csf = ax.pcolormesh(x, y, sq(v), vmin=vmin,vmax=vmax,cmap=cmap)
    return csf

def plotstressstrain(myax,e1,e2,sig1,sig2,conc):
    #
    cond=conc>0.15
    R = sig1[cond]/sig2[cond]
    phi=np.arctan2((e1-e2)[cond],(e1+e2)[cond])
    # turn angles < 0 in angle+pi
#    phi[phi<0] += np.pi
#    phi = abs(phi)
    H,xedges,yedges=np.histogram2d(R.flatten(),phi.flatten()/np.pi,
                                   range=[[-1,1],[0,1]],bins=200)
    myax.pcolormesh(xedges,yedges,np.log10(sq(H)).T,cmap=cm.matter); 
    myax.plot([0.3,0.3],[0.,1.],'k--')
    myax.grid(True, which = 'both')
    myax.set_xlabel('$\sigma_1/\sigma_2$')
    myax.set_ylabel('tan$^{-1}(\epsilon_{II}/\epsilon_{I})/\pi$')
    return 

def plotstressstressII(myax,sig1,sig2,conc):
    #
    cond=np.logical_and(conc>0.15,sig1!=0)
    R = sig1[cond]/sig2[cond]
    sII = (sig1-sig2)[cond]*0.1
    H,xedges,yedges=np.histogram2d(R.flatten(),sII.flatten(),bins=200)
    myax.pcolormesh(xedges,yedges,np.log10(sq(H)).T,cmap=cm.matter); 
    myax.plot([0.3,0.3],[0.,20.],'k--')
    myax.grid(True, which = 'both')
    myax.set_xlabel('$\sigma_1/\sigma_2$')
    myax.set_ylabel('$\sigma_{II}$ X $10^{4}$ kg\,s$^{-2}$')
    return 

def plotstressstressI(myax,sig1,sig2,conc):
    #
    cond=np.logical_and(conc>0.15,sig1!=0)
    R = sig1[cond]/sig2[cond]
    sI = (sig1+sig2)[cond]*0.1
    H,xedges,yedges=np.histogram2d(R.flatten(),sI.flatten(),bins=200)
    myax.pcolormesh(xedges,yedges,np.log10(sq(H)).T,cmap=cm.matter); 
    myax.plot([0.3,0.3],[5.,-70.],'k--')
    myax.grid(True, which = 'both')
    myax.set_xlabel('$\sigma_1/\sigma_2$')
    myax.set_ylabel('$\sigma_{I}$ X $10^{4}$ kg\,s$^{-2}$')
    return 

def nsteps(mye):
    deltat = -1
    if mye.find('4km') >= 0: deltat = 240
    elif mye.find('9km') >= 0: deltat = 900
    elif mye.find('18km') >= 0: deltat = 1200
    elif mye.find('36km') >= 0: deltat = 1800
    else: print "invalid experiment :"+mye
    return 86400./deltat

refdate = datetime.datetime(1992, 1, 1, 0, 0)
def iter2date( myiter, mye ):
    mydate  = refdate + datetime.timedelta(myiter/nsteps(mye)-1)
    return mydate

def date2iter( mydate, mye ):
    myiter = ((mydate-refdate).days+1)*nsteps(mye)
    return myiter

def readgrid(gdir):
    xg   = rdmds(os.path.join(gdir,"XG"))
    yg   = rdmds(os.path.join(gdir,"YG"))
    dxg, dyg = rdmds(os.path.join(gdir,"DXG")), rdmds(os.path.join(gdir,"DYG"))
    return xg,yg,dxg,dyg

def readstate(rundir,it):
    conc = rdmds(os.path.join(rundir,'SIarea'),it)
    heff = rdmds(os.path.join(rundir,'SIheff'),it)
    ui   = rdmds(os.path.join(rundir,'SIuice'),it)
    vi   = rdmds(os.path.join(rundir,'SIvice'),it)
    return conc, heff, ui, vi

def stresstensor(heff,conc,e11,e22,e12sq,DeltaMin,kT):
    ep   = e11+e22
    em   = e11-e22
    delt = Delta(ep,em,e12sq)
    delr = DeltaReg(delt,DeltaMin)
    zeta = 0.5*ice_strength(heff,conc)*(1.+kT)/delr
    pr   = ice_strength(heff,conc)#*delt/delr
    sig1 = 2.*zeta * ep - pr*(1.-kT)
    sig2 = 2.*zeta*recip_e2 * em
    sig12sq = (2.*zeta*recip_e2)**2 * e12sq
    return sig1, sig2, sig12sq, pr

def strainrates2( ui, vi, dxc, dyc ):
    dxv = 0.5*(dxc+np.roll(dxc,1,0))
    dyu = 0.5*(dyc+np.roll(dyc,1,1))
    e11=(np.roll(ui,-1,1)-ui)/dxc
    e22=(np.roll(vi,-1,0)-vi)/dyc
    e12=0.25*( (ui-np.roll(ui,1,0))/dyu + (vi-np.roll(vi,1,1))/dxv )
    tmp = 0.25*(e12 + np.roll(e12,-1,0) + np.roll(e12,-1,1) 
                + np.roll(np.roll(e12,-1,0),-1,1))
    e12sq=tmp**2
    return e11, e22, e12sq

mydate = datetime.datetime(2001,12,16)
mydate = datetime.datetime(2002,01,20)
shareaxes = True
fig, ax =plt.subplots(2,2,sharex=True,sharey=False,squeeze=True)
fig.set_size_inches(14, 12, forward=True)
for ie in range(len(exps)):
    
    mye = exps[ie]
    myr = rdir[ie]
    
    it    = date2iter(mydate,mye)
    xg,yg,dxg,dyg = readgrid(os.path.join(bdir,mye,"grid"))
    conc,heff,ui,vi = readstate(os.path.join(bdir,mye,myr),it)

    p = ice_strength(heff,conc)
    e11,e22,e12sq = strainrates( ui, vi, dxg, dyg )
    sigp, sigm, sig12sq, pr = stresstensor(heff,conc,e11,e22,e12sq,DeltaMin,0.)
    sig1, sig2 = principal_stress(sigp,sigm,sig12sq)
    e1, e2 = principal_stress(e11+e22,e11-e22,e12sq)

    plotstressstrain(ax[0,ie],e1,e2,sig1,sig2,conc)
    plotstressstressII(ax[1,ie],sig1,sig2,conc)
    plotstressstressI(ax[1,1],sig1,sig2,conc)
    
#    x,y = m(xg,yg)
    csf=[]
    # csf.append(myplotbase(m,ax[ie,0],x,y,isf,0.,1.,cm.thermal))
    # ax[ie,0].title.set_text('fastice distribution')

# for k,mycsf in enumerate(csf):
#     pos1 = ax[ie,k].get_position()
#     cbax = fig.add_axes([pos1.x0, 0.05, pos1.x1-pos1.x0, 0.01]) 
#     cbh = plt.colorbar(mycsf, cax = cbax, orientation = 'horizontal',
#                        extend=myext[k])
#     cbh.ax.set_title(cblabs[k])

fig.suptitle(mydate.strftime('%Y/%m/%d'))
fig.show()
