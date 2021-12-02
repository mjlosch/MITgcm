import numpy as np

def strainrates( ui, vi, dxg, dyg ):
    dxv = 0.5*(dxg+np.roll(dxg,1,1))
    dyu = 0.5*(dyg+np.roll(dyg,1,0))
    dxf = 0.5*(dxg+np.roll(dxg,-1,0))
    dyf = 0.5*(dyg+np.roll(dyg,-1,1))
    e11=(np.roll(ui,-1,1)-ui)/dxf
    e22=(np.roll(vi,-1,0)-vi)/dyf
    e12sql=0.25*( (ui-np.roll(ui,1,0))/dyu + (vi-np.roll(vi,1,1))/dxv )**2
    e12sq = 0.25*(e12sql + np.roll(e12sql,-1,0) + np.roll(e12sql,-1,1) 
                  + np.roll(np.roll(e12sql,-1,0),-1,1))
    return e11, e22, e12sq

def deformation( e1, e2, e12sq ):
    
    divergnc = e1
    sheardef = np.sqrt(e2**2 + 4.*e12sq)

    return divergnc, sheardef

def sheardef( e2, e12sq): 
    return np.sqrt(e2**2+4.*e12sq)

def ice_strength(heff,conc,pstar,cstar):
    return pstar*heff*np.exp(-cstar*(1.-conc))

def Delta(e1,e2,e12sq,recip_e2):
    return np.sqrt(e1**2 + (e2**2 + 4.*e12sq)*recip_e2)

def DeltaReg(Delta,DeltaMin):
    return np.where(Delta<DeltaMin,DeltaMin,Delta)

def principal_stress(sigp,sigm,sig12sq):
    sig1 = .5*(sigp + np.sqrt( sigm**2 + 4.*sig12sq ))
    sig2 = .5*(sigp - np.sqrt( sigm**2 + 4.*sig12sq ))
    return sig1, sig2

def isfastice(ui,vi,ar):
    spcrit = 5.e-4 # m/s corresponds to displacement of 600m in 2 weeks 
                   # (Lemieux et al. 2015), used for 2 week averages
    spd = np.sqrt( (ui+np.roll(ui,-1,1))**2+(vi+np.roll(vi,-1,0))**2 )
    spd = np.ma.masked_array( spd, np.logical_or(ar <= 0.0, spd==0.) )
    return spd<=spcrit

