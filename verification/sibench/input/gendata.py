import numpy as np
import random
import sys
sys.path.append('/home/csys/mlosch/python')
from myutils import *
#from matplotlib.mlab import find
import matplotlib.pyplot as plt
ieee='b'
accuracy='float64'

f0=1.4e-4
gravity=9.81
Ho=1000
c_min = 25e3
Lx=512e3
Ly=Lx
dx = 8e3;
dy = dx;
nx,ny,nz=int(Lx/dx)+1,int(Ly/dy)+1,3
x = (np.arange(nx,dtype = accuracy)+0.5)*dx;
y = (np.arange(ny,dtype = accuracy)+0.5)*dy;
xx,yy = np.meshgrid(x,y);
xu,yu = np.meshgrid(x-.5*dx,y);
xv,yv = np.meshgrid(x,y-.5*dy);

dxstr = "%i"%dx

# Flat bottom at z=-Ho
h=-Ho*np.ones((ny,nx),dtype = accuracy);
# channnel walls
h[:,-1]=0;
h[-1,:]=0;
writefield('bathy_3c_'+dxstr+'.bin',h)

#noisy cohesion
f = open("cohesion_2000.bin")
c = np.fromfile(f, ">f8")
c = np.reshape(c, (257,257))
c = c[::4,::4]
h= np.ones((ny,nx),dtype = accuracy)
i=0
while i < c.shape[0]:
    j = 0
    while j < c.shape[1]:
        h[i,j] = h[i,j] * c[i,j]
        j = j + 1
    i = i + 1

writefield("cohesion_random_"+dxstr+".bin",h)

variableWindField = True
shiftWindField = True
if variableWindField:
# variable wind field
    period = 16 # days
    if shiftWindField: period = 4 # days
    writeFreq = 3 # hours
    t = np.arange(0,period*24/writeFreq)/(24./writeFreq) # time in days
    vmax = 15.0; # maximale windgeschwindigkeit in m/s

    if shiftWindField:
        t = t + 4
        mx = Lx*.1 + Lx*.1*t
        my = Ly*.1 + Ly*.1*t
    else:
        tP=np.mod(t,period/2)
        tP[t>=period/2.]=period/2.-tP[t>=period/2.]
        tP = tP/(0.5*period)
        oLx=150.e3
        oLy=oLx
        mx = -oLx+(2*oLx+Lx)*tP
        my = -oLy+(2*oLy+Ly)*tP

    alpha0= np.pi/2. - np.pi/2./5. # 90 grad ist ohne Konvergenz oder Divergenz
    alpha = np.pi/2. - np.pi/2./5.*np.maximum(np.sign(np.roll(mx,-1)-mx),0.) \
            -np.pi/2./10.*np.maximum(np.sign(mx-np.roll(mx,-1)),0.)

    uwind = np.zeros((t.shape[0],xx.shape[0],xx.shape[1]))
    vwind = np.zeros((t.shape[0],yy.shape[0],yy.shape[1]))
    for k,myt in enumerate(t):
        wx =  np.cos(alpha[k])*(xx-mx[k]) + np.sin(alpha[k])*(yy-my[k])
        wy = -np.sin(alpha[k])*(xx-mx[k]) + np.cos(alpha[k])*(yy-my[k])
        r = np.sqrt((mx[k]-xx)*(mx[k]-xx)+(my[k]-yy)*(my[k]-yy))
        s = 1.0/50.e3*np.exp(-r/100.e3)
        if shiftWindField:
            w = np.tanh(myt*(8.0-myt)/2.)
        else:
            if myt<8:
                w = np.tanh(myt*(period/2.-myt)/2.)
            elif myt>=8 and myt<16:
                w = -np.tanh((myt-period/2.)*(period-myt)/2.)
        #    w = np.sin(2.0*np.pi*myt/period)

        # reset scaling factor w to one
        w = 1.
        uwind[k,:,:] = -wx*s*w*vmax;
        vwind[k,:,:] = -wy*s*w*vmax;

        spd=np.sqrt(uwind[k,:,:]**2+vwind[k,:,:]**2)
        div=uwind[k,1:-1,2:]-uwind[k,1:-1,:-2] \
             +vwind[k,2:,1:-1]-vwind[k,:-2,1:-1]
        # if spd.max() > 0:
        #     plt.clf();
        #     plt.subplot(211)
        #     pcol(xx/1e3,yy/1e3,sq(spd),vmax=vmax,vmin=0.)
        #     plt.axis('image')
        #     plt.colorbar();
        #     plt.quiver(xx/1e3,yy/1e3,uwind[k,:,:],vwind[k,:,:],pivot='middle')
        #     plt.title('time = '+str(myt))
        #     plt.subplot(212)
        #     pcol(xx[1:-1,1:-1]/1e3,yy[1:-1,1:-1]/1e3,sq(div),vmin=-1,vmax=1)
        #     plt.axis('image')
        #     plt.colorbar()
        #     plt.show();
        #     plt.pause(.01)

    if shiftWindField:
        writefield('Uwindfield_shifted_'+dxstr+'.bin',uwind)
        writefield('Vwindfield_shifted_'+dxstr+'.bin',vwind)
    else:
        writefield('Uwindfield_'+dxstr+'.bin',uwind)
        writefield('Vwindfield_'+dxstr+'.bin',vwind)

# ocean
uo = +0.01*(2*yy-Ly)/Ly
vo = -0.01*(2*xx-Lx)/Lx
writefield('uVel_'+dxstr+'.bin',uo)
writefield('vVel_'+dxstr+'.bin',vo)

# initial thickness:
hice = 0.3 + 0.005*np.sin(500*xx) + 0.005*np.sin(500*yy)
writefield('thickness_'+dxstr+'.bin',hice)
# symmetrize
hices = 0.5*(hice + hice.transpose())
writefield('thickness_sym_'+dxstr+'.bin',hices)
# initial thickness with random noise
hice = 0.3 + np.random.normal(scale=0.003,size=xx.shape)
writefield('noisy_thickness_'+dxstr+'.bin',hice)
# initial thickness for comparison with:
hice = 0.3 + 0.005*(np.sin(60./1000.e3*xx) + np.sin(30./1000.e3*yy))
writefield('thickness_aniso_'+dxstr+'.bin',hice)
hices = 0.5*(hice + hice.transpose())
writefield('thickness_aniso_sym_'+dxstr+'.bin',hices)

# constant
writefield('const_00_'+dxstr+'.bin',np.zeros(hice.shape))


# more initial conditions
r  = np.abs(0.2**2 - (xx/(2.*Lx)-0.25)**2 - (yy/(2.*Ly)-0.25)**2)
r1 = 2*(xx/(2.*Lx))**2 - 2* yy/(2.*Ly)
a0 = 1. - 0.5*np.exp(-800.*r) \
    - 0.4*np.exp(-90.*np.abs(r1+0.1)) \
    - 0.4*np.exp(-90*np.abs(r1+0.8))

a0[a0<0]=0.

hice = 0.3*a0
writefield('area_circle_'+dxstr+'.bin',a0)
writefield('thickness_circle_'+dxstr+'.bin',hice)

writefield('const_5.00_'+dxstr+'.bin',5.*np.ones(xx.shape))
writefield('const_0.8_'+dxstr+'.bin',.8*np.ones(xx.shape))
writefield('const_1.0_'+dxstr+'.bin',1.*np.ones(xx.shape))
