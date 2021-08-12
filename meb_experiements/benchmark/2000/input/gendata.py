import numpy as np
from myutils import *
from matplotlib.mlab import find
import matplotlib.pyplot as plt
ieee='b'
accuracy='float64'

f0=1.4e-4
gravity=9.81
Ho=1000
Lx=500e3
Ly=Lx
dx = 2e3;
dy = dx;
nx,ny,nz=Lx/dx,Ly/dy,3
x = (np.arange(nx,dtype = accuracy)+0.5)*dx;
y = (np.arange(ny,dtype = accuracy)+0.5)*dy;
xx,yy = np.meshgrid(x,y);

dxstr = "%i"%dx

# Flat bottom at z=-Ho
h=-Ho*np.ones((ny,nx),dtype = accuracy);
# channnel walls
h[:,0]=0;
h[0,:]=0;
writefield('bathy_3c_'+dxstr+'.bin',h)

# variable wind field
period = 16 # days
writeFreq = 3 # hours
t = np.arange(0,period*24/writeFreq)/(24./writeFreq) # time in days
vmax = 15.0; # maximale windgeschwindigkeit in m/s

tP=np.mod(t,period/2)
tP[t>=period/2.]=period/2.-tP[t>=period/2.]
tP = tP/(0.5*period)
oLx=150.e3
oLy=oLx
mx = -oLx+(2*oLx+Lx)*tP
my = -oLy+(2*oLy+Ly)*tP
alpha = np.pi/2. - np.pi/2./5. # 90 grad ist ohne Konvergenz oder Divergenz
alpha = np.pi/2. - np.pi/2./5.*np.maximum(np.sign(np.roll(mx,-1)-mx),0.) \
        -np.pi/2./10.*np.maximum(np.sign(mx-np.roll(mx,-1)),0.)
uwind = np.zeros((t.shape[0],xx.shape[0],xx.shape[1]))
vwind = np.zeros((t.shape[0],yy.shape[0],yy.shape[1]))
for k,myt in enumerate(t):
    wx =  np.cos(alpha[k])*(xx-mx[k]) + np.sin(alpha[k])*(yy-my[k])
    wy = -np.sin(alpha[k])*(xx-mx[k]) + np.cos(alpha[k])*(yy-my[k])
    r = np.sqrt((mx[k]-xx)*(mx[k]-xx)+(my[k]-yy)*(my[k]-yy))
    s = 1.0/50.0*np.exp(-r/100.e3)
    if myt<8: w = np.tanh(myt*(period/2.-myt)/2.)
    elif myt>=8 and myt<16: w = -np.tanh((myt-period/2.)*(period-myt)/2.)
#    w = np.sin(2.0*np.pi*myt/period)

    uwind[k,:,:] = -wx*1.e-3*s*w*vmax;
    vwind[k,:,:] = -wy*1.e-3*s*w*vmax;

    spd=np.sqrt(uwind[k,:,:]**2+vwind[k,:,:]**2)
    div=uwind[k,1:-1,2:]-uwind[k,1:-1,:-2]+vwind[k,2:,1:-1]-vwind[k,:-2,1:-1]
    if spd.max() > 0:
        plt.clf(); 
        plt.subplot(211)
        pcol(xx/1e3,yy/1e3,sq(spd),vmax=vmax,vmin=0.)
        plt.axis('image')
        plt.colorbar(); 
        plt.quiver(xx/1e3,yy/1e3,uwind[k,:,:],vwind[k,:,:],pivot='middle')
        plt.title('time = '+str(myt))
        plt.subplot(212)
        pcol(xx[1:-1,1:-1]/1e3,yy[1:-1,1:-1]/1e3,sq(div))
        plt.colorbar()
        plt.show();
        plt.pause(.01)

#writefield('Uwindfield_'+dxstr+'.bin',uwind)
#writefield('Vwindfield_'+dxstr+'.bin',vwind)

# ocean 
uo = +0.01*(2*yy-Ly)/Ly
vo = -0.01*(2*xx-Lx)/Lx
writefield('uVel_'+dxstr+'.bin',uo)
writefield('vVel_'+dxstr+'.bin',vo)

# initial thickness:
hice = 0.3 + 0.005*np.sin(500*xx) + 0.005*np.sin(500*yy)
writefield('thickness_'+dxstr+'.bin',hice)
# initial thickness with random noise
hice = 0.3 + np.random.normal(scale=0.003,size=xx.shape)
writefield('noisy_thickness_'+dxstr+'.bin',hice)

# constant
writefield('const_00_'+dxstr+'.bin',np.zeros(hice.shape))
