#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-
######################## -*- coding: utf-8 -*-
"""Usage: plotres.py variable INPUTFILE(S)
"""

import sys, os
from getopt import gnu_getopt as getopt
import matplotlib.pyplot as plt
import numpy as np
import datetime
import glob

sys.path.append('/albedo/home/mlosch/MITgcm/MITgcm/utils/python/MITgcmutils')
sys.path.append('/albedo/home/mlosch/python')

from myutils import *
import cmocean.cm as cmo

bdir = '/albedo/work/projects/p_idemix_tr181/llc90'
files16  = glob.glob(os.path.join(bdir,'run16/stdout.*'))
files16s = glob.glob(os.path.join(bdir,'run16_stormtide/stdout.*'))

def get_output (fnames, mystring):
    """parse fname and get some numbers out"""
    timev = []
    myvar = []
    pp    = []
    for fname in fnames:
        try:
            f=open(fname)
        except:
            print(fname + " does not exist, continuing")
        else:
            for line in f:
                if "time_secondsf" in line:
                    ll = line.split()
                    timev.append(float(ll[-1].replace('D','e')))
                    myvar.append(np.NaN)

                if mystring in line:
                    ll = line.split()
                    myvar[-1] = float(ll[-1].replace('D','e'))

            f.close()

    timevs=np.asarray(timev)
    myvars=np.asarray(myvar)
    timevu, iunique = np.unique(timevs,return_index=True)
    isort = np.argsort(timevu)
    timevu=timevu[isort]
    myvaru=myvars[iunique][isort]

    return timevu, myvaru
# done

globArea = 3.579796597066180e+14
globVol  = 1.3349554967164e+18
Cp = 3994.
rho0=1035.

mystrs = ['dynstat_theta_mean','dynstat_salt_mean','seaice_heff_mean',
          'seaice_area_mean','dynstat_sst_mean','dynstat_sss_mean']
units = ['J','kg',r'km$^3$',r'km$^2$',r'$^\circ$C',' ']
varfac = [rho0*Cp*globVol,globVol,globArea*1e-9,globArea*1e-6,1,1]

refdate = datetime.datetime(780,1,1)
h16     = []
xday16  = []
h16s    = []
xday16s = []
for mystr in mystrs:
    print(mystr)
    timesec, h = get_output(files16, mystr)
    h16.append(h)
    timeday = np.asarray(timesec)/86400.
    xday16.append(
        np.array([refdate + datetime.timedelta(days=i) for i in timeday]))
    timesec, h = get_output(files16s, mystr)
    h16s.append(h)
    timeday = np.asarray(timesec)/86400.
    xday16s.append(
        np.array([refdate + datetime.timedelta(days=i) for i in timeday]))


# save extracted fields to ascii files
import csv

header = mystrs.copy()
header.insert(0,'days')
unitsh = units.copy()
unitsh.insert(0,' ')
varfch = varfac.copy()
varfch.insert(0,' ')

data = h16.copy()
data.insert(0,xday16[0])
fname = os.path.join(bdir,'timeseries_CTRL.csv')
with open(fname,'w') as f:
    write = csv.writer(f)
    write.writerow(header)
    write.writerow(unitsh)
    write.writerow(varfch)
    write.writerows(np.asarray(data).transpose())

data = h16s.copy()
data.insert(0,xday16s[0])
fname = os.path.join(bdir,'timeseries_STORMTIDE.csv')
with open(fname,'w') as f:
    write = csv.writer(f)
    write.writerow(header)
    write.writerow(unitsh)
    write.writerow(varfch)
    write.writerows(np.asarray(data).transpose())

fig,ax = plt.subplots(3,2,sharex=True,figsize=(12,9))

ax=ax.ravel()
for k, mystr in enumerate(mystrs):
    ax[k].plot(xday16[k],h16[k]*varfac[k],label='CTRL',linewidth=0.5)
    ax[k].plot(xday16s[k],h16s[k]*varfac[k],label='STORMTIDE',linewidth=0.5)
    ax[k].set_title("%s / %s"%(mystr,units[k]))

for b in ax[-2:]: b.set_xlabel('calendar year')

for b in ax: b.grid()

ax[-1].legend()

fig.savefig(os.path.join(bdir,'timeseries'))
plt.show()
