#Analysis using XMITGCM and XGCM

import xgcm
from xmitgcm import open_mdsdataset
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
data_dir = "../run_2"
ds = open_mdsdataset(data_dir)
print("Es wurden ", ds.time.size, "Zeitschritte gemacht von ", pd.to_timedelta(ds.time[0].values), "zu ", pd.to_timedelta(ds.time[-1].values))
index = int(input("Plots nach wie vielen Zeitschritten?"))
index = index -1

grid = xgcm.Grid(ds, periodic=["X","Y"], boundary="extrapolate")
metrics = {
    ('X',): ['dxC', 'dxG'], # X distances 
    ('Y',): ['dyC', 'dyG'], # Y distances 
    ('X', 'Y'): ['rA', 'rAz', 'rAs', 'rAw'] # Areas 
}
grid = xgcm.Grid(ds, metrics=metrics)

plt.figure()

#Divergenz finite differenzen methode
div_1 = grid.diff(ds.SIvice, 'Y') / ds.dyF + grid.diff(ds.SIuice, 'X') / ds.dxF

plt.pcolormesh(div_1[index], cmap="coolwarm", vmin= -0.5e-5, vmax=0.5e-5)
plt.xlabel(ds.XC.long_name)
plt.ylabel(ds.YC.long_name)
plt.colorbar(label="Divergenz")
title = "iter =  "+ str(ds.iter[index].values) + ", time: " + str(pd.to_timedelta(ds.time[index].values))
plt.title(title)
plt.tight_layout()
plt.savefig("/home/csys/mbourget/Desktop/Plots/div_xgcm_1.pdf")
plt.cla()
#
#v_ice = grid.interp(ds.SIvice, "Y")
#u_ice = grid.interp(ds.SIuice, "X")
#u_deriv = grid.derivative(u_ice, "X") 
#v_deriv = grid.derivative(v_ice, "Y")
#div_2 = u_deriv.values + v_deriv.values 
#print(div_2.ndim)
#plt.pcolormesh(div_2, cmap="warmcold")
#plt.xlabel(ds.XC.long_name)
#plt.ylabel(ds.YC.long_name)
#plt.colorbar(label="Divergenz")
#title = "iter =  "+ str(ds.iter[index].values) + ", time: " + str(pd.to_timedelta(ds.time[index].values))
#plt.title(title)
#plt.tight_layout()
#plt.show()
#plt.savefig("/home/csys/mbourget/Desktop/Plots/div_xgcm_2.pdf")
