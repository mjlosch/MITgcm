#Analysis using XMITGCM

from xmitgcm import open_mdsdataset
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
data_dir = "../run_2"
ds = open_mdsdataset(data_dir)
print("Es wurden ", ds.time.size, "Zeitschritte gemacht von ", pd.to_timedelta(ds.time[0].values), "zu ", pd.to_timedelta(ds.time[-1].values))
index = int(input("Plots nach wie vielen Zeitschritten?"))
index = index -1


#divergence finite differenzen metode

#mit diff-Funktion
diff_u = ds.SIuice.diff("XG") # YC:65, XG:64
div_u = diff_u/ 8000
#div_1 = diff_2.dot(1/ds.dxF, dims=["YC"])

diff_v = ds.SIvice.diff("YG")
div_v = diff_v/ 8000

div_1 = np.resize(div_u[index].values, [65,65]) #resizing klappt an dieser Stelle nicht. Besser vorher fixen
div_2 = np.resize(div_v[index].values, [65,65])
div_mit_diff_ =  div_1 + div_2

#mit np.roll()
div_1 = (np.roll(ds.SIuice[index].values, -1,1) -  ds.SIuice[index].values)/ds.dxF.values + (np.roll(ds.SIvice[index].values, -1,0) - ds.SIvice[index].values)/ds.dyF.values


plt.figure()

plt.pcolormesh(div_1, cmap="coolwarm", vmin= -0.5e-5, vmax=0.5e-5)
plt.xlabel(ds.XC.long_name)
plt.ylabel(ds.YC.long_name)
plt.colorbar(label="Divergenz")
title = "iter =  "+ str(ds.iter[index].values) + ", time: " + str(pd.to_timedelta(ds.time[index].values))
plt.title(title)
plt.tight_layout()
plt.savefig("/home/csys/mbourget/Desktop/Plots/div_1.pdf")
plt.show()
plt.cla()


plt.pcolormesh(div_mit_diff_, cmap="coolwarm", vmin= -0.5e-5, vmax=0.5e-5)
plt.xlabel(ds.XC.long_name)
plt.ylabel(ds.YC.long_name)
#plt.colorbar(label="Divergenz")
title = "iter =  "+ str(ds.iter[index].values) + ", time: " + str(pd.to_timedelta(ds.time[index].values))
plt.title(title)
plt.tight_layout()
plt.savefig("/home/csys/mbourget/Desktop/Plots/div_diff.pdf")
plt.show()
plt.cla()

#divergenz finite volumen methode
div_2 = (np.roll(ds.SIuice[index].values, -1,1) * np.roll(ds.dyG.values, -1,0) - ds.SIuice[index].values * ds.dyG.values + np.roll(ds.SIvice[index].values, -1,0) * np.roll(ds.dxG.values, -1,1) - ds.SIvice[index].values * ds.dxG.values)/(ds.dxF.values * ds.dyF.values)

plt.pcolormesh(div_2, cmap="coolwarm", vmin= -0.5e-5, vmax=0.5e-5)
plt.xlabel(ds.XC.long_name)
plt.ylabel(ds.YC.long_name)
plt.colorbar(label="Divergenz")
title = "iter =  "+ str(ds.iter[index].values) + ", time: " + str(pd.to_timedelta(ds.time[index].values))
plt.title(title)
plt.tight_layout()
plt.savefig("/home/csys/mbourget/Desktop/Plots/div_2.pdf")
plt.show()
plt.cla()


#Shear
#ds.SIshear[index].plot()
#plt.tight_layout()
#plt.show()
#plt.cla()
#
##effective hight
#ds.SIheff[index].plot()
#plt.tight_layout()
#plt.show()
#plt.cla()
#
##fractional ice cover
#ds.SIarea[index].plot()
#plt.tight_layout()
#plt.show()
#plt.cla()



#plt.contourf(ds.XC/1000, ds.YC/1000, ds.Eta.isel(time=-1), np.linspace(-0.02, 0.05,8), cmap='hot_r')
#plt.colorbar()
#plt.show()

#Analysis using MITgcmutils
#from MITgcmutils import mds
#import numpy as np
#import matplotlib.pyplot as plt
#import os
#experiment = "/scratch/users/mbourget/MITgcm/verification/martin_8000/run_2"
#XC = mds.rdmds(os.path.join(experiment,'XC'))
#YC = mds.rdmds(os.path.join(experiment,'YC'))
#Eta = mds.rdmds('Eta', 7200)
#plt.contourf(XC/1000, YC/1000, Eta, np.linspace(-0.02, 0.05,8), cmap='hot_r')
#plt.colorbar()
#plt.show()

#Beide oben verwendeten Analysen laufen. Ich finde die obere mit der Verwendung von xarays angenehmer
#diagnostic files fehlen noch