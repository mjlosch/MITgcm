#Analysis using XMITGCM

from xmitgcm import open_mdsdataset
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import numpy.ma as ma
import pandas as pd
data_dir_1 = "../plane_symmetric_vp"
data_dir_2 = "../plane_symmetric_meb_e8"
data_dir_3 = "../benchmark_meb_softer"
data_dir_4 = "../benchmark_meb_harder"
ds_1 = open_mdsdataset(data_dir_1)
ds_2 = open_mdsdataset(data_dir_2)
ds_3 = open_mdsdataset(data_dir_3)
ds_4 = open_mdsdataset(data_dir_4)
print("Es wurden ", ds_1.time.size, "Zeitschritte gemacht von ", pd.to_timedelta(ds_1.time[0].values), "zu ", pd.to_timedelta(ds_1.time[-1].values))
index = int(input("Plots nach wie vielen Zeitschritten?"))
index = index -1

#Mask
mask = np.logical_not(ds_1.hFacC)

fig = plt.figure()
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)


#divergenz finite volumen methode
div_vp = (np.roll(ds_1.SIuice[index].values, -1,1) * np.roll(ds_1.dyG.values, -1,0) - ds_1.SIuice[index].values * ds_1.dyG.values + np.roll(ds_1.SIvice[index].values, -1,0) * np.roll(ds_1.dxG.values, -1,1) - ds_1.SIvice[index].values * ds_1.dxG.values)/(ds_1.dxF.values * ds_1.dyF.values)
div_meb = (np.roll(ds_2.SIuice[index].values, -1,1) * np.roll(ds_2.dyG.values, -1,0) - ds_2.SIuice[index].values * ds_2.dyG.values + np.roll(ds_2.SIvice[index].values, -1,0) * np.roll(ds_2.dxG.values, -1,1) - ds_2.SIvice[index].values * ds_2.dxG.values)/(ds_2.dxF.values * ds_2.dyF.values)
#div_meb_softer = (np.roll(ds_3.SIuice[index].values, -1,1) * np.roll(ds_3.dyG.values, -1,0) - ds_3.SIuice[index].values * ds_3.dyG.values + np.roll(ds_3.SIvice[index].values, -1,0) * np.roll(ds_3.dxG.values, -1,1) - ds_3.SIvice[index].values * ds_3.dxG.values)/(ds_3.dxF.values * ds_3.dyF.values)
#div_meb_harder = (np.roll(ds_4.SIuice[index].values, -1,1) * np.roll(ds_4.dyG.values, -1,0) - ds_4.SIuice[index].values * ds_4.dyG.values + np.roll(ds_4.SIvice[index].values, -1,0) * np.roll(ds_4.dxG.values, -1,1) - ds_4.SIvice[index].values * ds_4.dxG.values)/(ds_4.dxF.values * ds_4.dyF.values)
div_vp = ma.masked_array(div_vp, mask)
div_meb = ma.masked_array(div_meb, mask)
#div_meb_softer = ma.masked_array(div_meb_softer, mask)
#div_meb_harder = ma.masked_array(div_meb_harder, mask)
#print(div_meb)
title = "iter =  "+ str(ds_1.iter[index].values) + ", time: " + str(pd.to_timedelta(ds_1.time[index].values))
fig.suptitle(title)

plot1 =ax1.pcolormesh(div_vp, cmap="coolwarm", vmin= -0.5e-5, vmax=0.5e-5)
#plot1 =ax1.pcolormesh(div_vp, cmap="coolwarm", vmin= np.max([div_vp.min(), div_meb_harder.min()]), vmax= np.min([div_vp.max(), div_meb_harder.max()]))
ax1.set_title("VP-Rheology")
ax1.set_xlabel(ds_1.XC.long_name)
ax1.set_ylabel(ds_1.YC.long_name)
ax1.axis("equal")
#print("Max und Min vp: ", div_vp.min(), " ", div_vp.max())
#print("Max und Min meb: ", div_meb.min(), " ", div_meb.max())
#print("Max und Min meb_softer: ", div_meb_softer.min(), " ", div_meb_softer.max())
#print("Max und Min meb_harder: ", div_meb_harder.min(), " ", div_meb_harder.max())
divider = make_axes_locatable(ax1)
cax = divider.append_axes('right', size='5%', pad=0.05)
cb = fig.colorbar(plot1, cax=cax, orientation="vertical", label="Divergenz")
#cb.remove()

plot2 =ax2.pcolormesh(div_meb, cmap="coolwarm",  vmin= -0.5e-5, vmax=0.5e-5)
ax2.set_title("MEB-Rheology (e8)")
ax2.set_xlabel(ds_2.XC.long_name)
ax2.set_ylabel(ds_2.YC.long_name)
ax2.axis("equal")

divider = make_axes_locatable(ax2)
cax = divider.append_axes('right', size='5%', pad=0.05)
cb = fig.colorbar(plot2,  cax=cax, orientation='vertical', label="Divergenz")

plt.tight_layout()
plt.savefig("/home/csys/mbourget/Desktop/Plots/simple_setup/div.pdf")
cb.remove()
plt.clf()

#effective hight
fig = plt.figure()
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
title = "iter =  "+ str(ds_1.iter[index].values) + ", time: " + str(pd.to_timedelta(ds_1.time[index].values))
fig.suptitle(title)

heff_vp = ma.masked_array(ds_1.SIheff[index], mask)
heff_meb = ma.masked_array(ds_2.SIheff[index], mask)

plot1 =ax1.pcolormesh(heff_vp, cmap="coolwarm", vmin = 0.2, vmax = 0.425)
ax1.set_title("VP-Rheology")
ax1.set_xlabel(ds_1.XC.long_name)
ax1.set_ylabel(ds_1.YC.long_name)
ax1.axis("equal")

divider = make_axes_locatable(ax1)
cax = divider.append_axes('right', size='5%', pad=0.05)
cb = fig.colorbar(plot1, cax=cax, orientation="vertical", label="Effective ice thickness")

plot2 =ax2.pcolormesh(heff_meb, cmap="coolwarm", vmin = 0.2, vmax = 0.425)
ax2.set_title("MEB-Rheology (e8)")
ax2.set_xlabel(ds_2.XC.long_name)
ax2.set_ylabel(ds_2.YC.long_name)
ax2.axis("equal")

divider = make_axes_locatable(ax2)
cax = divider.append_axes('right', size='5%', pad=0.05)
cb = fig.colorbar(plot2,  cax=cax, orientation='vertical', label="Effective ice thickness")

plt.tight_layout()
plt.savefig("/home/csys/mbourget/Desktop/Plots/simple_setup/heff.pdf")

cb.remove()
plt.clf()

#fractional ice cover
fig = plt.figure()
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
title = "iter =  "+ str(ds_1.iter[index].values) + ", time: " + str(pd.to_timedelta(ds_1.time[index].values))
fig.suptitle(title)

area_vp = ma.masked_array(ds_1.SIarea[index], mask)
area_meb = ma.masked_array(ds_2.SIarea[index], mask)

plot1 =ax1.pcolormesh(area_vp, cmap="coolwarm", vmin = 0.98)
ax1.set_title("VP-Rheology")
ax1.set_xlabel(ds_1.XC.long_name)
ax1.set_ylabel(ds_1.YC.long_name)
ax1.axis("equal")

divider = make_axes_locatable(ax1)
cax = divider.append_axes('right', size='5%', pad=0.05)
cb = fig.colorbar(plot1, cax=cax, orientation="vertical", label="Fractional ice-covered area")
#cb.remove()

plot2 =ax2.pcolormesh(area_meb, cmap="coolwarm", vmin = 0.98)
ax2.set_title("MEB-Rheology (e8)")
ax2.set_xlabel(ds_2.XC.long_name)
ax2.set_ylabel(ds_2.YC.long_name)
ax2.axis("equal")

divider = make_axes_locatable(ax2)
cax = divider.append_axes('right', size='5%', pad=0.05)
cb = fig.colorbar(plot2,  cax=cax, orientation='vertical', label="Fractional ice-covered area")

plt.tight_layout()
plt.savefig("/home/csys/mbourget/Desktop/Plots/simple_setup/area.pdf")

cb.remove()
plt.clf()
#
##Shear
#fig = plt.figure()
#ax1 = fig.add_subplot(121)
#ax2 = fig.add_subplot(122)
#title = "iter =  "+ str(ds_1.iter[index].values) + ", time: " + str(pd.to_timedelta(ds_1.time[index].values))
#fig.suptitle(title)
#
#shear_vp  = ma.masked_array(ds_1.SIshear[index], mask)
#shear_meb = ma.masked_array(ds_4.SIshear[index], mask)
#
#Z1 = ds_1.SIshear[index].values
#plot1 =ax1.pcolormesh(shear_vp, norm=colors.LogNorm(vmin=1e-13, vmax=Z1.max()), cmap="coolwarm")
#ax1.set_title("VP-Rheology")
#ax1.set_xlabel(ds_1.XC.long_name)
#ax1.set_ylabel(ds_1.YC.long_name)
#ax1.axis("equal")
#
#divider = make_axes_locatable(ax1)
#cax = divider.append_axes('right', size='5%', pad=0.05)
#cb = fig.colorbar(plot1, cax=cax, orientation="vertical", label="Shear deformation rate")
#
#Z2 = ds_4.SIshear[index].values
#plot2 =ax2.pcolormesh(shear_meb, norm=colors.LogNorm(vmin=1e-13, vmax=Z2.max()), cmap="coolwarm")
#ax2.set_title("MEB-Rheology (5e10)")
#ax2.set_xlabel(ds_2.XC.long_name)
#ax2.set_ylabel(ds_2.YC.long_name)
#ax2.axis("equal")
#
#divider = make_axes_locatable(ax2)
#cax = divider.append_axes('right', size='5%', pad=0.05)
#cb = fig.colorbar(plot2,  cax=cax, orientation='vertical', label="Shear deformation rate")
#plt.tight_layout()
#plt.savefig("/home/csys/mbourget/Desktop/Plots/constwind_2days/bm_shear_harder.pdf")
#
#cb.remove()
#plt.clf()
#
#Damage
plt.figure()
damage  = ma.masked_array(ds_2.SIdamage[index], mask)
plt.pcolormesh(damage, cmap="coolwarm",vmin = 0.4, vmax = 1.0)
plt.xlabel(ds_2.XC.long_name)
plt.ylabel(ds_2.YC.long_name)
cb = plt.colorbar(label="Damage")
title = "iter =  "+ str(ds_4.iter[index].values) + ", time: " + str(pd.to_timedelta(ds_4.time[index].values))
plt.title(title)
plt.tight_layout()
plt.savefig("/home/csys/mbourget/Desktop/Plots/simple_setup/meb_damage.pdf")

#Forcing

fig = plt.figure()
ax1 = fig.add_subplot(111)
title = "iter =  "+ str(ds_1.iter[index].values) + ", time: " + str(pd.to_timedelta(ds_1.time[index].values))
fig.suptitle(title)

forcing_x  = ma.masked_array(ds_1.SItaux[index], mask)
forcing_y  = ma.masked_array(ds_1.SItauy[index], mask)


plot1 = ax1.quiver(forcing_x, forcing_y, cmap="coolwarm")
ax1.set_title("zonal wind stress")
ax1.set_xlabel(ds_1.XC.long_name)
ax1.set_ylabel(ds_1.YC.long_name)
ax1.axis("equal")


plt.tight_layout()
plt.savefig("/home/csys/mbourget/Desktop/Plots/simple_setup/forcing.pdf")
plt.clf()


##Mask
##print(ds_1.hFacC.values)
#
#
#plt.pcolormesh(div_1, cmap="coolwarm", vmin= -0.5e-5, vmax=0.5e-5)
#plt.xlabel(ds.XC.long_name)
#plt.ylabel(ds.YC.long_name)
#cb = plt.colorbar(label="Divergenz")
#title = "iter =  "+ str(ds.iter[index].values) + ", time: " + str(pd.to_timedelta(ds.time[index].values))
#plt.title(title)
#plt.tight_layout()
#plt.savefig("/home/csys/mbourget/Desktop/Plots/div_1_softer.pdf")
##plt.show()
#plt.cla()
#cb.remove()
#
#plt.pcolormesh(div_mit_diff_, cmap="coolwarm", vmin= -0.5e-5, vmax=0.5e-5)
#plt.xlabel(ds.XC.long_name)
#plt.ylabel(ds.YC.long_name)
#cb = plt.colorbar(label="Divergenz")
#title = "iter =  "+ str(ds.iter[index].values) + ", time: " + str(pd.to_timedelta(ds.time[index].values))
#plt.title(title)
#plt.tight_layout()
#plt.savefig("/home/csys/mbourget/Desktop/Plots/div_diff_softer.pdf")
##plt.show()
#plt.cla()
#cb.remove()

#divergence finite differenzen metode

##mit diff-Funktion
#diff_u = ds.SIuice.diff("XG") # YC:65, XG:64
#div_u = diff_u/ 8000
##div_1 = diff_2.dot(1/ds.dxF, dims=["YC"])
#
#diff_v = ds.SIvice.diff("YG")
#div_v = diff_v/ 8000
#
#div_1 = np.resize(div_u[index].values, [65,65]) #resizing klappt an dieser Stelle nicht. Besser vorher fixen
#div_2 = np.resize(div_v[index].values, [65,65])
#div_mit_diff_ =  div_1 + div_2
#
##mit np.roll()
#div_1 = (np.roll(ds.SIuice[index].values, -1,1) -  ds.SIuice[index].values)/ds.dxF.values + (np.roll(ds.SIvice[index].values, -1,0) - ds.SIvice[index].values)/ds.dyF.values



#plt.contourf(ds.XC/1000, ds.YC/1000, ds.Eta.isel(time=-1), np.linspace(-0.02, 0.05,8), cmap='hot_r')
#plt.colorbar()
#plt.show()

#Analysis using MITgcmutil
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