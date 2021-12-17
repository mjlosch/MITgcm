from xmitgcm import open_mdsdataset
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.animation as animation
import numpy as np
import numpy.ma as ma
import pandas as pd
data_dir_1 = "../plane_symmetric_vp"
data_dir_2 = "../plane_symmetric_meb_e8"
data_dir_3 = "../benchmark_meb_softer"
data_dir_4 = "../benchmark_meb_harder"
data_dir_5 = "../test"

ds_1 = open_mdsdataset(data_dir_1)
ds_2 = open_mdsdataset(data_dir_2)
ds_3 = open_mdsdataset(data_dir_3)
ds_4 = open_mdsdataset(data_dir_4)
ds_5 = open_mdsdataset(data_dir_5)

#Mask
mask = np.logical_not(ds_1.hFacC)
mycmap = plt.cm.coolwarm

##Shear
#fig = plt.figure()
#ax1 = fig.add_subplot(121)
#ax2 = fig.add_subplot(122)
#
#k= 0
#Z_1 = ds_1.SIshear[0].values
#shear_vp = ds_1.SIshear.where(ds_1.SIshear!=0)
#Z_2 = ds_2.SIshear[0].values
#shear_meb = ds_2.SIshear.where(ds_2.SIshear!=0)
#
#pc_1 = ax1.pcolormesh(shear_vp[k],norm=colors.LogNorm(vmin=1e-13, vmax=np.max([Z_2.max(), Z_1.max()])), cmap="coolwarm")
#ax1.set_title("VP-Rheology")
#ax1.set_xlabel(ds_1.XC.long_name)
#ax1.set_ylabel(ds_1.YC.long_name)
#ax1.axis("equal")
#
#divider = make_axes_locatable(ax1)
#cax = divider.append_axes('right', size='5%', pad=0.05)
#cb = fig.colorbar(pc_1, cax=cax, orientation="vertical", label="Shear deformation rate")
#
#pc_2 = ax2.pcolormesh(shear_meb[k], norm=colors.LogNorm(vmin=1e-13, vmax=np.max([Z_2.max(), Z_1.max()])), cmap="coolwarm")
#ax2.set_title("MEB-Rheology")
#ax2.set_xlabel(ds_2.XC.long_name)
#ax2.set_ylabel(ds_2.YC.long_name)
#ax2.axis("equal")
#
#divider = make_axes_locatable(ax2)
#cax = divider.append_axes('right', size='5%', pad=0.05)
#cb = fig.colorbar(pc_2,  cax=cax, orientation='vertical', label="Shear deformation rate")
#plt.tight_layout()
#
#def animate(k):
#    pc_1.set_array(shear_vp[k].data.ravel())
#    pc_2.set_array(shear_meb[k].data.ravel())
#    title = "iter =  "+ str(ds_2.iter[k].values) + ", time: " + str(pd.to_timedelta(ds_2.time[k].values))
#    fig.suptitle(title)
#
#
#anim = animation.FuncAnimation(
#    fig, animate, interval=100, frames=ds_1.iter.values.max()-1)
# 
#anim.save('/home/csys/mbourget/Desktop/Plots/simple_setup/shear_compared.gif')
#plt.clf()

#one shear
#fig, ax = plt.subplots()
#Z = ds_5.SIshear[0].values
#shear_meb  = ma.masked_array(ds_5.SIshear[0], mask)
#pc = ax.pcolormesh(shear_meb,norm=colors.LogNorm(vmin=1e-13, vmax=Z.max()), cmap="coolwarm")
#fig.colorbar(pc, label="Shear deformation rate")
#fig.suptitle("MEB-Rheology")
#
#def animate(i):
#    Z = ds_5.SIshear[i].values
#    shear_meb = ma.masked_array(ds_5.SIshear[i], mask)
#    pc = ax.pcolormesh(shear_meb,norm=colors.LogNorm(vmin=1e-13, vmax=Z.max()), cmap="coolwarm")
#    return pc,
#
#im_ani = animation.FuncAnimation(fig, animate, frames=ds_5.iter.values.max()-1, interval=50, blit=True)
#im_ani.save('/home/csys/mbourget/Desktop/Plots/simple_setup/meb_shear_without_healing.gif')

#one damage
fig, ax = plt.subplots()
shear_meb  = ma.masked_array(ds_5.SIdamage[0], mask)
pc = ax.pcolormesh(shear_meb, cmap="coolwarm")
fig.colorbar(pc, label="Damage")
fig.suptitle("MEB-Rheology")

def animate(i):
    Z = ds_5.SIdamage[i].values
    shear_meb = ma.masked_array(ds_5.SIdamage[i], mask)
    pc = ax.pcolormesh(shear_meb, cmap="coolwarm")
    return pc,

im_ani = animation.FuncAnimation(fig, animate, frames=ds_5.iter.values.max()-1, interval=50, blit=True)
im_ani.save('/home/csys/mbourget/Desktop/Plots/simple_setup/meb_damage_without-healing.gif')

##Forcing
#fig, ax = plt.subplots()
#taux_vp = ma.masked_array(ds_1.SItaux[0], mask)
#pc = ax.pcolormesh(taux_vp, cmap="coolwarm")
#fig.colorbar(pc, label="Wind stress (u)")
#fig.suptitle("VP-Rheology")
#
#def animate(i):
#    taux_vp = ma.masked_array(ds_1.SItaux[i], mask)
#    pc = ax.pcolormesh(taux_vp, cmap="coolwarm")
#    return pc,
#
#im_ani = animation.FuncAnimation(fig, animate, frames=ds_1.iter.values.max()-1, interval=50, blit=True)
#im_ani.save('/home/csys/mbourget/Desktop/Plots/vp_taux_softer.gif')
#plt.clf()
#
##compare Ice strength parameter shear
#fig = plt.figure()
#ax0 = fig.add_subplot(232)
#ax1 = fig.add_subplot(234)
#ax2 = fig.add_subplot(235)
#ax3 = fig.add_subplot(236)
#
#k= 0
#Z_0 = ds_1.SIshear[0].values
#shear_0 = ds_1.SIshear.where(ds_1.SIshear!=0)
#Z_1 = ds_2.SIshear[0].values
#shear_1 = ds_2.SIshear.where(ds_2.SIshear!=0)
#Z_2 = ds_3.SIshear[0].values
#shear_2 = ds_3.SIshear.where(ds_3.SIshear!=0)
#Z_3 = ds_4.SIshear[0].values
#shear_3 = ds_4.SIshear.where(ds_4.SIshear!=0)
#
#pc_0 = ax0.pcolormesh(shear_0[k],norm=colors.LogNorm(vmin=1e-13, vmax=np.max([Z_0.max(), Z_2.max(), Z_1.max(), Z_3.max()])), cmap="coolwarm")
#ax0.set_title("VP")
#ax0.set_xlabel(ds_3.XC.long_name)
#ax0.set_ylabel(ds_3.YC.long_name)
#ax0.axis("equal")
#
#divider = make_axes_locatable(ax0)
#cax = divider.append_axes('right', size='5%', pad=0.05)
#cb = fig.colorbar(pc_0, cax=cax, orientation="vertical", label="Shear")
#
#pc_1 = ax1.pcolormesh(shear_2[k],norm=colors.LogNorm(vmin=1e-13, vmax=np.max([Z_0.max(), Z_2.max(), Z_1.max(), Z_3.max()])), cmap="coolwarm")
#ax1.set_title(r"$5e6$")
#ax1.set_xlabel(ds_3.XC.long_name)
#ax1.set_ylabel(ds_3.YC.long_name)
#ax1.axis("equal")
#
#divider = make_axes_locatable(ax1)
#cax = divider.append_axes('right', size='5%', pad=0.05)
#cb = fig.colorbar(pc_1, cax=cax, orientation="vertical", label="Shear")
#
#pc_2 = ax2.pcolormesh(shear_1[k],norm=colors.LogNorm(vmin=1e-13, vmax=np.max([Z_0.max(), Z_2.max(), Z_1.max(), Z_3.max()])), cmap="coolwarm")
#ax2.set_title(r"$5e8$")
#ax2.set_xlabel(ds_2.XC.long_name)
#ax2.set_ylabel(ds_2.YC.long_name)
#ax2.axis("equal")
#
#divider = make_axes_locatable(ax2)
#cax = divider.append_axes('right', size='5%', pad=0.05)
#cb = fig.colorbar(pc_2,  cax=cax, orientation='vertical', label="Shear")
#
#pc_3 = ax3.pcolormesh(shear_3[k],norm=colors.LogNorm(vmin=1e-13, vmax=np.max([Z_0.max(), Z_2.max(), Z_1.max(), Z_3.max()])), cmap="coolwarm")
#ax3.set_title(r"$5e10$")
#ax3.set_xlabel(ds_4.XC.long_name)
#ax3.set_ylabel(ds_4.YC.long_name)
#ax3.axis("equal")
#
#divider = make_axes_locatable(ax3)
#cax = divider.append_axes('right', size='5%', pad=0.05)
#cb = fig.colorbar(pc_3,  cax=cax, orientation='vertical', label="Shear")
#title = "iter =  "+ str(ds_2.iter[k].values) + ", time: " + str(pd.to_timedelta(ds_2.time[k].values))
#fig.suptitle(title)
#
#plt.tight_layout()
#
#def animate(k):
#    pc_0.set_array(shear_0[k].data.ravel())
#    pc_1.set_array(shear_2[k].data.ravel())
#    pc_2.set_array(shear_1[k].data.ravel())
#    pc_3.set_array(shear_3[k].data.ravel())
#    title = "iter =  "+ str(ds_2.iter[k].values) + ", time: " + str(pd.to_timedelta(ds_2.time[k].values))
#    fig.suptitle(title)
#
#anim = animation.FuncAnimation(
#    fig, animate, interval=100, frames=ds_2.iter.values.max()-1)
# 
#anim.save('/home/csys/mbourget/Desktop/Plots/constwind_2days/shear_compared.gif')

##compare Ice strength parameter damage
#fig = plt.figure()
#ax1 = fig.add_subplot(131)
#ax2 = fig.add_subplot(132)
#ax3 = fig.add_subplot(133)
#
#k= 0
#
#damage_1 = ds_2.SIdamage.where(ds_2.SIdamage!=0)
#damage_2 = ds_3.SIdamage.where(ds_3.SIdamage!=0)
#damage_3 = ds_4.SIdamage.where(ds_4.SIdamage!=0)
#
#pc_1 = ax1.pcolormesh(damage_2[k],vmin = 0.0, vmax = 1.0, cmap="coolwarm")
#ax1.set_title(r"$5e6$")
#ax1.set_xlabel(ds_3.XC.long_name)
#ax1.set_ylabel(ds_3.YC.long_name)
#ax1.axis("equal")
#
#divider = make_axes_locatable(ax1)
#cax = divider.append_axes('right', size='5%', pad=0.05)
#cb = fig.colorbar(pc_1, cax=cax, orientation="vertical", label="Damage")
#
#pc_2 = ax2.pcolormesh(damage_1[k],vmin = 0.0, vmax = 1.0, cmap="coolwarm")
#ax2.set_title(r"$5e8$")
#ax2.set_xlabel(ds_2.XC.long_name)
#ax2.set_ylabel(ds_2.YC.long_name)
#ax2.axis("equal")
#
#divider = make_axes_locatable(ax2)
#cax = divider.append_axes('right', size='5%', pad=0.05)
#cb = fig.colorbar(pc_2,  cax=cax, orientation='vertical', label="Damage")
#
#pc_3 = ax3.pcolormesh(damage_3[k],vmin = 0.0, vmax = 1.0, cmap="coolwarm")
#ax3.set_title(r"$5e10$")
#ax3.set_xlabel(ds_4.XC.long_name)
#ax3.set_ylabel(ds_4.YC.long_name)
#ax3.axis("equal")
#
#divider = make_axes_locatable(ax3)
#cax = divider.append_axes('right', size='5%', pad=0.05)
#cb = fig.colorbar(pc_3,  cax=cax, orientation='vertical', label="Damage")
#title = "iter =  "+ str(ds_2.iter[k].values) + ", time: " + str(pd.to_timedelta(ds_2.time[k].values))
#fig.suptitle(title)
#
#plt.tight_layout()
#
#def animate(k):
#    pc_1.set_array(damage_2[k].data.ravel())
#    pc_2.set_array(damage_1[k].data.ravel())
#    pc_3.set_array(damage_3[k].data.ravel())
#    title = "iter =  "+ str(ds_2.iter[k].values) + ", time: " + str(pd.to_timedelta(ds_2.time[k].values))
#    fig.suptitle(title)
#
#anim = animation.FuncAnimation(
#    fig, animate, interval=100, frames=ds_2.iter.values.max()-1)
# 
#anim.save('/home/csys/mbourget/Desktop/Plots/constwind_2days/damage_compared.gif')
#Damage
#fig, ax = plt.subplots()
#damage = ma.masked_array(ds_2.SIdamage[0], mask)
#pc = ax.pcolormesh(damage, cmap="coolwarm")
#fig.colorbar(pc, label="Damage")
#fig.suptitle("MEB-Rheology")
#
#def animate(i):
#    damage = ma.masked_array(ds_2.SIdamage[i], mask)
#    pc = ax.pcolormesh(damage, cmap="coolwarm")
#    return pc,
#
#im_ani = animation.FuncAnimation(fig, animate, frames=ds_1.iter.values.max()-1, interval=50, blit=True)
#im_ani.save('/home/csys/mbourget/Desktop/Plots/constwind_2days/meb_damage_softer.gif')
#plt.clf()

#compare Ice strength parameter
fig = plt.figure()
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

k= 0
Z_1 = ds_2.SIdamage[0].values
damage_1 = ds_2.SIdamage.where(ds_2.SIdamage!=0)
Z_2 = ds_3.SIdamage[0].values
damage_2 = ds_3.SIdamage.where(ds_3.SIdamage!=0)

pc_1 = ax1.pcolormesh(damage_1[k], cmap="coolwarm", vmin = 0.0, vmax = 1.0)
ax1.set_title(r"$5e8$")
ax1.set_xlabel(ds_1.XC.long_name)
ax1.set_ylabel(ds_1.YC.long_name)
ax1.axis("equal")

divider = make_axes_locatable(ax1)
cax = divider.append_axes('right', size='5%', pad=0.05)
cb = fig.colorbar(pc_1, cax=cax, orientation="vertical", label="Damage")

pc_2 = ax2.pcolormesh(damage_2[k], cmap="coolwarm",vmin = 0.0, vmax = 1.0)
ax2.set_title(r"$5e6$")
ax2.set_xlabel(ds_2.XC.long_name)
ax2.set_ylabel(ds_2.YC.long_name)
ax2.axis("equal")

divider = make_axes_locatable(ax2)
cax = divider.append_axes('right', size='5%', pad=0.05)
cb = fig.colorbar(pc_2,  cax=cax, orientation='vertical', label="Damage")
plt.tight_layout()

def animate(k):
    pc_1.set_array(damage_1[k].data.ravel())
    pc_2.set_array(damage_2[k].data.ravel())
    #ax1[0].title.set_text('iter = %i, time = %s'%(theta0.iter[k],np.datetime_as_string(theta0.time[k], unit='s')))

anim = animation.FuncAnimation(
    fig, animate, interval=100, frames=ds_1.iter.values.max()-1)
 
anim.save('/home/csys/mbourget/Desktop/Plots/constwind_2days/damage_compared.gif')