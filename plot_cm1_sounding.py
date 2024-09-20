# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 12:52:41 2022

@author: Tony Reinhart

Plots CM1 base state sounding
"""

### Load modules ###

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import netCDF4 as nc
import metpy as met
from metpy.plots import SkewT,Hodograph
import metpy.calc as mc
from metpy.units import units
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#%% Read CM1 output

fp = 'D:/cm1out/merger/merger2-500m/'
#fp = '/Users/morgan.schneider/Documents/cm1/merger/ctl-new/'
outfile = fp+'cm1out_000001.nc'

umove = 6.
vmove = 0.

ds = nc.Dataset(outfile)
ix = np.where(ds.variables['xf'][:].data >= 20)[0][0]
iy = np.where(ds.variables['yf'][:].data >= 50)[0][0]

cm1sh = ds.variables['z'][:].data
cm1th = ds.variables['th'][:].data[0,:,iy,ix]
cm1sp = ds.variables['prs'][:].data[0,:,iy,ix]
cm1sq = ds.variables['qv'][:].data[0,:,iy,ix]
cm1su = ds.variables['uinterp'][:].data[0,:,iy,ix] + umove
cm1sv = ds.variables['vinterp'][:].data[0,:,iy,ix] + vmove
cm1st = np.array([cm1th[i]*(cm1sp[i]/100000.)**0.286 for i in range(len(cm1th))])

cm1th0 = ds.variables['th0'][:].data[0,:,iy,ix]
cm1sp0 = ds.variables['prs0'][:].data[0,:,iy,ix]
cm1sq0 = ds.variables['qv0'][:].data[0,:,iy,ix]
cm1su0 = ds.variables['u0'][:].data[0,:,iy,ix] + umove
cm1sv0 = ds.variables['v0'][:].data[0,:,iy,ix] + vmove
cm1st0 = np.array([cm1th0[i]*(cm1sp0[i]/100000.)**0.286 for i in range(len(cm1th0))])

ds.close()

#%%

#sndfile = fp+'input_sounding'
sndfp = '/Users/morgan.schneider/Documents/PERiLS_LIDAR_Soundings/IOP2/'
sndfile = sndfp+'NSSL1_MW41_output_20220330_195946.csv'

df = pd.read_csv(sndfile,header=2)

sndst = np.array(list(df['Filtered Temperature (K)']))
sndtd = np.array(list(df['Filtered Dewpoint (K)']))
sndsp = np.array(list(df['Filtered Pressure (mb)']))
sndwd = np.array(list(df['Filtered Wind Dir']))
sndws = np.array(list(df['Filtered Wind Spd (m/s)']))
sndsh = np.array(list(df['Filtered Altitude (m)']))

wnd = mc.wind_components(sndws*units('m/s'), sndwd*units.deg)
sndsu = wnd[0].magnitude[:]
sndsv = wnd[1].magnitude[:]


# fp = '/Users/morgan.schneider/Documents/python/cm1_skewt_plots/'
# sndsh = np.genfromtxt(fp+'cm1soundingSynthobHeights',skip_header=2,usecols=1)
# sndst = np.genfromtxt(fp+'cm1soundingSynthobThermo',skip_header=2,usecols=3)
# sndsp = np.genfromtxt(fp+'cm1soundingSynthobPres',skip_header=2,usecols=1)
# sndsq = np.genfromtxt(fp+'cm1soundingSynthobMoist',skip_header=2,usecols=3)
# sndsu = np.genfromtxt(fp+'cm1soundingSynthobWind',skip_header=2,usecols=1)
# sndsv = np.genfromtxt(fp+'cm1soundingSynthobWind',skip_header=2,usecols=2)

#%%

e = (cm1sq*cm1sp/100) / (0.622+cm1sq)
dewpt = 243.5 / ((17.67/(np.log(e/6.112)))-1)

e0 = (cm1sq0*cm1sp0/100) / (0.622+cm1sq0)
dewpt0 = 243.5 / ((17.67/(np.log(e0/6.112)))-1)

#%% Plot the skew-T

#plt.close('all')

### Plot CM1 point profile ###
plt.rcParams['figure.figsize']=(9,9)

skew = SkewT()
skew.plot(cm1sp/100., dewpt, '-g', linewidth=3)
skew.plot(cm1sp/100., (cm1st-273.15), '-r', linewidth=3)
# Add the relevant special lines
skew.plot_dry_adiabats()
skew.plot_moist_adiabats()
skew.plot_mixing_lines()
skew.ax.set_ylim(1000, 100)
skew.ax.set_xlim(-40, 30)
plt.title('CM1 point profile')
# Create a hodograph
ax_hod = inset_axes(skew.ax, '42%', '42%', loc=1)
h = Hodograph(ax_hod, component_range=40.)
h.add_grid(increment=20)
h.plot(cm1su, cm1sv, color='k', linewidth=1.8)

plt.show()
#plt.savefig(fp+'cm1_profile.png',dpi=300)


### Plot CM1 base state ###
plt.rcParams['figure.figsize']=(9,9)

skew = SkewT()
skew.plot(cm1sp0/100., dewpt0, '-g', linewidth=3)
skew.plot(cm1sp0/100., (cm1st0-273.15), '-r', linewidth=3)
# Add the relevant special lines
skew.plot_dry_adiabats()
skew.plot_moist_adiabats()
skew.plot_mixing_lines()
skew.ax.set_ylim(1000, 100)
skew.ax.set_xlim(-40, 30)
plt.title('CM1 base state')
# Create a hodograph
ax_hod = inset_axes(skew.ax, '42%', '42%', loc=1)
h = Hodograph(ax_hod, component_range=40.)
h.add_grid(increment=20)
h.plot(cm1su0, cm1sv0, color='k', linewidth=1.8)

plt.show()
#plt.savefig(fp+'cm1_profile.png',dpi=300)


'''
### Plot raw input sounding ###
plt.rcParams['figure.figsize']=(9,9)

skew2 = SkewT()
skew2.plot(cm1sp/100., dewpt, '-g', linewidth=3)
skew2.plot(cm1sp/100., (cm1st-273.15), '-r', linewidth=3)
skew2.plot(sndsp, (sndtd-273.15), '-k', linewidth=1.5)
skew2.plot(sndsp, (sndst-273.15), '-k', linewidth=1.5)
# Add special lines
skew2.plot_dry_adiabats()
skew2.plot_moist_adiabats()
skew2.plot_mixing_lines()
skew2.ax.set_ylim(1000, 100)
skew2.ax.set_xlim(-40, 30)
plt.title('CM1 base state and input sounding')
# Create hodograph
ax_hod = inset_axes(skew2.ax, '42%', '42%', loc=1)
h2 = Hodograph(ax_hod, component_range=40.)
h2.add_grid(increment=20)
h2.plot(cm1su, cm1sv, color='m', linewidth=1.8)
h2.plot(sndsu, sndsv, color='k', linewidth=1)
h2.plot(cm1su-umove, cm1sv-vmove, color='C5', linewidth=1.5)

plt.show()
#plt.savefig(fp+'cm1+snd_profile.png',dpi=300)
'''










