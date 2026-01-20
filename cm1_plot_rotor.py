# -*- coding: utf-8 -*-
"""
Created on Thu Jan 15 10:50:15 2026
Python version 3.11
@author: mschne28

Code to investigate the rotor. Will probably look at
1) Cross sections of model fields (horizontal vorticity, wind vectors, etc) around when parcels are in the rotor
    -Horizontal vorticity fields
    -Wind vectors
    -???
2) Parcel trajectories through the rotor (a la reviewer 3 comments)

--FILES NEEDED--
1) cm1out_000044?-000053.nc -- these are on bigbang under /merger/merger-125m
2) (maybe) cm1out_pdata.nc -- also on bigbang under /merger/merger-125m
3) traj_MV1.pkl -- Trajectory data of all parcels in the MV every 5 min. The pickles might all be on bigbang under /merger/merger-125m/pickles?
4) traj_clusters_220min_v2.pkl -- Indices of parcels in the MV at 220 min clustered by source region
5) storm_motion.pkl -- Estimated storm motion at every model output time (every 1 min) from my MV boxes
6) hvort_traj_220min.pkl -- Trajectory horizontal vorticity of parcels in the MV at 220 min interpolated from model output

"""

#%% Load modules

import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc
import pickle
from scipy.interpolate import interp1d
import pyart

#%% User inputs

# Path to wherever you're accessing the files from
fp = '/Volumes/Promise_Pegasus_70TB/merger/merger-125m/'
fp = 'C:/Users/mschne28/Documents/stuff_for_matt/'

plottime = 220 #Time to plot model fields (min) -- anywhere from probably 215-220 to look at the rotor development


#%% Load data

fnum = plottime - 167 # CM1 output file number
mvtime = 220


# ds = nc.Dataset(fp + f"cm1out_{fnum:06d}.nc")
# time = ds.variables['time'][:].data[0]
# xh = ds.variables['xh'][:].data
# yh = ds.variables['yh'][:].data
# zh = ds.variables['z'][:].data

# xvort = ds.variables['xvort'][:].data[0,:,:,:]
# yvort = ds.variables['yvort'][:].data[0,:,:,:]
# hvort = np.sqrt(xvort**2 + yvort**2)

# thrpert = ds.variables['thrpert'][:].data[:,:,:]

# uinterp = ds.variables['uinterp'][:].data[0,:,:,:]
# vinterp = ds.variables['vinterp'][:].data[0,:,:,:]
# winterp = ds.variables['winterp'][:].data[0,:,:,:]

# ds.close()


dbfile = open(fp + "coords.pkl", 'rb')
coords = pickle.load(dbfile)
xh = coords['xh']
yh = coords['yh']
zh = coords['zh']
dbfile.close()





dbfile = open(fp + "traj_MV1.pkl", 'rb')
traj = pickle.load(dbfile)
xp = traj[f"{mvtime}min"]['x']
yp = traj[f"{mvtime}min"]['y']
zp = traj[f"{mvtime}min"]['z']
up = traj[f"{mvtime}min"]['u']
vp = traj[f"{mvtime}min"]['v']
wp = traj[f"{mvtime}min"]['w']
dbfile.close()


dbfile = open(fp + f"traj_clusters_{mvtime}min_v2.pkl", 'rb')
dat = pickle.load(dbfile)
cc = dat['mv1']
x_ml = xp[:,(cc==1)]
y_ml = yp[:,(cc==1)]
z_ml = zp[:,(cc==1)]
u_ml = up[:,(cc==1)]
v_ml = vp[:,(cc==1)]
w_ml = wp[:,(cc==1)]
dbfile.close()


dbfile = open(fp + "storm_motion.pkl", 'rb')
dat = pickle.load(dbfile)
u_storm = dat['u_storm']
v_storm = dat['v_storm']
dbfile.close()


mtime = np.linspace(10800, 14400, 61) # CM1 output times
ptime = np.linspace(10800, 14400, 241) # Parcel output times

ti = np.where(ptime == mvtime*60)[0][0]


# Interpolate storm motion from model output timesto parcel output times (60 s to 15 s)
fu = interp1d(mtime, u_storm)
u_storm_prcl = fu(ptime)
fv = interp1d(mtime, v_storm)
v_storm_prcl = fv(ptime)


# Get storm-relative trajectories starting from MV
x_sr = np.zeros(shape=x_ml[:ti+1,:].shape, dtype=float)
y_sr = np.zeros(shape=y_ml[:ti+1,:].shape, dtype=float)
x_sr[ti,:] = x_ml[ti,:]
y_sr[ti,:] = y_ml[ti,:]
inds = np.linspace(ti-1, 0, ti)
for i in inds:
    i = int(i)
    dt = ptime[i+1] - ptime[i]
    delta_x = np.sum(u_storm_prcl[i:ti] * dt)
    delta_y = np.sum(v_storm_prcl[i:ti] * dt)
    x_sr[i,:] = x_ml[i] + delta_x
    y_sr[i,:] = y_ml[i] + delta_y




# Load saved parcel trajectory horizontal vorticity (interpolated from model fields)
dbfile = open(fp + f"hvort_traj_{mvtime}min.pkl", 'rb')
vort_traj = pickle.load(dbfile)
xvort_ml = vort_traj['xvort_ml']
yvort_ml = vort_traj['yvort_ml']
hvort_ml = np.sqrt(xvort_ml**2 + yvort_ml**2)
dbfile.close()



# Median parcel trajectory values
x_median = np.median(x_sr[:ti+1,:], axis=1)/1000
y_median = np.median(y_sr[:ti+1,:], axis=1)/1000
z_median = np.median(z_ml[:ti+1,:], axis=1)/1000
u_median = np.median(u_ml[:ti+1,:], axis=1)
v_median = np.median(v_ml[:ti+1,:], axis=1)
w_median = np.median(w_ml[:ti+1,:], axis=1)
xvort_median = np.median(xvort_ml[:ti+1,:], axis=1)
yvort_median = np.median(yvort_ml[:ti+1,:], axis=1)
hvort_median = np.median(hvort_ml[:ti+1,:], axis=1)


#%% Make plots to investigate rotor -> model output cross sections

ti = np.where(ptime == plottime*60)[0][0]


ix = np.argmin(np.abs(xh - x_median[ti]))
iy = np.argmin(np.abs(yh - y_median[ti]))

ds = nc.Dataset(fp + f"cm1out_{fnum:06d}.nc")
# X-Z cross section
thrpert_xz = ds.variables['thrpert'][:].data[:,iy,:]
u_xz = ds.variables['uinterp'][:].data[0,:,iy,:]
v_xz = ds.variables['vinterp'][:].data[0,:,iy,:]
w_xz = ds.variables['winterp'][:].data[0,:,iy,:]
xvort_xz = ds.variables['xvort'][:].data[0,:,iy,:]
yvort_xz = ds.variables['yvort'][:].data[0,:,iy,:]
hvort_xz = np.sqrt(xvort_xz**2 + yvort_xz**2)

# Y-Z cross section
thrpert_yz = ds.variables['thrpert'][:].data[:,:,ix]
u_yz = ds.variables['uinterp'][:].data[0,:,:,ix]
v_yz = ds.variables['vinterp'][:].data[0,:,:,ix]
w_yz = ds.variables['winterp'][:].data[0,:,:,ix]
xvort_yz = ds.variables['xvort'][:].data[0,:,:,ix]
yvort_yz = ds.variables['yvort'][:].data[0,:,:,ix]
hvort_yz = np.sqrt(xvort_yz**2 + yvort_yz**2)
ds.close()


levs = np.linspace(0, 0.02, 21) #contour levels for hvort

# Plot cross sections

# X-Z cross section -- Horizontal vorticity magnitude
fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
c = ax.contourf(xh, zh, hvort_xz, levels=levs, vmin=levs[0], vmax=levs[-1], cmap='Reds', antialiased=True)
c.set_edgecolor('face')
ax.quiver(xh[::4], zh[::3], u_xz[::3,::4]+6, w_xz[::3,::4], color='k', scale=500, width=0.003, pivot='tail')
ax.scatter(x_median[:ti+1], z_median[:ti+1], s=30, c='k') #parcel median SR trajectory
ax.set_xlim(-15, 0)
ax.set_ylim(0, 3)
ax.set_title(f"\u03c9$_H$ X-Z cross section, {plottime} min", fontsize=14)
ax.set_xlabel('x (km)', fontsize=12)
ax.set_ylabel('z (km)', fontsize=12)
cb = plt.colorbar(c, ax=ax, extend='max')
cb.set_label("\u03c9$_H$ (s$^{-1}$)", fontsize=12)
cb.formatter.set_powerlimits((0,0))
plt.show()
# plt.savefig(fp + f"hvort_xz_cross_section_{plottime}min.png", dpi=300)



# Y-Z cross section -- Horizontal vorticity magnitude
fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
c = ax.contourf(yh, zh, hvort_yz, levels=levs, vmin=levs[0], vmax=levs[-1], cmap='Reds', antialiased=True)
c.set_edgecolor('face')
# ax.quiver(yh[::4], zh[::3], v_yz[::3,::4], w_yz[::3,::4], color='k', scale=500, width=0.003, pivot='tail')
ax.scatter(y_median[:ti+1], z_median[:ti+1], s=30, c='k') #parcel median SR trajectory
ax.set_xlim(-60, -45)
ax.set_ylim(0, 3)
ax.set_title(f"\u03c9$_H$ Y-Z cross section, {plottime} min", fontsize=14)
ax.set_xlabel('y (km)', fontsize=12)
ax.set_ylabel('z (km)', fontsize=12)
cb = plt.colorbar(c, ax=ax, extend='max')
cb.set_label("\u03c9$_H$ (s$^{-1}$)", fontsize=12)
cb.formatter.set_powerlimits((0,0))
plt.show()
# plt.savefig(fp + f"hvort_yz_cross_section_{plottime}min.png", dpi=300)




#%% Combined cross sections

from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1 import make_axes_locatable


ti = np.where(ptime == plottime*60)[0][0]


ix = np.argmin(np.abs(xh - x_median[ti]))
iy = np.argmin(np.abs(yh - y_median[ti]))
iz = np.argmin(np.abs(zh - z_median[ti]))


ds = nc.Dataset(fp + f"cm1out_{fnum:06d}.nc")
# Plan view - closest to median parcel height
thrpert = ds.variables['thrpert'][:].data[iz,:,:]
uinterp = ds.variables['uinterp'][:].data[0,iz,:,:]
vinterp = ds.variables['vinterp'][:].data[0,iz,:,:]
xvort = ds.variables['xvort'][:].data[0,iz,:,:]
yvort = ds.variables['yvort'][:].data[0,iz,:,:]
hvort = np.sqrt(xvort**2 + yvort**2)

# X-Z cross section
thrpert_xz = ds.variables['thrpert'][:].data[:,iy,:]
u_xz = ds.variables['uinterp'][:].data[0,:,iy,:]
w_xz = ds.variables['winterp'][:].data[0,:,iy,:]
xvort_xz = ds.variables['xvort'][:].data[0,:,iy,:]
yvort_xz = ds.variables['yvort'][:].data[0,:,iy,:]
hvort_xz = np.sqrt(xvort_xz**2 + yvort_xz**2)

# Y-Z cross section
thrpert_yz = ds.variables['thrpert'][:].data[:,:,ix]
v_yz = ds.variables['vinterp'][:].data[0,:,:,ix]
w_yz = ds.variables['winterp'][:].data[0,:,:,ix]
xvort_yz = ds.variables['xvort'][:].data[0,:,:,ix]
yvort_yz = ds.variables['yvort'][:].data[0,:,:,ix]
hvort_yz = np.sqrt(xvort_yz**2 + yvort_yz**2)
ds.close()



xlim = [-15, 0]
ylim = [-60, -45]
zlim = [0, 3]



fig,ax = plt.subplots(figsize=(11,8), layout='constrained')

divider = make_axes_locatable(ax)
top_ax = divider.append_axes('top', 1.85, pad=0.2, sharex=ax)
side_ax = divider.append_axes('left', 2, pad=0.2, sharey=ax)

top_ax.xaxis.set_tick_params(labelbottom=False)
ax.yaxis.set_tick_params(labelleft=False)

c = ax.contourf(xh, yh, hvort, levels=levs, vmin=levs[0], vmax=levs[-1], cmap='Reds', antialiased=True)
c.set_edgecolor('face')
# ax.quiver(xh[::4], yh[::4], uinterp[::4,::4]+6, vinterp[::4,::4], color='k', scale=500, width=0.003, pivot='tail')
ax.scatter(x_median[:ti+1], y_median[:ti+1], s=30, c='k') #parcel median SR trajectory
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_xlabel('x (km)', fontsize=12)
ax.axhline(yh[iy], color='k', linewidth=0.75, linestyle='--')
ax.axvline(xh[ix], color='k', linewidth=0.75, linestyle='--')
cb = plt.colorbar(c, ax=ax, extend='max')
cb.set_label("\u03c9$_H$ (s$^{-1}$)", fontsize=12)
cb.formatter.set_powerlimits((0,0))

c1 = top_ax.contourf(xh, zh, hvort_xz, levels=levs, vmin=levs[0], vmax=levs[-1], cmap='Reds', antialiased=True)
c1.set_edgecolor('face')
# top_ax.quiver(xh[::4], zh[::4], u_xz[::4,::4], w_xz[::4,::4], color='k', scale=500, width=0.003, pivot='tail')
top_ax.scatter(x_median[:ti+1], z_median[:ti+1], s=30, c='k')
top_ax.set_xlim(xlim)
top_ax.set_ylim(zlim)
top_ax.set_ylabel('z (km)', fontsize=12)

c2 = side_ax.contourf(zh, yh, hvort_yz.transpose(), levels=levs, vmin=levs[0], vmax=levs[-1], cmap='Reds', antialiased=True)
c2.set_edgecolor('face')
# side_ax.quiver(zh[::4], yh[::4], w_yz[::4,::4].transpose(), v_yz[::4,::4].transpose(), color='k', scale=500, width=0.003, pivot='tail')
side_ax.scatter(z_median[:ti+1], y_median[:ti+1], s=30, c='k')
side_ax.set_xlim(zlim)
side_ax.set_ylim(ylim)
side_ax.set_xlabel('z (km)', fontsize=12)
side_ax.set_ylabel('y (km)', fontsize=12)
side_ax.invert_xaxis()

plt.suptitle(f"\u03c9$_H$ cross sections, {plottime} min", fontsize=12)
plt.show()
# plt.savefig(fp + f"hvort_cross_sections_{plottime}min.png", dpi=300)









