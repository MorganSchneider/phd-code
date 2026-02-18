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
2) Parcel trajectories through the rotor?

--FILES NEEDED--
1) cm1out_000044?-000053 or 58.nc -- these are on bigbang under /merger/merger-125m
2) cm1out_000013 -- also on bigbang under /merger/merger-125m. Needed for prspert
3) traj_MV1.pkl -- Trajectory data of all parcels in the MV every 5 min. The pickles might all be on bigbang under /merger/merger-125m/pickles?
4) traj_clusters_210min_v2.pkl, traj_clusters_220min_v2.pkl -- Indices of parcels in the MV at 210/220 min clustered by source region
5) storm_motion.pkl -- Estimated storm motion at every model output time (every 1 min) from my MV boxes

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

plottime = 220 #Time to plot model fields (min) -- anywhere from probably 215-220 or 225 to look at the rotor development


#%% Load parcel data

fnum = plottime - 167 # CM1 output file number


ds = nc.Dataset(fp + f"cm1out_{fnum:06d}.nc")
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
zh = ds.variables['z'][:].data
ds.close()



dbfile = open(fp + "traj_MV1.pkl", 'rb')
traj = pickle.load(dbfile)
xp = traj["220min"]['x']
yp = traj["220min"]['y']
zp = traj["220min"]['z']
up = traj["220min"]['u']
vp = traj["220min"]['v']
wp = traj["220min"]['w']
dbfile.close()


dbfile = open(fp + "traj_clusters_220min_v2.pkl", 'rb')
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

ti = np.where(ptime == 220*60)[0][0]


# Interpolate storm motion from model output times to parcel output times (60 s to 15 s)
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
dbfile = open(fp + "hvort_traj_220min.pkl", 'rb')
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

# x/y indices of median parcel trajectory position at t=plottime
ix = np.argmin(np.abs(xh - x_median[ti]))
iy = np.argmin(np.abs(yh - y_median[ti]))

# Load data
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


levs = np.linspace(0, 0.02, 21) #contour levels for hvort - change if you plot a different variable

# plot limits - change as you want
xlim = [-15, 0]
ylim = [-60, -45]
zlim = [0, 3]


# X-Z cross section -- Horizontal vorticity magnitude
fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
c = ax.contourf(xh, zh, hvort_xz, levels=levs, vmin=levs[0], vmax=levs[-1], cmap='Reds', antialiased=True)
c.set_edgecolor('face')
ax.quiver(xh[::4], zh[::3], u_xz[::3,::4]+6, w_xz[::3,::4], color='k', scale=500, width=0.003, pivot='tail')
ax.scatter(x_median[:ti+1], z_median[:ti+1], s=30, c='k') #parcel median SR trajectory
ax.set_xlim(xlim)
ax.set_ylim(zlim)
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
ax.set_xlim(ylim)
ax.set_ylim(zlim)
ax.set_title(f"\u03c9$_H$ Y-Z cross section, {plottime} min", fontsize=14)
ax.set_xlabel('y (km)', fontsize=12)
ax.set_ylabel('z (km)', fontsize=12)
cb = plt.colorbar(c, ax=ax, extend='max')
cb.set_label("\u03c9$_H$ (s$^{-1}$)", fontsize=12)
cb.formatter.set_powerlimits((0,0))
plt.show()
# plt.savefig(fp + f"hvort_yz_cross_section_{plottime}min.png", dpi=300)





#%% NEW SECTION: Remake Fig 18 cross section figures for 210 min and 220 min

import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc
import pickle
import pyart
from scipy.interpolate import interp1d
from matplotlib.ticker import MultipleLocator

###################
### USER INPUTS ###
###################
fp = '/Volumes/Promise_Pegasus_70TB/merger/merger-125m/'

mvtime = 220

if mvtime == 210:
    fn = 41
elif mvtime == 220:
    fn = 51

###################
### Everything else below should be automated... unless the wind vector spacing needs to be adjusted, in which case change the values for qix and qiz
###################


### Load data ###

ptime = np.linspace(10800, 14400, 241) # Parcel times

ti = np.where(ptime == mvtime*60)[0][0]
ti0 = np.where(ptime == 180*60)[0][0]

# Load MV parcel data
dbfile = open(fp + 'traj_MV1.pkl', 'rb')
traj = pickle.load(dbfile)
pids_mv = traj[f"{mvtime}min"]['pids']
x_mv = traj[f"{mvtime}min"]['x']
y_mv = traj[f"{mvtime}min"]['y']
z_mv = traj[f"{mvtime}min"]['z']
u_mv = traj[f"{mvtime}min"]['u']
v_mv = traj[f"{mvtime}min"]['v']
w_mv = traj[f"{mvtime}min"]['w']
dbfile.close()

# Load parcel source region clusters
dbfile = open(fp + f"traj_clusters_{mvtime}min_v2.pkl", 'rb')
c = pickle.load(dbfile)
cc = c['mv1']
dbfile.close()


# Mid-level source (0 -> LL inflow, 1 -> ML inflow, 2 -> supercell outflow, 3 -> QLCS outflow)
x_ml = x_mv[:,(cc==1)]
y_ml = y_mv[:,(cc==1)]
z_ml = z_mv[:,(cc==1)]
u_ml = u_mv[:,(cc==1)]
v_ml = v_mv[:,(cc==1)]
w_ml = w_mv[:,(cc==1)]


# Load saved storm motion estimates
dbfile = open(fp + 'storm_motion.pkl', 'rb')
sm = pickle.load(dbfile)
u_storm = sm['u_storm']
v_storm = sm['v_storm']
dbfile.close()

stimes = np.linspace(10800, 14400, 61)


# Interpolate storm motion from model output timesto parcel output times (60 s to 15 s)
fu = interp1d(stimes, u_storm)
u_storm_interp = fu(ptime)
fv = interp1d(stimes, v_storm)
v_storm_interp = fv(ptime)


# Get storm-relative trajectories starting from MV
x_sr = np.zeros(shape=x_ml[:ti+1,:].shape, dtype=float)
y_sr = np.zeros(shape=y_ml[:ti+1,:].shape, dtype=float)
x_sr[ti,:] = x_ml[ti,:]
y_sr[ti,:] = y_ml[ti,:]
inds = np.linspace(ti-1, 0, ti)
for i in inds:
    i = int(i)
    dt = ptime[i+1] - ptime[i]
    delta_x = np.sum(u_storm_interp[i:ti] * dt)
    delta_y = np.sum(v_storm_interp[i:ti] * dt)
    # delta_x = u_storm_interp[i] * dt
    # delta_y = v_storm_interp[i] * dt
    x_sr[i,:] = x_ml[i] + delta_x
    y_sr[i,:] = y_ml[i] + delta_y


# Median parcel trajectory
x_median = np.median(x_sr[ti0:ti+1,:], axis=1)/1000
y_median = np.median(y_sr[ti0:ti+1,:], axis=1)/1000
z_median = np.median(z_ml[ti0:ti+1,:], axis=1)/1000
u_median = np.median(u_ml[ti0:ti+1,:], axis=1)
v_median = np.median(v_ml[ti0:ti+1,:], axis=1)
w_median = np.median(w_ml[ti0:ti+1,:], axis=1)


# Set x/z limits to reduce memory needed for loading variables
xl = [-30,15] #[-55,25]
zl = [0,3.5]

# Load model fields
ds = nc.Dataset(fp + f"cm1out_{fn:06d}.nc")
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
zh = ds.variables['z'][:].data
iz = slice(0, np.where(zh >= zl[1])[0][1])
ix = slice(np.where(xh >= xl[0])[0][0], np.where(xh >= xl[1])[0][1])
iy = (np.abs(yh - np.median(y_ml[ti-9,:]/1000))).argmin()

thrpert = ds.variables['thrpert'][:].data[iz,iy,ix]
prs = ds.variables['prs'][:].data[0,iz,iy,ix]
uinterp = ds.variables['uinterp'][:].data[0,iz,iy,ix]
winterp = ds.variables['winterp'][:].data[0,iz,iy,ix]
ds.close()

# Get base state pressure to calculate prspert
ds = nc.Dataset(fp+'cm1out_000013.nc')
prspert = (prs - ds.variables['prs0'][:].data[0,iz,iy,ix])/100
ds.close()


# Custom colormap for the parcel w colorbar because I'm super normal about colorbars
from matplotlib.colors import LinearSegmentedColormap
import matplotlib as mpl
parcel_cm = mpl.colors.LinearSegmentedColormap.from_list('parcel_cm',
                    np.vstack((pyart.graph.cmweather.cm_colorblind.ChaseSpectral(np.linspace(0.1,0.5,154)),
                    pyart.graph.cmweather.cm_colorblind.ChaseSpectral(np.linspace(0.5,1,102)))))



##################
### MAKE PLOTS ###
##################

figsave = False


wl = [-15,10] # Parcel w datalims

qix = 4 # x spacing of wind vectors
qiz = 3 # z spacing of wind vectors


### End of user inputs

fig,ax = plt.subplots(1,1,figsize=(10,4), layout='constrained')

c = ax.contourf(xh[ix], zh[iz], thrpert, levels=np.linspace(-16,16,17), vmin=[-16], vmax=[16], cmap='PuOr_r', antialiased=True)
c.set_edgecolor('face')
ax.contour(xh[ix], zh[iz], prspert, levels=[-3,-2,-1,1,2,3], colors=['b','b','b','r','r','r'], linewidths=0.75, linestyles='-')
ax.quiver(xh[ix][::qix], zh[iz][::qiz], uinterp[::qiz,::qix], winterp[::qiz,::qix], scale=250, width=0.002, pivot='tail')
p = ax.scatter(x_median, z_median, s=60, c=w_median, cmap=parcel_cm, edgecolor='none', linewidth=0.25, vmin=wl[0], vmax=wl[1])
ax.scatter(x_median[::20], z_median[::20], s=60, c=w_median[::20], cmap=parcel_cm, edgecolor='k', linewidth=1.5, vmin=wl[0], vmax=wl[1])
cb = plt.colorbar(c, ax=ax, extend='both')
cb.set_label("\u03B8'\u1D68 (K)", fontsize=13)
cb2 = plt.colorbar(p, ax=ax, extend='both') #, ticks=np.linspace(wl[0],wl[1],11))
cb2.set_label("Parcel $w$ (m s$^{-1}$)", fontsize=13)
for i in range(len(x_median[::20])):
    t = ptime[ti0:ti+1][::20][i]/60
    if mvtime == 210:
        ax.text(x_median[::20][i]+0.3, z_median[::20][i]+0.05, f"{t:.0f}", fontsize=12, fontweight='bold')
    elif mvtime == 220:
        if (i==1) | (i==3):
            x_offset = -1; z_offset = -0.25
        elif (i == 6):
            x_offset = -0.2; z_offset = 0.1
        elif i >= 7:
            x_offset = -2.25; z_offset = 0.1
        else:
            x_offset = -1; z_offset = 0.1
        ax.text(x_median[::20][i]+x_offset, z_median[::20][i]+z_offset, f"{t:.0f}", fontsize=12, fontweight='bold')
ax.set_xlabel('x (km)', fontsize=14)
ax.set_ylabel('z (km)', fontsize=14)
ax.set_xlim([-25,0])
ax.set_ylim([0,3.5])
# ax.set_title(f"Mid-level parcels in the MV at {mv_time} min", fontsize=14)

if figsave:
    plt.savefig(fp + f"traj_thr_cross_{mvtime}min_SR.png", dpi=300)



#%% NEW SECTION: Make new cross section figure of y-vorticity -- probably do panels for 217 and 218 min


import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc
import pickle
import pyart
from scipy.interpolate import interp1d
from matplotlib.ticker import MultipleLocator

###################
### USER INPUTS ###
###################

fp = '/Volumes/Promise_Pegasus_70TB/merger/merger-125m/'

plottime = 217 #do 217, 218, and 219

###################
### Everything else below should be automated...
### ... unless the wind vector spacing, contour levels/limits, or axis limits need adjusting
###################


### Load data

mvtime = 220
fn = plottime - 167


ptime = np.linspace(10800, 14400, 241) # set parcel times

ti = np.where(ptime == mvtime*60)[0][0]
ti0 = np.where(ptime == 180*60)[0][0]

# Load MV parcel data
dbfile = open(fp + 'traj_MV1.pkl', 'rb')
traj = pickle.load(dbfile)
pids_mv = traj[f"{mvtime}min"]['pids']
x_mv = traj[f"{mvtime}min"]['x']
y_mv = traj[f"{mvtime}min"]['y']
z_mv = traj[f"{mvtime}min"]['z']
u_mv = traj[f"{mvtime}min"]['u']
v_mv = traj[f"{mvtime}min"]['v']
w_mv = traj[f"{mvtime}min"]['w']
dbfile.close()

# Load parcel source region clusters
dbfile = open(fp + f"traj_clusters_{mvtime}min_v2.pkl", 'rb')
c = pickle.load(dbfile)
cc = c['mv1']
dbfile.close()


# Mid-level source (0 -> LL inflow, 1 -> ML inflow, 2 -> supercell outflow, 3 -> QLCS outflow)
x_ml = x_mv[:,(cc==1)]
y_ml = y_mv[:,(cc==1)]
z_ml = z_mv[:,(cc==1)]
u_ml = u_mv[:,(cc==1)]
v_ml = v_mv[:,(cc==1)]
w_ml = w_mv[:,(cc==1)]


# Load saved storm motion estimates
dbfile = open(fp + 'storm_motion.pkl', 'rb')
sm = pickle.load(dbfile)
u_storm = sm['u_storm']
v_storm = sm['v_storm']
dbfile.close()

stimes = np.linspace(10800, 14400, 61)


# Interpolate storm motion from model output timesto parcel output times (60 s to 15 s)
fu = interp1d(stimes, u_storm)
u_storm_interp = fu(ptime)
fv = interp1d(stimes, v_storm)
v_storm_interp = fv(ptime)


# Get storm-relative trajectories starting from MV
x_sr = np.zeros(shape=x_ml[:ti+1,:].shape, dtype=float)
y_sr = np.zeros(shape=y_ml[:ti+1,:].shape, dtype=float)
x_sr[ti,:] = x_ml[ti,:]
y_sr[ti,:] = y_ml[ti,:]
inds = np.linspace(ti-1, 0, ti)
for i in inds:
    i = int(i)
    dt = ptime[i+1] - ptime[i]
    delta_x = np.sum(u_storm_interp[i:ti] * dt)
    delta_y = np.sum(v_storm_interp[i:ti] * dt)
    # delta_x = u_storm_interp[i] * dt
    # delta_y = v_storm_interp[i] * dt
    x_sr[i,:] = x_ml[i] + delta_x
    y_sr[i,:] = y_ml[i] + delta_y


# Median parcel trajectory
x_median = np.median(x_sr[ti0:ti+1,:], axis=1)/1000
y_median = np.median(y_sr[ti0:ti+1,:], axis=1)/1000
z_median = np.median(z_ml[ti0:ti+1,:], axis=1)/1000
u_median = np.median(u_ml[ti0:ti+1,:], axis=1)
v_median = np.median(v_ml[ti0:ti+1,:], axis=1)
w_median = np.median(w_ml[ti0:ti+1,:], axis=1)


# Set x and z limits to reduce memory needed to load data
xl = [-30,15] #[-55,25]
zl = [0,3.5]

ti = np.where(ptime == plottime*60)[0][0] # index for parcel positions at t=plottime

# Load model fields
ds = nc.Dataset(fp + f"cm1out_{fn:06d}.nc")
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
zh = ds.variables['z'][:].data

ix = slice(np.where(xh >= xl[0])[0][0], np.where(xh >= xl[1])[0][1])
iy = np.argmin(np.abs(yh - y_median[ti]))
iz = slice(0, np.where(zh >= zl[1])[0][1])

yvort = ds.variables['yvort'][:].data[0,iz,iy,ix]
uinterp = ds.variables['uinterp'][:].data[0,iz,iy,ix]
winterp = ds.variables['winterp'][:].data[0,iz,iy,ix]
ds.close()



##################
### MAKE PLOTS ###
##################

# Horizontal and vertical spacing of wind vectors
qix = 4
qiz = 3




fig,ax = plt.subplots(1, 1, figsize=(9,4), layout='constrained')

c = ax.contourf(xh[ix], zh[iz], yvort, levels=np.linspace(-0.1, 0.1, 21), vmin=-0.1, vmax=0.1, cmap='balance', antialiased=True)
c.set_edgecolor('face')
ax.quiver(xh[ix][::qix], zh[iz][::qiz], uinterp[::qiz,::qix], winterp[::qiz,::qix], scale=250, width=0.002, pivot='tail')
p = ax.scatter(x_median, z_median, s=60, edgecolor='k', facecolor='w', linewidth=0.25)
ax.scatter(x_median[::20], z_median[::20], edgecolor='k', facecolor='w', linewidth=0.25)
cb = plt.colorbar(c, ax=ax, extend='both', ticks=np.linspace(-0.1,0.1,11))
cb.set_label("\u03B7 (s$^{-1}$)", fontsize=13)
for i in range(len(x_median[::20])):
    t = ptime[ti0:ti+1][::20][i]/60
    if mvtime == 210:
        ax.text(x_median[::20][i]+0.3, z_median[::20][i]+0.05, f"{t:.0f}", fontsize=12, fontweight='bold')
    elif mvtime == 220:
        if (i==1) | (i==3):
            x_offset = -1; z_offset = -0.25
        elif (i == 6):
            x_offset = -0.2; z_offset = 0.1
        elif i >= 7:
            x_offset = -2.25; z_offset = 0.1
        else:
            x_offset = -1; z_offset = 0.1
        ax.text(x_median[::20][i]+x_offset, z_median[::20][i]+z_offset, f"{t:.0f}", fontsize=12, fontweight='bold')
ax.set_xlabel('x (km)', fontsize=14)
ax.set_ylabel('z (km)', fontsize=14)
ax.set_xlim([-25,0])
ax.set_ylim([0,3.5])
# ax.set_title(f"Mid-level parcels in the MV at {mv_time} min", fontsize=14)

if figsave:
    plt.savefig(fp + f"traj_yvort_cross_{plottime}min.png", dpi=300)




