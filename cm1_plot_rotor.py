# -*- coding: utf-8 -*-
"""
Created on Thu Jan 15 10:50:15 2026

@author: mschne28
"""

#%% Load modules

import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc
import pickle
from scipy.interpolate import interp1d

#%% User inputs

# Path to wherever you're accessing the files from
fp = '/Volumes/Promise_Pegasus_70TB/merger/merger-125m/'

# Time to plot things
mvtime = 220


#%% Load data

fnum = mvtime - 167 # CM1 file number


ds = nc.Dataset(fp + "cm1out_{fnum:06d}.nc")
time = ds.variables['time'][:].data[0]
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
zh = ds.variables['z'][:].data

xvort = ds.variables['xvort'][:].data[0,:,:,:]
yvort = ds.variables['yvort'][:].data[0,:,:,:]
hvort = np.sqrt(xvort**2 + yvort**2)

thrpert = ds.variables['thrpert'][:].data[:,:,:]

uinterp = ds.variables['uinterp'][:].data[0,:,:,:]
vinterp = ds.variables['vinterp'][:].data[0,:,:,:]
winterp = ds.variables['winterp'][:].data[0,:,:,:]

ds.close()





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











