# -*- coding: utf-8 -*-
"""
Created on Fri Sep 26 13:14:53 2025

@author: mschne28
"""

# on bigbang, I *think* my output files are in /raid/students/mschneider/merger/merger-125m
# I'll send the pickle files for parcel positions etc.

import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc
import cmocean
import pyart #need an earlier version of xarray -> 0.20.2 or earlier
import pickle

#%% Set user-defined things

# Set file paths here
fp_sim = '/path/to/cm1/output/files/'
fp_pkl = '/path/to/pickle/files/' #I will send you these files

fp_pkl = '/Users/mschne28/Documents/merger/merger-125m/'

# This is the time we analyze the mesocyclone at
# Run this for 210 and 220
time = 210

plot_xz = True #plot X-Z cross sections of dbz?
plot_yz = False #plot Y-Z cross sections of dbz?



#%% Load the data and make the plots


##### Load data #####

# i saved the model domain xh/yh/zh and parcel time to their own file
dbfile = open(fp_pkl+'coords.pkl', 'rb')
coords = pickle.load(dbfile)
xh = coords['xh'] #model grid x (in km)
yh = coords['yh'] #model grid y
zh = coords['zh'] #model grid z
ptime = coords['ptime'] #parcel times - every 15s from 10800s (180 min) to 14400s (240 min)
dbfile.close()


# This is basically an array of the source regions for mesocyclone parcels
# E.g. wherever cc=0, the parcel ID at that index is in the low-level inflow source
# cc=1 is mid-level inflow, cc=2 is supercell outflow, cc=3 is qlcs outflow
dbfile = open(fp_pkl+f"traj_clusters_{time}min_v2.pkl", 'rb')
tmp = pickle.load(dbfile)
cc = tmp['mv1'] #vector of source region assignments
dbfile.close()


# This is the parcel data for only parcels in the mesocyclone at various times
# saved variables are parcel x, y, z, w, b, zvort, u, v, and pid
# I am also going to pull data for only the mid-level source parcels (wherever cc=1)
dbfile = open(fp_pkl+f"traj_MV1.pkl", 'rb')
traj = pickle.load(dbfile)
xp = traj[f"{time}min"]['x'][:,(cc==1)] #mid-level source parcel x (in m)
yp = traj[f"{time}min"]['y'][:,(cc==1)] #mid-level source parcel y
zp = traj[f"{time}min"]['z'][:,(cc==1)] #mid-level source parcel z
dbfile.close()

# Get median mid-level parcel positions
xp_ml = np.median(xp, axis=1)
yp_ml = np.median(yp, axis=1)
zp_ml = np.median(zp, axis=1)



# CM1 output file numbers for times when parcels are in the downdraft
if time == 210:
    # times = np.arange(203,209) #want to look at 203-208 min for the 210 min meso
    fnums = np.arange(36,42) #output files for 203-208 min
elif time == 220:
    # times = np.arange(214,220) #want to look at 214-219 min for the 220 min meso
    fnums = np.arange(47,53) #output files for 214-219 min


# Setting x limits to load less data and just load within the limits I'll plot in
xl = [-40, 20] #x limits for loading less data
yl = [-90, -30] #y limits for loading less data
zl = [0, 4] # z limits
jx = slice(np.argmin(np.abs(xh - xl[0])), np.argmin(np.abs(xh - xl[1]))) #get indices
jy = slice(np.argmin(np.abs(yh - yl[0])), np.argmin(np.abs(yh - yl[1])))
jz = slice(0, np.argmin(zh - zl[1]))

# Preset arrays for dbz cross sections
dbz_xz = np.zeros(shape=(len(fnums), len(zh[jz]), len(xh[jx])), dtype=float)
dbz_yz = np.zeros(shape=(len(fnums), len(zh[jz]), len(yh[jy])), dtype=float)





##### Loop through downdraft times, load data, and make plots #####
for i in range(len(fnums)):
    ds = nc.Dataset(fp_sim+f"cm1out_{fnums[i]:06d}.nc")
    stime = ds.variables['time'][:].data[0] #model time
    
    # Index finding
    ti = np.where(ptime == stime)[0][0] #index where parcel time = model output time
    xi = np.argmin(np.abs(xh*1000 - xp_ml[ti])) #index where model x = median parcel x at stime
    yi = np.argmin(np.abs(yh*1000 - yp_ml[ti])) #index where model y = median parcel y at stime
    
    # Load dbz cross sections and read into preset arrays
    dbz_xz[i,:,:] = ds.variables['dbz'][:].data[0,jz,yi,jx] #x-z cross section at median parcel y
    dbz_yz[i,:,:] = ds.variables['dbz'][:].data[0,jz,jy,xi] #y-z cross section at median parcel x
    
    ds.close()
    
    
    ### Make plots
    
    # Plot X-Z dbz cross section
    if plot_xz:
        fig,ax = plt.subplots(1, 1, figsize=(10,4), layout='constrained')
        
        c = ax.contourf(xh[jx], zh[jz], dbz_xz, levels=np.arange(0,65,5), vmin=0, vmax=60,
                        cmap='NWSRef', antialiased=True)
        c.set_edgecolor('face')
        cb = plt.colorbar(c, ax=ax, extend='both')
        cb.set_label('dBZ', fontsize=12)
        ax.set_title(f"DBZ X-Z cross section at {stime/60} min", fontsize=14)
        ax.set_xlabel('x (km)', fontsize=14)
        ax.set_ylabel('z (km)', fontsize=14)
        ax.set_xlim(xl) #change these as needed
        ax.set_ylim(zl)
        
        plt.show()
    
    
    # Plot Y-Z dbz cross section
    if plot_yz:
        fig,ax = plt.subplots(1, 1, figsize=(10,4), layout='constrained')
        
        c = ax.contourf(yh[jy], zh[jz], dbz_yz, levels=np.arange(0,65,5), vmin=0, vmax=60,
                        cmap='NWSRef', antialiased=True)
        c.set_edgecolor('face')
        cb = plt.colorbar(c, ax=ax, extend='both')
        cb.set_label('dBZ', fontsize=12)
        ax.set_title(f"DBZ Y-Z cross section at {stime/60} min", fontsize=14)
        ax.set_xlabel('y (km)', fontsize=14)
        ax.set_ylabel('z (km)', fontsize=14)
        ax.set_xlim(yl) #change these as needed
        ax.set_ylim(zl)
        
        plt.show()
    
    
    
    
    
    