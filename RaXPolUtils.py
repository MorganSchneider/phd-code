# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 12:43:23 2023
@author: morgan.schneider

Functions for reading and plotting RaXPol data
"""

####################
### Load modules ###
####################

import matplotlib.pyplot as plt
import pandas as pd
import netCDF4 as nc
import numpy as np
import pyart
import matplotlib.dates as mdates
from glob import glob
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
from datetime import datetime
import cmocean
import wrf
# import cartopy.crs as ccrs

import warnings
warnings.filterwarnings('ignore', category=UserWarning)

#%%
########################
### Define colormaps ###
########################

# nwsdmap = LinearSegmentedColormap.from_list('nwsdmap',
#         np.vstack((plt.cm.gray(np.linspace(0,1,12)), pyart.graph.cm.RefDiff(np.linspace(0,1,12)))))
nwsdmap = LinearSegmentedColormap.from_list('nwsdmap',
        np.vstack((plt.cm.gray(np.linspace(0,1,12)), pyart.graph.cm_colorblind.HomeyerRainbow(np.linspace(0,1,12)))))

# Define colormaps for plotting radar variables
cmaps = {
    'dbz': {'cm': 'pyart_NWSRef', 'label': "$Z_H$ (dBZ)"},
    'vel': {'cm': 'pyart_Carbone42', 'label': "$v$ (m s$^{-1}$)"},
    'sw':  {'cm': 'pyart_NWS_SPW', 'label': "SW (m s$^{-1}$)"},
    'zdr': {'cm': nwsdmap, 'label': "$Z_{DR}$ (dB)"},
    'rhohv': {'cm': 'pyart_LangRainbow12', 'label': "\u03C1$_{HV}$"},
    'temp': {'cm': 'pyart_HomeyerRainbow', 'label': "$T$ (C)"},
    'dewpt': {'cm': cmocean.cm.algae, 'label': "$T_D$ (C)"},
    'theta': {'cm': 'pyart_HomeyerRainbow', 'label': "\u03B8 (K)"},
    'vort': {'cm': 'pyart_ChaseSpectral', 'label': "Pseudovorticity (s$^{-1}$)"},
    'div': {'cm': 'pyart_balance', 'label': "Pseudodivergence (s$^{-1}$)"}
}

#%%

# Calculate a moving-block average
def movmean(data, npts):
    # data: 1-D vector of data
    # npts: number of points to average
    data_mean = np.convolve(data, np.ones(npts), 'same')/npts
    return data_mean


# Convert lat/lon coordinates to x/y distances relative to an origin point
def latlon2xy(lat, lon, lat_o, lon_o):
    # lat, lon:     1-D vectors of lat/lon in decimal deg N/deg E
    # lat_o, lon_o: lat/lon of origin in decimal deg N/deg E
    
    r_earth = 6378.1
    
    thy = lat_o*np.pi/180
    thz = -lon_o*np.pi/180
    
    Ry = [[np.cos(thy), 0, np.sin(thy)],
          [0, 1, 0],
          [-np.sin(thy), 0, np.cos(thy)]]
    
    Rz = [[np.cos(thz), -np.sin(thz), 0],
          [np.sin(thz), np.cos(thz), 0],
          [0, 0, 1]]
    
    R = np.matmul(Ry,Rz)
    xyz = r_earth * np.array([np.cos(lat*np.pi/180) * np.cos(lon*np.pi/180),
                np.cos(lat*np.pi/180) * np.sin(lon*np.pi/180),
                np.sin(lat*np.pi/180)])
    
    posx = np.matmul(R[1],xyz)
    posy = np.matmul(R[2],xyz)
    
    return posx,posy


# Read RaXPol .nc files and load data into struct
def read_raxpol(fn):
    # fn: full path and filename
    
    ds = nc.Dataset(fn)
    
    vol_num = ds.variables['volume_number'][:].data # volume index
    swp_start = ds.variables['time_coverage_start'][:].data #UTC
    swp_end = ds.variables['time_coverage_end'][:].data #UTC
    swp_num = ds.variables['sweep_number'][:].data #relative to volume start
    fixed_angle = ds.variables['fixed_angle'][:].data #target elevation
    
    lat = ds.variables['latitude'][:].data #deg N
    lon = ds.variables['longitude'][:].data #deg E
    alt = ds.variables['altitude'][:].data #m
    
    time = ds.variables['time'][:].data #seconds since volume start
    r = ds.variables['range'][:].data #m
    az = ds.variables['azimuth'][:].data
    el = ds.variables['elevation'][:].data
    
    Z = ds.variables['DBZ'][:].data
    vr = ds.variables['VEL'][:].data
    sw = ds.variables['WIDTH'][:].data
    zdr = ds.variables['ZDR'][:].data
    rhohv = ds.variables['RHOHV'][:].data
    phidp = ds.variables['PHIDP'][:].data
    
    pulse_width = list(set(ds.variables['pulse_width'][:].data))[0]/1000 #s
    prt = list(set(ds.variables['prt'][:].data))[0] #s
    va = list(set(ds.variables['nyquist_velocity'][:].data))[0] #m/s
    
    ds.close()
    
    dr = r[1] - r[0]
    elev = np.mean(el).round(1)
    num_gates = len(r)
    num_az = len(az)
    iaz = az.argsort()
    
    az = az[iaz]
    Z = Z[iaz,:]
    vr = vr[iaz,:]
    sw = sw[iaz,:]
    zdr = zdr[iaz,:]
    rhohv = rhohv[iaz,:]
    phidp = phidp[iaz,:]
    
    r_km = np.linspace(dr, dr*num_gates, num_gates)/1000
    az_rad = az*np.pi/180
    
    az_mat = np.transpose(np.tile(az_rad, (num_gates, 1)))
    r_mat = np.tile(r_km, (num_az, 1))
    
    xx = r_mat * np.sin(az_mat) * np.cos(elev*np.pi/180)
    yy = r_mat * np.cos(az_mat) * np.cos(elev*np.pi/180)
    zz = r_mat * np.sin(elev*np.pi/180)
    
    dat = {'dbz':Z, 'vel':vr, 'sw':sw, 'zdr':zdr, 'rhohv':rhohv, 'phidp':phidp,
           'lat':lat, 'lon':lon, 'r':r_km, 'az':az, 'elev':elev, 'xx':xx, 'yy':yy, 'zz':zz,
           'vol_num':vol_num, 'swp_start':swp_start, 'swp_num':swp_num, 'num_gates':num_gates,
           'num_az':num_az, 'dr':dr, 'va':va, 'tau':pulse_width, 'prt':prt}
    
    return dat


# Read mobile mesonet .dat files and load data into struct
def read_MM(fn):
    # 3 - Derived_WS (m/s)
    # 4 - Derived_WD (deg)
    # 5 - SlowT (C)
    # 7 - UT_FastT (C)
    # 8 - Tdc (C)
    # 10 - Pressure (mb)
    # 11 - CompassDir (deg)
    # 12 - GPSDate (DDMMYY)
    # 13 - GPSTime (HHMMSS)
    # 15 - GPS_Lat (deg N)
    # 16 - GPS_Lon (deg E)
    # 17 - GPS_Alt (m)
    # 18 - GPS_Spd (m/s)
    # 19 - GPS_Dir (deg)
    # 24 - Flux_Dir (deg)
    
    cols = [3,4,5,7,8,10,11,12,13,15,16,17,18,19]
    df = pd.read_fwf(fn, sep='\t', usecols=cols, engine='python')
    
    #datestr = [str(dat.date[i]) if len(str(dat.date[i]))==6 else '0'+str(dat.date[i]) for i in range(len(df))]
    #timestr = [str(dat.time[i]) if len(str(dat.time[i]))==6 else '0'+str(dat.time[i]) for i in range(len(df))]
    
    compass_dir = df['CompassDir'].values
    gps_dir = df['GPS_Dir'].values
    gps_spd = df['GPS_Spd'].values
    
    veh_dir = np.array([gps_dir[i] if gps_spd[i]>=1 else compass_dir[i] for i in range(len(df))])
    veh_spd = np.array([gps_spd[i] if gps_spd[i]>=1 else 0 for i in range(len(df))])
    
    # Observations
    uwind = -df['Derived_WS'].values * np.sin(df['Derived_WD'].values*np.pi/180)
    vwind = -df['Derived_WS'].values * np.cos(df['Derived_WD'].values*np.pi/180)
    
    dat = {'date':df['GPSDate'].values, 'time':df['GPSTime'].values, 'veh_dir':veh_dir,
           'veh_spd':veh_spd, 'lat':df['GPS_Lat'].values, 'lon':df['GPS_Lon'].values,
           'alt':df['GPS_Alt'].values, 'SlowT':df['SlowT'].values, 'Utube':df['UT_FastT'].values,
           'Dewpt':df['Tdc'].values, 'Pres':df['Pressure'].values, 'wspd':df['Derived_WS'].values,
           'wdir':df['Derived_WD'].values, 'uwind':uwind, 'vwind':vwind}
    
    if "Probe_1" in fn: # P1 Utube temp is biased high
        dat.update({'Utube':dat['Utube']-1.23})
        
    # Time averaging
    av_period = [3, 10, 60]
    for n in av_period:
        dat.update({f"SlowT_{n}s":movmean(dat['SlowT'],n), f"Utube_{n}s":movmean(dat['Utube'],n),
                    f"Dewpt_{n}s":movmean(dat['Dewpt'],n), f"Pres_{n}s":movmean(dat['Pres'],n),
                    f"wspd_{n}s":movmean(dat['wspd'],n), f"uwind_{n}s":movmean(dat['uwind'],n),
                    f"vwind_{n}s":movmean(dat['vwind'],n)})
    
    if True:
        T_K = dat['Utube'] + 273.15
        Td_K = dat['Dewpt'] + 273.15
        e = 6.11*np.exp(2.5e6/461.5 * (1/273.15 - 1/Td_K))
        qv = 0.622 * e/(dat['Pres'] - e)
        
        theta = T_K * (1000/dat['Pres'])**0.286
        theta_v = theta * (1 + 0.61*qv)
        theta_e = (T_K + (2.5e6*qv)/(1004.5*T_K)) * (1000/dat['Pres'])**0.286
        
        dat.update({'Theta':theta, 'ThetaV':theta_v, 'ThetaE':theta_e})
        for n in av_period:
            dat.update({f"Theta_{n}s":movmean(theta,n), f"ThetaV_{n}s":movmean(theta_v,n),
                        f"ThetaE_{n}s":movmean(theta_e,n)})
    
    return dat


# Function to create 2D colorfill plots
def plot_cfill(x, y, data, field, ax, datalims=None, xlims=None, ylims=None,
               cmap=None, cbar=True, **kwargs):
    if cmap is None:
        cm, cb_label = cmaps[field]['cm'], cmaps[field]['label']
    else:
        cm, cb_label = cmap, cmaps[field]['label']
    
    if datalims is None:
        datamin = None
        datamax = None
    else:
        datamin = datalims[0]
        datamax = datalims[1]
    
    # Create the plot
    c = ax.pcolormesh(x, y, data, vmin=datamin, vmax=datamax, cmap=cm, **kwargs)

    # Format the colorbar
    # c.cmap.set_bad('grey', 1.0)
    if cbar is True:
        cb = plt.colorbar(c, ax=ax)
        cb.set_label(cb_label)
    
    if xlims is not None:
        ax.set_xlim(xlims[0], xlims[1])
    if ylims is not None:
        ax.set_ylim(ylims[0], ylims[1])
    
    return c


# Wrapper function for contourf
def plot_contourf(x, y, data, field, ax, levels=None, datalims=None, xlims=None, ylims=None,
                  cmap=None, cbar=True, cbfs=None, **kwargs):
    if cmap is None:
        cm, cb_label = cmaps[field]['cm'], cmaps[field]['label']
    else:
        cm, cb_label = cmap, cmaps[field]['label']
    
    if levels is None:
        levs = None
    else:
        levs = levels
    
    if datalims is None:
        datamin = None
        datamax = None
    else:
        datamin = datalims[0]
        datamax = datalims[1]
    
    c = ax.contourf(x, y, data, levels=levs, vmin=datamin, vmax=datamax, cmap=cm, antialiased=True, **kwargs)
    # ax.contour(x, y, data, levels=levs, vmin=datamin, vmax=datamax, cmap=cm, antialiased=True, **kwargs)
    # ax.contourf(x, y, data, levels=levs, vmin=datamin, vmax=datamax, cmap=cm, antialiased=True, **kwargs)
    # ax.contourf(x, y, data, levels=levs, vmin=datamin, vmax=datamax, cmap=cm, antialiased=True, **kwargs)
    c.set_edgecolor('face')
    
    if cbar:
        cb = plt.colorbar(c, ax=ax, extend='both')
        cb.set_label(cb_label)
        if np.max(np.abs(datalims)) < 0.1:
            cb.formatter.set_powerlimits((0,0))
        if cbfs is None:
            cb.set_label(cb_label)
        else:
            cb.set_label(cb_label, fontsize=cbfs)
    
    if xlims is not None:
        ax.set_xlim(xlims[0], xlims[1])
    if ylims is not None:
        ax.set_ylim(ylims[0], ylims[1])
    
    return c


































