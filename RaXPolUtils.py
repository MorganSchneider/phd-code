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
from datetime import datetime,timedelta
import cmocean
import wrf # current version needs python 3.10 or 3.11 and need to install scipy first if i'm using that too
import pickle
# import cartopy.crs as ccrs
from os.path import exists
# import scipy

import warnings
warnings.filterwarnings('ignore', category=UserWarning)

#%%
########################
### Define colormaps ###
########################

# I made a custom ZDR map because pyart didn't have the one from radarscope
nwsdmap = LinearSegmentedColormap.from_list('nwsdmap',
        np.vstack((plt.cm.gray(np.linspace(0,1,12)), pyart.graph.cmweather.cm_colorblind.HomeyerRainbow(np.linspace(0,1,12)))))

# Define colormaps for plotting radar variables
cmaps = {
    'dbz': {'cm': 'NWSRef', 'label': "$Z_H$ (dBZ)"},
    'vel': {'cm': 'balance', 'label': "$v$ (m s$^{-1}$)"},
    'sw':  {'cm': 'NWS_SPW', 'label': "SW (m s$^{-1}$)"},
    'zdr': {'cm': nwsdmap, 'label': "$Z_{DR}$ (dB)"},
    'rhohv': {'cm': 'LangRainbow12', 'label': "\u03C1$_{HV}$"},
    'temp': {'cm': 'HomeyerRainbow', 'label': "$T$ (C)"},
    'dewpt': {'cm': cmocean.cm.algae, 'label': "$T_D$ (C)"},
    'theta': {'cm': 'HomeyerRainbow', 'label': "\u03B8 (K)"},
    'vort': {'cm': 'ChaseSpectral', 'label': "Pseudovorticity (s$^{-1}$)"},
    'div': {'cm': 'balance', 'label': "Pseudodivergence (s$^{-1}$)"}
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
    # lat, lon:     1-D vectors of lat/lon in decimal degrees N/deg E
    # lat_o, lon_o: lat/lon of origin in decimal degrees N/deg E
    
    r_earth = 6378.1 # km
    
    thy = lat_o*np.pi/180 # convert to radians
    thz = -lon_o*np.pi/180
    
    # transform matrices
    Ry = [[np.cos(thy), 0, np.sin(thy)],
          [0, 1, 0],
          [-np.sin(thy), 0, np.cos(thy)]]
    
    Rz = [[np.cos(thz), -np.sin(thz), 0],
          [np.sin(thz), np.cos(thz), 0],
          [0, 0, 1]]
    
    # i'm actually not sure exactly how this works, i just copied this function from some of Boonleng's code
    R = np.matmul(Ry,Rz)
    xyz = r_earth * np.array([np.cos(lat*np.pi/180) * np.cos(lon*np.pi/180),
                np.cos(lat*np.pi/180) * np.sin(lon*np.pi/180),
                np.sin(lat*np.pi/180)])
    
    # get x and y positions
    posx = np.matmul(R[1],xyz)
    posy = np.matmul(R[2],xyz)
    
    return posx,posy


# Convert x/y distances relative to an origin point back to lat/lon
def xy2latlon(x, y, lat_o, lon_o):
    # x, y: 1-D vectors of xy positions relative to an origin point in km
    # lat_o, lon_o: lat/lon of origin in decimal deg N/deg E
    
    r_earth = 6378.1 # km
    
    x_o = r_earth * np.cos(lat_o*np.pi/180) * np.sin(lon_o*np.pi/180)
    y_o = r_earth * np.sin(lat_o*np.pi/180)
    
    lat = np.arcsin(y / r_earth)
    lon = np.arcsin(x / (r_earth * np.cos(lat)))
    
    lat_deg = (lat * 180/np.pi) + lat_o
    lon_deg = (lon * 180/np.pi) + lon_o
    
    return lat_deg,lon_deg


# Compute distances of mobile mesonet data to meso (adapted from Tyler Pardun's code)
def compute_distances(meso_times, meso_lats, meso_lons, point_times, point_lats, point_lons, max_time_diff=60):
    """
    Compute distances from each point to the mesocyclone location at the closest time step.

    Parameters:
    - meso_times (array-like): Timestamps of mesocyclone locations.
    - meso_lats (array-like): Latitudes of mesocyclone locations.
    - meso_lons (array-like): Longitudes of mesocyclone locations.
    - point_times (array-like): Timestamps of points.
    - point_lats (array-like): Latitudes of points.
    - point_lons (array-like): Longitudes of points.
    - max_time_diff (int): Maximum allowed time difference in seconds.

    Returns:
    - distances (array): Distances (km) between mesocyclone and points.
    - matched_times (array): Matched mesocyclone times for each point.
    """

    distances = []
    matched_times = []

    # Convert times to pandas datetime for easy matching
    meso_times = pd.to_datetime(meso_times)
    point_times = pd.to_datetime(point_times)

    for p_time, p_lat, p_lon in zip(point_times, point_lats, point_lons):
        # Find the mesocyclone time closest to the point time
        time_diffs = np.abs((meso_times - p_time).total_seconds())
        closest_idx = np.argmin(time_diffs)
        
        # Check if the time difference is within the allowed window
        if time_diffs[closest_idx] <= max_time_diff:
            matched_times.append(meso_times[closest_idx])
            
            # Compute distance
            meso_loc = (meso_lats[closest_idx], meso_lons[closest_idx])
            point_loc = (p_lat, p_lon)
            # distance_km = get_dx_dy(meso_lons[closest_idx],meso_lats[closest_idx],p_lon,p_lat)
            dx_km,dy_km = latlon2xy(p_lat,p_lon,meso_lats[closest_idx],meso_lons[closest_idx])
            distance_km = np.array([dx_km, dy_km])
            distances.append(distance_km)
        else:
            # If no close time match, store NaN
            matched_times.append(np.nan)
            distances.append([np.nan,np.nan])

    return np.array(distances)



# Grid radar data to a Cartesian grid using fast Barnes interpolation
def grid_PPI(file):
    from fastbarnes.interpolation import barnes
    
    

# Read RaXPol CFradial files and load data into dict because it's easier than digging for the data in pyart objects
def read_raxpol(file):
    # file: full file path
    ds = nc.Dataset(file)
    
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


# Read mobile mesonet .dat files and load data into dict
def read_MM_dat(file):
    # file = full file path
    
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
    df = pd.read_fwf(file, sep='\t', usecols=cols, engine='python')
    
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
           'alt':df['GPS_Alt'].values, 'Tslow':df['SlowT'].values, 'Tfast':df['UT_FastT'].values,
           'Dewpt':df['Tdc'].values, 'Pres':df['Pressure'].values, 'wspd':df['Derived_WS'].values,
           'wdir':df['Derived_WD'].values, 'uwind':uwind, 'vwind':vwind}
    
    if "Probe_1" in file: # P1 temp is biased high
        dat.update({'Tfast':dat['Tfast']-1.23})
        
    # Time averaging
    av_period = [3, 10, 60]
    for n in av_period:
        dat.update({f"Tslow_{n}s":movmean(dat['Tslow'],n), f"Tfast_{n}s":movmean(dat['Tfast'],n),
                    f"Dewpt_{n}s":movmean(dat['Dewpt'],n), f"Pres_{n}s":movmean(dat['Pres'],n),
                    f"wspd_{n}s":movmean(dat['wspd'],n), f"uwind_{n}s":movmean(dat['uwind'],n),
                    f"vwind_{n}s":movmean(dat['vwind'],n)})
    
    if True:
        T_K = dat['Tfast'] + 273.15
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


# Wrapper function for pcolormesh because i like my code to be readable
def plot_cfill(x, y, data, field, ax, datalims=None, xlims=None, ylims=None,
               cmap=None, cbar=True, cbfs=None, cbts=None, **kwargs):
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
        if np.max(np.abs(datalims)) < 0.1:
            cb.formatter.set_powerlimits((0,0))
        if cbfs is None:
            cb.set_label(cb_label)
        else:
            cb.set_label(cb_label, fontsize=cbfs)
        if cbts is not None:
            cb.ax.tick_params(labelsize=cbts)
    
    if xlims is not None:
        ax.set_xlim(xlims[0], xlims[1])
    if ylims is not None:
        ax.set_ylim(ylims[0], ylims[1])
    
    return c


# Wrapper function for contourf because i like my code to be readable
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


# Automate saving data to new or existing pickle file because i'm insane
def save_to_pickle(data, pkl_fname, new_pkl=False):
    # data: dict of variables to save
    # pkl_fname: filename to save data to (includes path)
    # new_pkl: True or False (if True, will overwrite any existing file)
    if (not exists(pkl_fname)) | new_pkl:
        dbfile = open(pkl_fname, 'wb')
        pickle.dump(data, dbfile)
        dbfile.close()
    elif exists(pkl_fname):
        dbfile = open(pkl_fname, 'rb')
        save_data = pickle.load(dbfile)
        dbfile.close()
        
        save_data.update(data)
        dbfile = open(pkl_fname, 'wb')
        pickle.dump(save_data, dbfile)
        dbfile.close()































