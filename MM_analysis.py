#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  4 14:52:22 2025

@author: morgan.schneider
"""

from RaxpolUtils import *
# from MM_functions import *

#% Load MM data and vortex locations

# Probe 1 - Tyler's bias-corrected files
ds = nc.Dataset('/Users/morgan.schneider/Documents/perils2023/iop2/mesonet/Probe_1_IOP2_QC_all.nc')
P1 = dict(time=ds.variables['time'][:].data+21600,
          lat=ds.variables['lat'][:].data,
          lon=ds.variables['lon'][:].data,
          elev=ds.variables['elev'][:].data,
          temp_raw=ds.variables['temp'][:].data,
          pres_raw=ds.variables['pres'][:].data,
          rh_raw=ds.variables['rh'][:].data,
          wspd_raw=ds.variables['wspd'][:].data,
          wdir_raw=ds.variables['wdir'][:].data,
          wspd_corr=ds.variables['wspd_corr'][:].data,
          wdir_corr=ds.variables['wdir_corr'][:].data,
          u_corr=ds.variables['u_corr'][:].data,
          v_corr=ds.variables['v_corr'][:].data,
          temp_corr=ds.variables['temp_unbiased'][:].data,
          rh_corr=ds.variables['rh_unbiased'][:].data,
          pres_corr=ds.variables['pres_unbiased'][:].data)
ds.close()

# Probe 2 - Tyler's bias-corrected files
ds = nc.Dataset('/Users/morgan.schneider/Documents/perils2023/iop2/mesonet/Probe_2_IOP2_QC_all.nc')
P2 = dict(time=ds.variables['time'][:].data+21600,
          lat=ds.variables['lat'][:].data,
          lon=ds.variables['lon'][:].data,
          elev=ds.variables['elev'][:].data,
          temp_raw=ds.variables['temp'][:].data,
          pres_raw=ds.variables['pres'][:].data,
          rh_raw=ds.variables['rh'][:].data,
          wspd_raw=ds.variables['wspd'][:].data,
          wdir_raw=ds.variables['wdir'][:].data,
          wspd_corr=ds.variables['wspd_corr'][:].data,
          wdir_corr=ds.variables['wdir_corr'][:].data,
          u_corr=ds.variables['u_corr'][:].data,
          v_corr=ds.variables['v_corr'][:].data,
          temp_corr=ds.variables['temp_unbiased'][:].data,
          rh_corr=ds.variables['rh_unbiased'][:].data,
          pres_corr=ds.variables['pres_unbiased'][:].data)
ds.close()


# P1['time'] = P1['time'] + 21600
P1_times = np.array( [ datetime.fromtimestamp(P1['time'][i]) for i in range(len(P1['time'])) ] )
P2_times = np.array( [ datetime.fromtimestamp(P2['time'][i]) for i in range(len(P2['time'])) ] )


fn = glob(f"/Users/morgan.schneider/Documents/perils2023/iop2/raxpol/CFradial_cal/*_081928_*.nc")[0]
rax = pyart.io.read(fn)
rax_lat = rax.latitude['data'][0]
rax_lon = rax.longitude['data'][0]


meso_x = np.array([-2.2, -1.8, -1.5, -1.2, -0.4, -0.05, 0.35, 0.65, 1.0, 1.3])
meso_y = np.array([1.5, 2.0, 2.3, 2.7, 3.0, 3.5, 3.9, 4.3, 4.7, 5.2])
meso_x = np.linspace(-2.2, 1.3, 11)
meso_y = np.linspace(1.5, 5.2, 11)
meso_lats,meso_lons = xy2latlon(meso_x, meso_y, rax_lat, rax_lon)

# convert all the vortex loc x/y to lat/lons, then figure out how to find and track a center point
# --> center on the visible couplet (vortex 3)? use the scans where it's visible to get couplet motion
# and then back out the position for the earlier times

meso_times = [datetime(2023,3,3,8,20,0), datetime(2023,3,3,8,20,30), datetime(2023,3,3,8,21,0),
              datetime(2023,3,3,8,21,30), datetime(2023,3,3,8,22,0), datetime(2023,3,3,8,22,30), datetime(2023,3,3,8,23,0),
              datetime(2023,3,3,8,23,30), datetime(2023,3,3,8,24,0), datetime(2023,3,3,8,24,30),
              datetime(2023,3,3,8,25,0)]


P1_distances = compute_distances(meso_times, meso_lats, meso_lons, P1_times, P1['lat'], P1['lon'])
P2_distances = compute_distances(meso_times, meso_lats, meso_lons, P2_times, P2['lat'], P2['lon'])

P1_dx = P1_distances[:,0]
P1_dy = P1_distances[:,1]
P2_dx = P2_distances[:,0]
P2_dy = P2_distances[:,1]

P1_x = P1_dx[(~np.isnan(P1_dx))]
P1_y = P1_dy[(~np.isnan(P1_dy))]
P1_temp = P1['temp_corr'][(~np.isnan(P1_dx))]
P1_u = P1['u_corr'][(~np.isnan(P1_dx))]
P1_v = P1['v_corr'][(~np.isnan(P1_dx))]
P1_wspd = P1['wspd_corr'][(~np.isnan(P1_dx))]
P1_wdir = P1['wdir_corr'][(~np.isnan(P1_dx))]
P2_x = P2_dx[(~np.isnan(P2_dx))]
P2_y = P2_dy[(~np.isnan(P2_dy))]
P2_temp = P2['temp_corr'][(~np.isnan(P2_dx))]
P2_u = P2['u_corr'][(~np.isnan(P2_dx))]
P2_v = P2['v_corr'][(~np.isnan(P2_dx))]
P2_wspd = P2['wspd_corr'][(~np.isnan(P2_dx))]
P2_wdir = P2['wdir_corr'][(~np.isnan(P2_dx))]


# saved raxpol data from circuit
if 'vol' not in locals():
    dbfile = open('/Users/morgan.schneider/Documents/perils2023/iop2/circuit_data.pkl', 'rb')
    dat = pickle.load(dbfile)
    vol = dat['Rax']
    dbfile.close()

dbfile = open('/Users/morgan.schneider/Documents/perils2023/iop2/raxpol_vortex_locs.pkl', 'rb')
locs = pickle.load(dbfile)
dbfile.close()


vi = 13
eli = 1
filetime = vol[vi]['scan_time'][eli]
# vortex_num = 3
# vortex 1: vi = 7-9
# vortex 2: vi = 8-13
# vortex 3: vi = 8-17
# vortex 4: vi = 9-12
# rotor: vi = 8

x_rot = np.append(locs[filetime]['vortex2']['x'], locs[filetime]['vortex3']['x'])
y_rot = np.append(locs[filetime]['vortex2']['y'], locs[filetime]['vortex3']['y'])


#%%

ip = '/Users/morgan.schneider/Documents/perils2023/iop2/figs/'

figsave = False


if True:
    xl = [-4, 4]
    yl = [0, 8]
    
    datalims = [19,21]
    
    fig,(ax1,ax2) = plt.subplots(1,2,figsize=(11.5,4), sharex=True, sharey=True, subplot_kw=dict(box_aspect=1), layout='constrained')
    
    plot_cfill(vol[vi]['xx'][eli,:,:], vol[vi]['yy'][eli,:,:], vol[vi]['dbz'][eli,:,:], 'dbz', ax1, datalims=[0,70], xlims=xl, ylims=yl)
    ax1.scatter(x_rot, y_rot, s=30, c='k', marker='.')
    ax1.set_title(f"{filetime} UTC {vol[vi]['elev'][eli].round(1)}\N{DEGREE SIGN} reflectivity", fontsize=14)
    ax1.set_xlabel('E-W distance from radar (km)', fontsize=12)
    ax1.set_ylabel('N-S distance from radar (km)', fontsize=12)
    b1 = ax1.barbs(P1_x[::30], P1_y[::30]+3.5, P1_u[::30], P1_v[::30], barbcolor='k', length=7)
    b2 = ax1.barbs(P2_x[::30], P2_y[::30]+3.5, P2_u[::30], P2_v[::30], barbcolor='k', length=7)
    s1 = ax1.scatter(P1_x[::30], P1_y[::30]+3.5, s=50, c=P1_temp[::30], cmap=cmaps['temp']['cm'],
                    vmin=datalims[0], vmax=datalims[1], marker='s', edgecolors='k')
    s2 = ax1.scatter(P2_x[::30], P2_y[::30]+3.5, s=50, c=P2_temp[::30], cmap=cmaps['temp']['cm'],
                    vmin=datalims[0], vmax=datalims[1], marker='^', edgecolors='k')
    ax1.scatter(0, 0, s=50, c='k')
    ax1.text(-1, 0.4, 'RaXPol', fontsize=12, fontweight='bold')
    plt.colorbar(s1,label='MM temperature (C)')
    ax1.legend([s1,s2], ['Probe 1', 'Probe 2'], loc='upper left')
    # ax1.plot(vol[vi]['xx'][eli,azi,:], vol[vi]['yy'][eli,azi,:], '--k', linewidth=1.25)
    
    plot_cfill(vol[vi]['xx'][eli,:,:], vol[vi]['yy'][eli,:,:], vol[vi]['vel'][eli,:,:], 'vel', ax2, datalims=[-30,30], xlims=xl, ylims=yl)
    ax2.scatter(x_rot, y_rot, s=30, c='k', marker='.')
    ax2.set_title(f"{filetime} UTC {vol[vi]['elev'][eli].round(1)}\N{DEGREE SIGN} radial velocity", fontsize=14)
    ax2.set_xlabel('E-W distance from radar (km)', fontsize=12)
    ax2.set_ylabel('N-S distance from radar (km)', fontsize=12)
    b1 = ax2.barbs(P1_x[::30], P1_y[::30]+3.5, P1_u[::30], P1_v[::30], barbcolor='k', length=7)
    b2 = ax2.barbs(P2_x[::30], P2_y[::30]+3.5, P2_u[::30], P2_v[::30], barbcolor='k', length=7)
    s1 = ax2.scatter(P1_x[::30], P1_y[::30]+3.5, s=50, c=P1_temp[::30], cmap=cmaps['temp']['cm'],
                    vmin=datalims[0], vmax=datalims[1], marker='s', edgecolors='k')
    s2 = ax2.scatter(P2_x[::30], P2_y[::30]+3.5, s=50, c=P2_temp[::30], cmap=cmaps['temp']['cm'],
                    vmin=datalims[0], vmax=datalims[1], marker='^', edgecolors='k')
    ax2.scatter(0, 0, s=50, c='k')
    ax2.text(-1, 0.4, 'RaXPol', fontsize=12, fontweight='bold')
    plt.colorbar(s1,label='MM temperature (C)')
    # plt.legend([s1,s2], ['Probe 1', 'Probe 2'], loc='upper left')
    # ax2.plot(vol[vi]['xx'][eli,azi,:], vol[vi]['yy'][eli,azi,:], '--k', linewidth=1.25)
    
    # plt.suptitle(f"{filetime} UTC, azimuth = {azimuth}\N{DEGREE SIGN}", fontsize=14)
    # plt.suptitle(f"{filetime} UTC", fontsize=14)
    if figsave:
        plt.savefig(ip+f"vol{vi}_{filetime}_PPI_MMCircuit.png", dpi=300)






