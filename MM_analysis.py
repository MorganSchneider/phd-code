#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  4 14:52:22 2025

@author: morgan.schneider
"""

from RaxpolUtils import *
# from MM_functions import *

#% Load MM data and vortex locations

# Temperature bias: +/- 1.12 C ---> difference betweem P1 and P2 unbiased temperatures is 0.87 C
# Pressure bias: +/- 0.2 hPa ---> difference between P1 and P2 unbiased pres is 0.56 hPa, but doesn't look systematic
# RH bias: +/- 0.05 ---> difference between P1 and P2 unbiased RH is 0.031
# u bias: +/- 0.16 m/s
# v bias: +/- 0.02 m/s


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
          temp_unbiased=ds.variables['temp_unbiased'][:].data,
          rh_unbiased=ds.variables['rh_unbiased'][:].data,
          pres_unbiased=ds.variables['pres_unbiased'][:].data)
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
          temp_unbiased=ds.variables['temp_unbiased'][:].data-0.88,
          rh_unbiased=ds.variables['rh_unbiased'][:].data+0.03,
          pres_unbiased=ds.variables['pres_unbiased'][:].data-0.56)
ds.close()







# P1['time'] = P1['time'] + 21600
P1_times = np.array( [ datetime.fromtimestamp(P1['time'][i]) for i in range(len(P1['time'])) ] )
P2_times = np.array( [ datetime.fromtimestamp(P2['time'][i]) for i in range(len(P2['time'])) ] )


fn = glob(f"/Users/morgan.schneider/Documents/perils2023/iop2/raxpol/CFradial_cal/*_081928_*.nc")[0]
rax = pyart.io.read(fn)
rax_lat = rax.latitude['data'][0]
rax_lon = rax.longitude['data'][0]


# based on reflectivity leading edge
# meso_x = np.array([-3.2, -2.7, -2.2, -1.8, -1.4, -1.0, -0.5, -0.1, 0.3, 0.7, 1.1, 1.4, 1.9])
# meso_y = np.array([ 0.2,  0.6,  1.0,  1.4,  1.8,  2.2,  2.6,  3.0, 3.3, 3.6, 3.9, 4.2, 4.6])

meso_x = np.array([-2.6, -1.95, -1.6, -1.3, -1.05, -0.15, 0.3, 0.6, 0.9, 1.3]) #couplet positions at volume center times
meso_y = np.array([0.6, 1.2, 1.6, 2.0, 2.5, 3.6, 4.1, 4.7, 5.2, 5.8])

# based on velocity couplet/vortex locs
# meso_x = np.array([-2.6, -2.35, -2.1, -1.8, -1.5, -1.3, -1.0, -0.4, -0.15, 0.25, 0.6,  0.9, 1.2])
# meso_y = np.array([ 0.5,  0.9,   1.3,  1.7,  2.0,  2.25, 2.6,  2.9,  3.4,  3.9,  4.25, 4.8, 5.4])

meso_lats,meso_lons = xy2latlon(meso_x, meso_y, rax_lat, rax_lon)

# convert all the vortex loc x/y to lat/lons, then figure out how to find and track a center point
# --> center on the visible couplet (vortex 3)? use the scans where it's visible to get couplet motion
# and then back out the position for the earlier times

meso_times = [datetime(2023,3,3,8,19,10), datetime(2023,3,3,8,20,10),
              datetime(2023,3,3,8,20,40), datetime(2023,3,3,8,21,10),
              datetime(2023,3,3,8,21,36), datetime(2023,3,3,8,22,40),
              datetime(2023,3,3,8,23,10), datetime(2023,3,3,8,23,40),
              datetime(2023,3,3,8,24,10), datetime(2023,3,3,8,24,40)]


P1_distances = compute_distances(meso_times, meso_lats, meso_lons, P1_times, P1['lat'], P1['lon'])
P2_distances = compute_distances(meso_times, meso_lats, meso_lons, P2_times, P2['lat'], P2['lon'])

P1_dx = P1_distances[:,0]
P1_dy = P1_distances[:,1]
P2_dx = P2_distances[:,0]
P2_dy = P2_distances[:,1]

P1_x = P1_dx[(~np.isnan(P1_dx))]
P1_y = P1_dy[(~np.isnan(P1_dy))]
P1_temp = P1['temp_raw'][(~np.isnan(P1_dx))]-1.25
P1_pres = P1['pres_unbiased'][(~np.isnan(P1_dx))]
P1_rh = P1['rh_unbiased'][(~np.isnan(P1_dx))]
P1_u = P1['u_corr'][(~np.isnan(P1_dx))]
P1_v = P1['v_corr'][(~np.isnan(P1_dx))]
P1_wspd = P1['wspd_corr'][(~np.isnan(P1_dx))]
P1_wdir = P1['wdir_corr'][(~np.isnan(P1_dx))]
P2_x = P2_dx[(~np.isnan(P2_dx))]
P2_y = P2_dy[(~np.isnan(P2_dy))]
P2_temp = P2['temp_raw'][(~np.isnan(P2_dx))]
P2_pres = P2['pres_unbiased'][(~np.isnan(P2_dx))]
P2_rh = P2['rh_unbiased'][(~np.isnan(P2_dx))]
P2_u = P2['u_corr'][(~np.isnan(P2_dx))]
P2_v = P2['v_corr'][(~np.isnan(P2_dx))]
P2_wspd = P2['wspd_corr'][(~np.isnan(P2_dx))]
P2_wdir = P2['wdir_corr'][(~np.isnan(P2_dx))]


P1_theta = (P1_temp+273.15) * (1000/P1_pres)**0.286
P2_theta = (P2_temp+273.15) * (1000/P2_pres)**0.286

P1_es = 6.11*np.exp(2.5e6/461.5 * (1/273.15 - 1/(P1_temp+273.15)))
P2_es = 6.11*np.exp(2.5e6/461.5 * (1/273.15 - 1/(P2_temp+273.15)))

P1_e = P1_es * P1_rh
P2_e = P2_es * P2_rh

P1_qv = 0.622 * P1_e/(P1_pres - P1_e)
P2_qv = 0.622 * P2_e/(P2_pres - P2_e)

P1_thetav = P1_theta * (1 + 0.61*P1_qv)
P2_thetav = P2_theta * (1 + 0.61*P2_qv)

P1_thpert = P1_theta - 294.3
P2_thpert= P2_theta - 294.3
P1_thvpert = P1_thetav - 296.7
P2_thvpert = P2_thetav - 296.7


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

# x_rot = np.append(locs[filetime]['vortex2']['x'], locs[filetime]['vortex3']['x'])
# y_rot = np.append(locs[filetime]['vortex2']['y'], locs[filetime]['vortex3']['y'])

x_rot = np.array([])
y_rot = np.array([])

for key in list(locs[filetime].keys()):
    if 'vortex' in key:
        x_rot = np.append(x_rot, locs[filetime][key]['x'])
        y_rot = np.append(y_rot, locs[filetime][key]['y'])


#%%

ip = '/Users/morgan.schneider/Documents/perils2023/iop2/figs/'

figsave = False



vi = 16
eli = 1
filetime = vol[vi]['scan_time'][eli]

x_rot = np.array([])
y_rot = np.array([])
for key in list(locs[filetime].keys()):
    if 'vortex' in key:
        x_rot = np.append(x_rot, locs[filetime][key]['x'])
        y_rot = np.append(y_rot, locs[filetime][key]['y'])

x_rot = locs[filetime]['vortex3']['x']
y_rot = locs[filetime]['vortex3']['y']

sep = 30

if vi == 9:
    idx = slice(1,7)
elif vi == 13:
    idx = slice(6,12)
elif vi == 16:
    idx = slice(11,17)


x1 = P1_x[::sep][idx]
y1 = P1_y[::sep][idx]+3.3
u1 = P1_u[::sep][idx]
v1 = P1_v[::sep][idx]
T1 = P1_temp[::sep][idx]

x2 = P2_x[::sep][idx]
y2 = P2_y[::sep][idx]+3.3
u2 = P2_u[::sep][idx]
v2 = P2_v[::sep][idx]
T2 = P2_temp[::sep][idx]


if False:
    xl = [-5, 5]
    yl = [0, 10]
    
    datalims = [18,21]
    
    
    fig,(ax1,ax2) = plt.subplots(1,2,figsize=(11.5,4), sharex=True, sharey=True, subplot_kw=dict(box_aspect=1), layout='constrained')
    
    plot_cfill(vol[vi]['xx'][eli,:,:], vol[vi]['yy'][eli,:,:], vol[vi]['dbz'][eli,:,:], 'dbz', ax1, datalims=[0,70], xlims=xl, ylims=yl)
    # ax1.scatter(x_rot, y_rot, s=30, c='k', marker='.')
    ax1.set_title(f"{filetime} UTC {vol[vi]['elev'][eli].round(1)}\N{DEGREE SIGN} reflectivity", fontsize=14)
    ax1.set_xlabel('E-W distance from radar (km)', fontsize=12)
    ax1.set_ylabel('N-S distance from radar (km)', fontsize=12)
    b1 = ax1.barbs(x1, y1, u1, v1, barbcolor='k', length=7)
    b2 = ax1.barbs(x2, y2, u2, v2, barbcolor='k', length=7)
    s1 = ax1.scatter(x1, y1, s=50, c=T1, cmap=cmaps['temp']['cm'],
                    vmin=datalims[0], vmax=datalims[1], marker='s', edgecolors='k')
    s2 = ax1.scatter(x2, y2, s=50, c=T2, cmap=cmaps['temp']['cm'],
                    vmin=datalims[0], vmax=datalims[1], marker='^', edgecolors='k')
    ax1.scatter(0, 0, s=50, c='k')
    ax1.text(-1, 0.4, 'RaXPol', fontsize=12, fontweight='bold')
    plt.colorbar(s1,label='MM temperature (C)')
    ax1.legend([s1,s2], ['Probe 1', 'Probe 2'], loc='upper left')
    # ax1.plot(vol[vi]['xx'][eli,azi,:], vol[vi]['yy'][eli,azi,:], '--k', linewidth=1.25)
    
    plot_cfill(vol[vi]['xx'][eli,:,:], vol[vi]['yy'][eli,:,:], vol[vi]['vel'][eli,:,:], 'vel', ax2, datalims=[-30,30], xlims=xl, ylims=yl)
    # ax2.scatter(x_rot, y_rot, s=30, c='k', marker='.')
    ax2.set_title(f"{filetime} UTC {vol[vi]['elev'][eli].round(1)}\N{DEGREE SIGN} radial velocity", fontsize=14)
    ax2.set_xlabel('E-W distance from radar (km)', fontsize=12)
    ax2.set_ylabel('N-S distance from radar (km)', fontsize=12)
    b1 = ax2.barbs(x1, y1, u1, v1, barbcolor='k', length=7)
    b2 = ax2.barbs(x2, y2, u2, v2, barbcolor='k', length=7)
    s1 = ax2.scatter(x1, y1, s=50, c=T1, cmap=cmaps['temp']['cm'],
                    vmin=datalims[0], vmax=datalims[1], marker='s', edgecolors='k')
    s2 = ax2.scatter(x2, y2, s=50, c=T2, cmap=cmaps['temp']['cm'],
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



if True:
    xl = [-5, 5]
    yl = [0, 10]
    
    P1_c = P1_temp; P2_c = P2_temp; datalims = [19,20.4]; cmap = cmaps['temp']['cm']
    P1_c = P1_thvpert; P2_c = P2_thvpert; datalims = [-1, 1]; cmap = 'PuOr_r' #was 'bwr'
    
    sep = 30
    i1 = slice(1,7)
    i2 = slice(6,12)
    i3 = slice(11,17)
    # (1,7), (6,12), (11,17) for sep=30
    # (1,9), (9,17), (17,25) for sep=20
    
    u_P1 = P1_u[::sep]
    v_P1 = P1_v[::sep]
    for i in range(len(P1_u[::sep])):
        ind = i * sep
        if np.isnan(P1_u[::sep][i]):
            isnan = True
            ct = 0
            while (isnan) & (ct<11):
                if np.isnan(P1_u[ind-ct]) & np.isnan(P1_u[ind+ct]):
                    isnan = True
                    ct += 1
                elif ~np.isnan(P1_u[ind-ct]):
                    u_P1[i] = P1_u[ind-ct]
                    v_P1[i] = P1_v[ind-ct]
                    isnan = False
                    ct += 1
                elif ~np.isnan(P1_u[ind+ct]):
                    u_P1[i] = P1_u[ind+ct]
                    v_P1[i] = P1_v[ind+ct]
                    isnan = False
                    ct += 1
                    
    u_P2 = P2_u[::sep]
    v_P2 = P2_v[::sep]
    for i in range(len(P2_u[::sep])):
        ind = i * sep
        if np.isnan(P2_u[::sep][i]):
            isnan = True
            ct = 0
            while (isnan) & (ct<11):
                if np.isnan(P2_u[ind-ct]) & np.isnan(P2_u[ind+ct]):
                    isnan = True
                    ct += 1
                elif ~np.isnan(P2_u[ind-ct]):
                    u_P2[i] = P2_u[ind-ct]
                    v_P2[i] = P2_v[ind-ct]
                    isnan = False
                    ct += 1
                elif ~np.isnan(P2_u[ind+ct]):
                    u_P2[i] = P2_u[ind+ct]
                    v_P2[i] = P2_v[ind+ct]
                    isnan = False
                    ct += 1
    
    
    fig,ax = plt.subplots(2, 2, figsize=(8.25,8), sharex=True, sharey=True, subplot_kw=dict(box_aspect=1), layout='constrained')
    
    c = plot_cfill(vol[9]['xx'][eli,:,:], vol[9]['yy'][eli,:,:], np.ma.masked_array(vol[9]['vel'][eli,:,:], vol[9]['dbz'][eli,:,:]<1), 'vel', ax[0,0], datalims=[-30,30], xlims=xl, ylims=yl, cmap='balance', cbar=False)
    # ax[0,0].scatter(x_rot, y_rot, s=30, c='k', marker='.')
    # ax[0,0].set_title(f"{vol[9]['scan_time'][eli]} UTC {vol[9]['elev'][eli].round(1)}\N{DEGREE SIGN} radial velocity", fontsize=12)
    ax[0,0].set_ylabel('N-S distance from radar (km)', fontsize=14)
    # b1 = ax[0,0].barbs(P1_x[::sep][i1], P1_y[::sep][i1]+3.3, P1_u[::sep][i1], P1_v[::sep][i1], barbcolor='k', length=7)
    # b2 = ax[0,0].barbs(P2_x[::sep][i1], P2_y[::sep][i1]+3.3, P2_u[::sep][i1], P2_v[::sep][i1], barbcolor='k', length=7)
    # ax[0,0].quiver(P1_x[::sep][i1], P1_y[::sep][i1]+3.3, P1_u[::sep][i1], P1_v[::sep][i1], color='k', scale=75, width=0.008, pivot='tail')
    # ax[0,0].quiver(P2_x[::sep][i1], P2_y[::sep][i1]+3.3, P2_u[::sep][i1], P2_v[::sep][i1], color='k', scale=75, width=0.008, pivot='tail')
    s1 = ax[0,0].scatter(P1_x[::sep][i1], P1_y[::sep][i1]+3.3, s=60, c=P1_c[::sep][i1], cmap=cmap,
                    vmin=datalims[0], vmax=datalims[1], marker='s', edgecolors='k', linewidth=0.75)
    s2 = ax[0,0].scatter(P2_x[::sep][i1], P2_y[::sep][i1]+3.3, s=60, c=P2_c[::sep][i1], cmap=cmap,
                    vmin=datalims[0], vmax=datalims[1], marker='^', edgecolors='k', linewidth=0.75)
    # ax[0,0].quiver(P1_x[::sep][i1], P1_y[::sep][i1]+3.3, P1_u[::sep][i1], P1_v[::sep][i1], color='k', scale=75, width=0.008, pivot='tail')
    # ax[0,0].quiver(P2_x[::sep][i1], P2_y[::sep][i1]+3.3, P2_u[::sep][i1], P2_v[::sep][i1], color='k', scale=75, width=0.008, pivot='tail')
    ax[0,0].scatter(0, 0.15, s=50, c='k')
    ax[0,0].text(-4.8, 9.0, 'a)', fontsize=20, fontweight='bold', color='k')
    ax[0,0].text(-1.5, 9.2, '0819-0821 UTC', fontsize=16, fontweight='bold')
    ax[0,0].text(-0.4, 0.4, 'RaXPol', fontsize=16, fontweight='bold')
    ax[0,0].legend([s1,s2], ['Probe 1', 'Probe 2'], loc='lower left', fontsize=12)
    
    plot_cfill(vol[13]['xx'][eli,:,:], vol[13]['yy'][eli,:,:], np.ma.masked_array(vol[13]['vel'][eli,:,:], vol[13]['dbz'][eli,:,:]<1), 'vel', ax[0,1], datalims=[-30,30], xlims=xl, ylims=yl, cmap='balance', cbar=False)
    # ax[0,0].scatter(x_rot, y_rot, s=30, c='k', marker='.')
    # ax[0,1].set_title(f"{vol[13]['scan_time'][eli]} UTC {vol[9]['elev'][eli].round(1)}\N{DEGREE SIGN} radial velocity", fontsize=12)
    # b1 = ax[0,1].barbs(P1_x[::sep][i2], P1_y[::sep][i2]+3.3, P1_u[::sep][i2], P1_v[::sep][i2], barbcolor='k', length=7)
    # b2 = ax[0,1].barbs(P2_x[::sep][i2], P2_y[::sep][i2]+3.3, P2_u[::sep][i2], P2_v[::sep][i2], barbcolor='k', length=7)
    # ax[0,1].quiver(P1_x[::sep][i2], P1_y[::sep][i2]+3.3, P1_u[::sep][i2], P1_v[::sep][i2], color='k', scale=75, width=0.008, pivot='tail')
    # ax[0,1].quiver(P2_x[::sep][i2], P2_y[::sep][i2]+3.3, P2_u[::sep][i2], P2_v[::sep][i2], color='k', scale=75, width=0.008, pivot='tail')
    ax[0,1].scatter(P1_x[::sep][i2], P1_y[::sep][i2]+3.3, s=60, c=P1_c[::sep][i2], cmap=cmap,
                    vmin=datalims[0], vmax=datalims[1], marker='s', edgecolors='k', linewidth=0.75)
    ax[0,1].scatter(P2_x[::sep][i2], P2_y[::sep][i2]+3.3, s=60, c=P2_c[::sep][i2], cmap=cmap,
                    vmin=datalims[0], vmax=datalims[1], marker='^', edgecolors='k', linewidth=0.75)
    # ax[0,1].quiver(P1_x[::sep][i2], P1_y[::sep][i2]+3.3, P1_u[::sep][i2], P1_v[::sep][i2], color='k', scale=75, width=0.008, pivot='tail')
    # ax[0,1].quiver(P2_x[::sep][i2], P2_y[::sep][i2]+3.3, P2_u[::sep][i2], P2_v[::sep][i2], color='k', scale=75, width=0.008, pivot='tail')
    ax[0,1].scatter(0, 0.15, s=50, c='k')
    ax[0,1].text(-4.8, 9.0, 'b)', fontsize=20, fontweight='bold', color='k')
    ax[0,1].text(-1.5, 9.2, '0821-0823 UTC', fontsize=16, fontweight='bold')
    ax[0,1].text(-0.4, 0.4, 'RaXPol', fontsize=16, fontweight='bold')
    
    plot_cfill(vol[16]['xx'][eli,:,:], vol[16]['yy'][eli,:,:], np.ma.masked_array(vol[16]['vel'][eli,:,:], vol[16]['dbz'][eli,:,:]<1), 'vel', ax[1,0], datalims=[-30,30], xlims=xl, ylims=yl, cmap='balance', cbar=False)
    # ax[0,0].scatter(x_rot, y_rot, s=30, c='k', marker='.')
    # ax[1,0].set_title(f"{vol[16]['scan_time'][eli]} UTC {vol[9]['elev'][eli].round(1)}\N{DEGREE SIGN} radial velocity", fontsize=12)
    ax[1,0].set_xlabel('E-W distance from radar (km)', fontsize=14)
    ax[1,0].set_ylabel('N-S distance from radar (km)', fontsize=14)
    # b1 = ax[1,0].barbs(P1_x[::sep][i3], P1_y[::sep][i3]+3.3, P1_u[::sep][i3], P1_v[::sep][i3], barbcolor='k', length=7)
    # b2 = ax[1,0].barbs(P2_x[::sep][i3], P2_y[::sep][i3]+3.3, P2_u[::sep][i3], P2_v[::sep][i3], barbcolor='k', length=7)
    # ax[1,0].quiver(P1_x[::sep][i3], P1_y[::sep][i3]+3.3, P1_u[::sep][i3], P1_v[::sep][i3], color='k', scale=75, width=0.008, pivot='tail')
    # ax[1,0].quiver(P2_x[::sep][i3], P2_y[::sep][i3]+3.3, P2_u[::sep][i3], P2_v[::sep][i3], color='k', scale=75, width=0.008, pivot='tail')
    ax[1,0].scatter(P1_x[::sep][i3], P1_y[::sep][i3]+3.3, s=60, c=P1_c[::sep][i3], cmap=cmap,
                    vmin=datalims[0], vmax=datalims[1], marker='s', edgecolors='k', linewidth=0.75)
    ax[1,0].scatter(P2_x[::sep][i3], P2_y[::sep][i3]+3.3, s=60, c=P2_c[::sep][i3], cmap=cmap,
                    vmin=datalims[0], vmax=datalims[1], marker='^', edgecolors='k', linewidth=0.75)
    # ax[1,0].quiver(P1_x[::sep][i3], P1_y[::sep][i3]+3.3, P1_u[::sep][i3], P1_v[::sep][i3], color='k', scale=75, width=0.008, pivot='tail')
    # ax[1,0].quiver(P2_x[::sep][i3], P2_y[::sep][i3]+3.3, P2_u[::sep][i3], P2_v[::sep][i3], color='k', scale=75, width=0.008, pivot='tail')
    ax[1,0].scatter(0, 0.15, s=50, c='k')
    ax[1,0].text(-4.8, 9.0, 'c)', fontsize=20, fontweight='bold', color='k')
    ax[1,0].text(-1.5, 9.2, '0823-0825 UTC', fontsize=16, fontweight='bold')
    ax[1,0].text(-0.4, 0.4, 'RaXPol', fontsize=16, fontweight='bold')
    
    plot_cfill(vol[13]['xx'][eli,:,:], vol[13]['yy'][eli,:,:], np.ma.masked_array(vol[13]['vel'][eli,:,:], vol[13]['dbz'][eli,:,:]<1), 'vel', ax[1,1], datalims=[-30,30], xlims=xl, ylims=yl, cmap='balance', cbar=False)
    # ax[0,0].scatter(x_rot, y_rot, s=30, c='k', marker='.')
    # ax[1,1].set_title(f"{vol[13]['scan_time'][eli]} UTC {vol[9]['elev'][eli].round(1)}\N{DEGREE SIGN} radial velocity", fontsize=12)
    ax[1,1].set_xlabel('E-W distance from radar (km)', fontsize=14)
    # b1 = ax[1,1].barbs(P1_x[::sep], P1_y[::sep]+3.3, P1_u[::sep], P1_v[::sep], barbcolor='k', length=7)
    # b2 = ax[1,1].barbs(P2_x[::sep], P2_y[::sep]+3.3, P2_u[::sep], P2_v[::sep], barbcolor='k', length=7)
    ax[1,1].scatter(P1_x[::sep], P1_y[::sep]+3.3, s=60, c=P1_c[::sep], cmap=cmap,
                    vmin=datalims[0], vmax=datalims[1], marker='s', edgecolors='k', linewidth=0.75)
    ax[1,1].scatter(P2_x[::sep], P2_y[::sep]+3.3, s=60, c=P2_c[::sep], cmap=cmap,
                    vmin=datalims[0], vmax=datalims[1], marker='^', edgecolors='k', linewidth=0.75)
    # ax[1,1].scatter(P1_x, P1_y+3.3, s=40, c=P1_c, cmap=cmap, vmin=datalims[0], vmax=datalims[1], marker='.')
    # ax[1,1].scatter(P2_x, P2_y+3.3, s=40, c=P2_c, cmap=cmap, vmin=datalims[0], vmax=datalims[1], marker='.')
    # ax[1,1].quiver(P1_x[::sep], P1_y[::sep]+3.3, P1_u[::sep], P1_v[::sep], color='k', scale=100, width=0.008, pivot='tail')
    # ax[1,1].quiver(P2_x[::sep], P2_y[::sep]+3.3, P2_u[::sep], P2_v[::sep], color='k', scale=100, width=0.008, pivot='tail')
    ax[1,1].quiver(P1_x[::sep], P1_y[::sep]+3.3, u_P1, v_P1, color='k', scale=100, width=0.008, pivot='tail')
    ax[1,1].quiver(P2_x[::sep], P2_y[::sep]+3.3, u_P2, v_P2, color='k', scale=100, width=0.008, pivot='tail')
    ax[1,1].quiver(-4.5, 0.5, 10, 0, color='k', scale=100, width=0.008, pivot='tail')
    ax[1,1].text(-4.5, 0.8, '10 m/s', fontsize=12, fontweight='bold')
    ax[1,1].scatter(0, 0.15, s=50, c='k')
    ax[1,1].text(-4.8, 9.0, 'd)', fontsize=20, fontweight='bold', color='k')
    # ax[1,1].text(-1.5, 9.2, '0819-0825 UTC', fontsize=16, fontweight='bold')
    ax[1,1].text(-0.2, 9.2, 'Full circuit', fontsize=18, fontweight='bold')
    ax[1,1].text(-0.4, 0.4, 'RaXPol', fontsize=16, fontweight='bold')
    
    cb1 = plt.colorbar(c, ax=[ax[1,0],ax[1,1]], extend='both', orientation='horizontal', aspect=30)
    cb1.set_label("RaXPol velocity (m s$^{-1}$)", fontsize=14)
    cb1.ax.tick_params(labelsize=12)
    cb2 = plt.colorbar(s1, ax=[ax[0,1],ax[1,1]], extend='both', aspect=30)
    # cb2.set_label("Temperature ($^{\circ}$C)", fontsize=16)
    cb2.set_label("Mesonet \u03B8'$_v$ (K)", fontsize=16)
    cb2.set_ticks(np.arange(-1, 1.2, 0.2))
    cb2.ax.tick_params(labelsize=12)
    
    if figsave:
        plt.savefig(ip+'circuit_thvpert.png', dpi=300)


# Very rough estimate of divergence at the NTV
u1 = P1_u[::sep][4] # i1[3]
v1 = P1_v[::sep][4]
x1 = P1_x[::sep][4] * 1000
y1 = P1_y[::sep][4] * 1000
u2 = P2_u[::sep][13] # i3[2]
v2 = P2_v[::sep][13]
x2 = P2_x[::sep][13] * 1000
y2 = P2_y[::sep][13] * 1000

dx = x1-x2; dy = y1-y2; du = u1-u2; dv = v1-v2

div = du/dx + dv/dy


#%%

fig,ax = plt.subplots(1, 1, figsize=(8,6))
l1, = ax.plot(P1_times[14253:14728], P1['lat'][14253:14728], 'k')
l2, = ax.plot(P2_times[11784:12261], P2['lat'][11784:12261], 'b')
# ax.scatter(P1_times[14255], P1['lat'][14255], s=20, c='k') # 14308 for 081854
# ax.scatter(P2_times[11786], P2['lat'][11786], s=20, c='b') # 11840 for 081854
# ax.text(P2_times[11786], P2['lat'][11786]+0.001, '081800 UTC') # 11840 for 081854
plt.legend(handles=[l1,l2], labels=['P1','P2'])
ax.set_ylabel('Latitude')
ax.set_xlabel('Time')
ax.set_title('Mesonet latitude')


fig,ax = plt.subplots(1, 1, figsize=(8,6))
l1, = ax.plot(P1_times[14253:14728], P1['temp_unbiased'][14253:14728], 'k')
l2, = ax.plot(P2_times[11784:12261], P2['temp_unbiased'][11784:12261], 'b')
# ax.scatter(P1_times[14255], P1['temp_unbiased'][14255], s=20, c='k') # 14308 for 081854
# ax.scatter(P2_times[11786], P2['temp_unbiased'][11786], s=20, c='b') # 11840 for 081854
# ax.text(P2_times[11786], P2['temp_unbiased'][11786]+0.04, '081800 UTC') # 11840 for 081854
plt.legend(handles=[l1,l2], labels=['P1','P2'])
ax.set_ylabel('Temp_unbiased')
ax.set_xlabel('Time')
ax.set_title('Mesonet unbiased temperature')


fig,ax = plt.subplots(1, 1, figsize=(8,6))
l1, = ax.plot(P1_times[14253:14728], P1['rh_unbiased'][14253:14728], 'k')
l2, = ax.plot(P2_times[11784:12261], P2['rh_unbiased'][11784:12261], 'b')
# ax.scatter(P1_times[14255], P1['rh_unbiased'][14255], s=20, c='k') # 14308 for 081854
# ax.scatter(P2_times[11786], P2['rh_unbiased'][11786], s=20, c='b') # 11840 for 081854
# ax.text(P2_times[11786], P2['rh_unbiased'][11786], '081800 UTC') # 11840 for 081854
plt.legend(handles=[l1,l2], labels=['P1','P2'])
ax.set_ylabel('RH_unbiased')
ax.set_xlabel('Time')
ax.set_title('Mesonet unbiased relative humidity')


fig,ax = plt.subplots(1, 1, figsize=(8,6))
l1, = ax.plot(P1_times[14253:14728], P1['pres_unbiased'][14253:14728], 'k')
l2, = ax.plot(P2_times[11784:12261], P2['pres_unbiased'][11784:12261], 'b')
# ax.scatter(P1_times[14255], P1['pres_unbiased'][14255], s=20, c='k') # 14308 for 081854
# ax.scatter(P2_times[11786], P2['pres_unbiased'][11786], s=20, c='b') # 11840 for 081854
# ax.text(P2_times[11786], P2['pres_unbiased'][11786], '081800 UTC') # 11840 for 081854
plt.legend(handles=[l1,l2], labels=['P1','P2'])
ax.set_ylabel('Pres_unbiased')
ax.set_xlabel('Time')
ax.set_title('Mesonet unbiased pressure')


fig,ax = plt.subplots(1, 1, figsize=(8,6))
l1, = ax.plot(P1_times[14253:14728], P1['u_corr'][14253:14728], 'k')
l2, = ax.plot(P2_times[11784:12261], P2['u_corr'][11784:12261], 'b')
# ax.scatter(P1_times[14255], P1['u_corr'][14255], s=20, c='k') # 14308 for 081854
# ax.scatter(P2_times[11786], P2['u_corr'][11786], s=20, c='b') # 11840 for 081854
# ax.text(P2_times[11786], P2['u_corr'][11786], '081800 UTC') # 11840 for 081854
plt.legend(handles=[l1,l2], labels=['P1','P2'])
ax.set_ylabel('U_corrected')
ax.set_xlabel('Time')
ax.set_title('Mesonet corrected u wind')


fig,ax = plt.subplots(1, 1, figsize=(8,6))
l1, = ax.plot(P1_times[14253:14728], P1['v_corr'][14253:14728], 'k')
l2, = ax.plot(P2_times[11784:12261], P2['v_corr'][11784:12261], 'b')
# ax.scatter(P1_times[14255], P1['v_corr'][14255], s=20, c='k') # 14308 for 081854
# ax.scatter(P2_times[11786], P2['v_corr'][11786], s=20, c='b') # 11840 for 081854
# ax.text(P2_times[11786], P2['v_corr'][11786], '081800 UTC') # 11840 for 081854
plt.legend(handles=[l1,l2], labels=['P1','P2'])
ax.set_ylabel('V_corrected')
ax.set_xlabel('Time')
ax.set_title('Mesonet corrected v wind')


fig,ax = plt.subplots(1, 1, figsize=(8,6))
l1, = ax.plot(P1_times[14253:14728], P1['wspd_corr'][14253:14728], 'k')
l2, = ax.plot(P2_times[11784:12261], P2['wspd_corr'][11784:12261], 'b')
# ax.scatter(P1_times[14255], P1['wspd_corr'][14255], s=20, c='k') # 14308 for 081854
# ax.scatter(P2_times[11786], P2['wspd_corr'][11786], s=20, c='b') # 11840 for 081854
# ax.text(P2_times[11786], P2['wspd_corr'][11786], '081800 UTC') # 11840 for 081854
plt.legend(handles=[l1,l2], labels=['P1','P2'])
ax.set_ylabel('Wspd_corrected')
ax.set_xlabel('Time')
ax.set_title('Mesonet corrected wind speed')

plt.show()



#%%

fig,ax = plt.subplots(1, 1, figsize=(8,6))
l1, = ax.plot(P1_times[13667:14256], P1['lat'][13667:14256], 'k')
l2, = ax.plot(P2_times[11187:11787], P2['lat'][11187:11787], 'b')
plt.legend(handles=[l1,l2], labels=['P1','P2'])
ax.set_ylabel('Latitude')
ax.set_xlabel('Time')
ax.set_title('Mesonet latitude')


fig,ax = plt.subplots(1, 1, figsize=(8,6))
l1, = ax.plot(P1_times[13667:14256], P1['temp_unbiased'][13667:14256], 'k')
l2, = ax.plot(P2_times[11187:11787], P2['temp_unbiased'][11187:11787], 'b')
plt.legend(handles=[l1,l2], labels=['P1','P2'])
ax.set_ylabel('Temp_unbiased')
ax.set_xlabel('Time')
ax.set_title('Mesonet unbiased temperature')


fig,ax = plt.subplots(1, 1, figsize=(8,6))
l1, = ax.plot(P1_times[13667:14256], P1['rh_unbiased'][13667:14256], 'k')
l2, = ax.plot(P2_times[11187:11787], P2['rh_unbiased'][11187:11787], 'b')
plt.legend(handles=[l1,l2], labels=['P1','P2'])
ax.set_ylabel('RH_unbiased')
ax.set_xlabel('Time')
ax.set_title('Mesonet unbiased relative humidity')


fig,ax = plt.subplots(1, 1, figsize=(8,6))
l1, = ax.plot(P1_times[13667:14256], P1['pres_unbiased'][13667:14256], 'k')
l2, = ax.plot(P2_times[11187:11787], P2['pres_unbiased'][11187:11787], 'b')
plt.legend(handles=[l1,l2], labels=['P1','P2'])
ax.set_ylabel('Pres_unbiased')
ax.set_xlabel('Time')
ax.set_title('Mesonet unbiased pressure')


fig,ax = plt.subplots(1, 1, figsize=(8,6))
l1, = ax.plot(P1_times[13667:14256], movmean(P1['u_corr'][13667:14256],10), 'k')
l2, = ax.plot(P2_times[11187:11787], movmean(P2['u_corr'][11187:11787],10), 'b')
plt.legend(handles=[l1,l2], labels=['P1','P2'])
ax.set_ylabel('U_corrected')
ax.set_xlabel('Time')
ax.set_title('Mesonet corrected u wind')


fig,ax = plt.subplots(1, 1, figsize=(8,6))
l1, = ax.plot(P1_times[13667:14256], movmean(P1['v_corr'][13667:14256],10), 'k')
l2, = ax.plot(P2_times[11187:11787], movmean(P2['v_corr'][11187:11787],10), 'b')
plt.legend(handles=[l1,l2], labels=['P1','P2'])
ax.set_ylabel('V_corrected')
ax.set_xlabel('Time')
ax.set_title('Mesonet corrected v wind')


fig,ax = plt.subplots(1, 1, figsize=(8,6))
l1, = ax.plot(P1_times[13667:14256], movmean(P1['wspd_corr'][13667:14256],10), 'k')
l2, = ax.plot(P2_times[11187:11787], movmean(P2['wspd_corr'][11187:11787],10), 'b')
plt.legend(handles=[l1,l2], labels=['P1','P2'])
ax.set_ylabel('Wspd_corrected')
ax.set_xlabel('Time')
ax.set_title('Mesonet corrected wind speed')

plt.show()














