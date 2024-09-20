# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 12:34:41 2023
@author: morgan.schneider

Plotting RaXPol data
"""

####################
### Load modules ###
####################

from RaXPolUtils import *
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#%% IOP2 single time with mobile mesonet

figsave = False

fp = 'Documents/perils2023/iop2/raxpol/CFradial/'
filetime = '082530'
fn = glob(fp+f"*_{filetime}_*.nc")[0]

rax = read_raxpol(fn)
raxpol = pyart.io.read(fn)

P1 = read_MM('Documents/perils2023/iop2/mesonet/Probe_1_IOP2_QC_all.dat')
P2 = read_MM('Documents/perils2023/iop2/mesonet/Probe_2_IOP2_QC_all.dat')
xx1,yy1 = latlon2xy(P1.lat, P1.lon, raxpol.latitude['data'][0], raxpol.longitude['data'][0])
xx2,yy2 = latlon2xy(P2.lat, P2.lon, raxpol.latitude['data'][0], raxpol.longitude['data'][0])
P1.SetAttr('xx', xx1)
P1.SetAttr('yy', yy1)
P2.SetAttr('xx', xx2)
P2.SetAttr('yy', yy2)

elev = np.mean(raxpol.elevation['data']).round(1)

i1 = np.where(P1.time == int(filetime))[0][0]
i2 = np.where(P2.time == int(filetime))[0][0]

T_lims = [18,21] # Temperature
Td_lims = [17,19] # Dewpoint
datalims = T_lims


c = pyart.graph.RadarDisplay(raxpol)

fig = plt.figure(figsize=(9,6))
ax = fig.add_subplot(111)
c.plot('DBZ', 0, vmin=0, vmax=70, cmap='pyart_NWSRef')
c.plot_range_rings([2.5, 5, 7.5, 10])
ax.set_xlim([-10,10])
ax.set_ylim([-10,10])
s1 = ax.scatter(P1.xx[i1], P1.yy[i1], s=30, c=P1.Utube[i1], cmap=cmaps['temp']['cm'],
                vmin=datalims[0], vmax=datalims[1], marker='s', edgecolors='k')
s2 = ax.scatter(P2.xx[i2], P2.yy[i2], s=30, c=P2.Utube[i2], cmap=cmaps['temp']['cm'],
                vmin=datalims[0], vmax=datalims[1], marker='^', edgecolors='k')
ax.scatter(0, 0, s=30, c='k')
plt.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')
plt.colorbar(s1,label='Sfc temperature (C)')
plt.legend([s1,s2], ['Probe 1', 'Probe 2'], loc='upper right')
plt.show()

if figsave:
    plt.savefig(f"Documents/perils2023/iop2/figs/dbz+T_{filetime}.png", dpi=400)


figsave = False

# fig = plt.figure(figsize=(8,6))
# ax = fig.add_subplot(111)
# plot_cfill(rax.xx, rax.yy, rax.vel, 'vel', ax, datalims=[-rax.va,rax.va], xlims=[-8,8], ylims=[-8,8])
# # s1 = ax.scatter(P1.xx[i1], P1.yy[i1], s=30, c='k', marker='s')
# # s2 = ax.scatter(P2.xx[i2], P2.yy[i2], s=30, c='k', marker='^')
# # ax.text(P1.xx[i1]+0.3, P1.yy[i1], 'P1', fontsize=12, fontweight='bold')
# # ax.text(P2.xx[i2]+0.3, P2.yy[i2], 'P2', fontsize=12, fontweight='bold')
# # ax.scatter(0, 0, s=30, c='k')
# # ax.text(-1, 0.3, 'RaXPol', fontsize=12, fontweight='bold')
# # plt.legend([s1,s2], ['Probe 1', 'Probe 2'], loc='upper right')
# ax.set_xlabel('E-W distance from radar (km)', fontsize=12)
# ax.set_ylabel('N-S distance from radar (km)', fontsize=12)
# ax.set_title(f"{rax.elev:.1f}\N{DEGREE SIGN} radial velocity - {filetime} UTC", fontsize=14, fontweight='bold')
# # plt.suptitle(f"{vol[vi].elev[eli].round(1)}\N{DEGREE SIGN} radial velocity - {filetime} UTC")
# if figsave:
#     plt.savefig(f"Documents/conferences/perilsmeeting/vol{vi}_{filetime}_az{azimuth}_PPI_V_noAz.png", dpi=400)


fig,(ax1,ax2) = plt.subplots(1,2,figsize=(14,5))
plot_cfill(rax.xx, rax.yy, rax.dbz, 'dbz', ax1, datalims=[0,70], xlims=[-rlim,rlim], ylims=[-rlim,rlim])
# ax1.plot(vol[vi].xx[eli,azi,:], vol[vi].yy[eli,azi,:], '--k', linewidth=1)
# s1 = ax1.scatter(P1.xx[i1], P1.yy[i1], s=50, c=P1.Utube[i1], cmap=cmaps['temp']['cm'],
#                 vmin=datalims[0], vmax=datalims[1], marker='s', edgecolors='k')
# s2 = ax1.scatter(P2.xx[i2], P2.yy[i2], s=50, c=P2.Utube[i2], cmap=cmaps['temp']['cm'],
#                 vmin=datalims[0], vmax=datalims[1], marker='^', edgecolors='k')
# ax1.scatter(0, 0, s=50, c='k')
# ax1.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')
# plt.colorbar(s1,label='Sfc temperature (C)')
# plt.legend([s1,s2], ['Probe 1', 'Probe 2'], loc='upper right')
ax1.set_title(f"{rax.elev.round(1)}\N{DEGREE SIGN} reflectivity", fontsize=14)
ax1.set_xlabel('E-W distance from radar (km)', fontsize=12)
ax1.set_ylabel('N-S distance from radar (km)', fontsize=12)

plot_cfill(vol[vi].xx[eli,:,:], vol[vi].yy[eli,:,:], vol[vi].vel[eli,:,:], 'vel', ax2, datalims=[-va,va], xlims=[-rlim,rlim], ylims=[-rlim,rlim])
# ax2.plot(vol[vi].xx[eli,azi,:], vol[vi].yy[eli,azi,:], '--k', linewidth=1)
ax2.set_title(f"{rax.elev.round(1)}\N{DEGREE SIGN} radial velocity", fontsize=14)
ax2.set_xlabel('E-W distance from radar (km)', fontsize=12)
ax2.set_ylabel('N-S distance from radar (km)', fontsize=12)

# plt.suptitle(f"{filetime} UTC, azimuth = {azimuth}\N{DEGREE SIGN}", fontsize=14)
plt.suptitle(f"{filetime} UTC", fontsize=14)
if figsave:
    # plt.savefig(f"Documents/perils2023/iop2/figs/rhis/vol{vi}_{filetime}_az{azimuth}_PPI.png", dpi=400)
    plt.savefig(f"Documents/conferences/perilsmeeting/vol{vi}_{filetime}_az{azimuth}_PPI.png", dpi=400)


#%% IOP5 single time with skyler

# fp = 'Documents/perils2023/iop5/raxpol/CFradial/'
# filetime = '160606'
# fn = glob(fp+f"*_{filetime}_*.nc")[0]

# 163100 - 1065
# 163500 - 1140

fp = 'Documents/perils2023/iop5/raxpol/CFradial/'
files = sorted(glob(fp+'*.nc'))
# fn = files[1140]
fn = fp+'cfrad.20230405_160705_RaXPol_v104_s4.nc'


raxpol = pyart.io.read(fn)

sky_lat = 35.4110
sky_lon = -90.9549
sky_xx,sky_yy = latlon2xy(sky_lat, sky_lon, raxpol.latitude['data'][0], raxpol.longitude['data'][0])

c = pyart.graph.RadarDisplay(raxpol)

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111)
c.plot('DBZ', 0, vmin=0, vmax=70, cmap='pyart_NWSRef')
c.plot_range_rings([10,20,30])
ax.set_xlim([-30,30])
ax.set_ylim([-30,30])
ax.scatter(0, 0, s=50, c='k')
ax.scatter(sky_xx, sky_yy, s=50, c='k')
plt.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')
plt.text(sky_xx-1, sky_yy+0.4, 'Skyler', fontsize=10, fontweight='bold')
plt.show()


fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111)
c.plot('VEL', 0, vmin=-30, vmax=30, cmap='pyart_Carbone42')
c.plot_range_rings([10,20,30])
ax.set_xlim([-30,30])
ax.set_ylim([-30,30])
ax.scatter(0, 0, s=50, c='k')
ax.scatter(sky_xx, sky_yy, s=50, c='k')
plt.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')
plt.text(sky_xx-1, sky_yy+0.4, 'Skyler', fontsize=10, fontweight='bold')
plt.show()






#%% IOP2 load data into volume structs

# from RaXPolUtils import *

# Pre-transect:       File 381 (index 380) 0.9 deg, 081758 UTC
# Transect start:     File 411 (index 410) 0.9 deg, 081858 UTC
# Outflow at P1:      File 462 (index 461) 2.6 deg, 082100 UTC
# Reflectivity at P1: File 495 (index 494) 0.9 deg, 082228 UTC
# Transect end:       File 600 (index 599) 0.9 deg, 082558 UTC
# Post-transect:      File (index ) 0.9 deg, 082658 UTC

vol_nums = [78, 81, 84, 87, 90, 93, 96, 102, 105, 108, 112, 115, 118, 121, 124, 127, 130, 133, 136, 139]
# v78 320-331 # 081558-081620
# v81 335-346 # 081628-081650
# v84 350-361 # 081658-081720
# v87 365-376 # 081728-181750
# v90 380-391 # 081758-081820
# v93 395-406 # 081828-181850
# v96 410-421 # 081858-181920 *Transect start - vol[6]
# v102 430-441 # 081958-082020
# v105 445-456 # 082028-082050
# v108 460-471 # 082058-082120 *Outflow at P1 - vol[9]
# v112 480-490 # 082158-082220
# v115 494-505 # 082228-082250 *Reflectivity at P1 - vol[11]
# v118 509-520 # 082258-082320
# v121 524-535 # 082328-082350
# v124 539-550 # 082358-082420
# v127 554-565 # 082428-082450
# v130 569-580 # 082458-082520
# v133 584-595 # 082528-082550
# v136 599-610 # 082558-082620 *Transect end - vol[18]
# v139 614-625 # 082628-082650

t_start = '081758'
transect_start = '081858'
outflow_at_MM = '082100'
precip_at_MM = '082228'
transect_end = '082558'
t_end = '082658'

fp = 'Documents/perils2023/iop2/raxpol/CFradial/'
files = sorted(glob(fp+'*.nc'))

vol = [struct() for i in range(len(vol_nums))]

for vn in range(len(vol_nums)):
    inds = [i for i,s in enumerate(files) if f"v{vol_nums[vn]}" in s]
    a = slice(inds[1], inds[-1]+2)
    f = files[a]
    
    dbz_tmp = np.zeros(shape=(len(f),360,1246,), dtype=float)
    vel_tmp = np.zeros(shape=(len(f),360,1246,), dtype=float)
    sw_tmp = np.zeros(shape=(len(f),360,1246,), dtype=float)
    zdr_tmp = np.zeros(shape=(len(f),360,1246,), dtype=float)
    rhohv_tmp = np.zeros(shape=(len(f),360,1246,), dtype=float)
    xx_tmp = np.zeros(shape=(len(f),360,1246,), dtype=float)
    yy_tmp = np.zeros(shape=(len(f),360,1246,), dtype=float)
    zz_tmp = np.zeros(shape=(len(f),360,1246,), dtype=float)
    el_tmp = np.zeros(shape=(len(f),), dtype=float)
    az_tmp = np.zeros(shape=(len(f),360,), dtype=float)
    time_tmp = []
    fname_tmp = []
    
    for ii in range(len(f)):
        fn = f[ii]
        d = read_raxpol(fn)
        
        r = d.r
        el_tmp[ii] = d.elev
        time_tmp.append(fn[57:63])
        fname_tmp.append(fn)
        
        for ix in np.linspace(0,359,360):
            if ix in d.az.round(0):
                ind = np.where(d.az.round(0) == ix)[0][0]
                
                az_tmp[ii,int(ix)] = d.az[ind]
                dbz_tmp[ii,int(ix),:] = d.dbz[ind,:]
                vel_tmp[ii,int(ix),:] = d.vel[ind,:]
                sw_tmp[ii,int(ix),:] = d.sw[ind,:]
                zdr_tmp[ii,int(ix),:] = d.zdr[ind,:]
                rhohv_tmp[ii,int(ix),:] = d.rhohv[ind,:]
                xx_tmp[ii,int(ix),:] = d.xx[ind,:]
                yy_tmp[ii,int(ix),:] = d.yy[ind,:]
                zz_tmp[ii,int(ix),:] = d.zz[ind,:]
            else:
                az_tmp[ii,int(ix)] = ix
                xx_tmp[ii,int(ix),:] = r * np.sin(ix*np.pi/180) * np.cos(d.elev*np.pi/180)
                yy_tmp[ii,int(ix),:] = r * np.cos(ix*np.pi/180) * np.cos(d.elev*np.pi/180)
                zz_tmp[ii,int(ix),:] = r * np.sin(d.elev*np.pi/180)
    
    vol[vn].SetAttr('dbz', dbz_tmp)
    vol[vn].SetAttr('vel', vel_tmp)
    vol[vn].SetAttr('sw', sw_tmp)
    vol[vn].SetAttr('zdr', zdr_tmp)
    vol[vn].SetAttr('rhohv', rhohv_tmp)
    vol[vn].SetAttr('xx', xx_tmp)
    vol[vn].SetAttr('yy', yy_tmp)
    vol[vn].SetAttr('zz', zz_tmp)
    vol[vn].SetAttr('az', az_tmp)
    vol[vn].SetAttr('elev', el_tmp)
    vol[vn].SetAttr('scan_time', time_tmp)
    vol[vn].SetAttr('vol_num', vol_nums[vn])
    vol[vn].SetAttr('filename', fname_tmp)
        
r = d.r
rax_lat = d.lat
rax_lon = d.lon
va = d.va

P1 = read_MM('Documents/perils2023/iop2/mesonet/Probe_1_IOP2_QC_all.dat')
P2 = read_MM('Documents/perils2023/iop2/mesonet/Probe_2_IOP2_QC_all.dat')
xx1,yy1 = latlon2xy(P1.lat, P1.lon, rax_lat, rax_lon)
xx2,yy2 = latlon2xy(P2.lat, P2.lon, rax_lat, rax_lon)
P1.SetAttr('xx', xx1)
P1.SetAttr('yy', yy1)
P2.SetAttr('xx', xx2)
P2.SetAttr('yy', yy2)


#%% IOP2 plot reconstructed RHIs

# Scatter plot of rotor(s) location on top of PPIs for each volume
# Azimuth-height pcolor plot of velocity through MV/rotor
# Do I need to advection correct for the reconstructed RHIs?? Check storm motion
# If couplet translation is due to advection, storm motion is about 15 m/s to the ~NE

# plt.close('all')
figsave = False

rlim = 8
zlim = 2.5

vi = 7
eli = 1
filetime = vol[vi].scan_time[eli]
azimuth = 305
azi = np.where(vol[ii].az[eli,:].round(0) == azimuth)[0][0]
rr = (vol[vi].xx[:,azi,:]**2 + vol[vi].yy[:,azi,:]**2)**0.5

# i1 = np.where(P1.time == int(filetime))[0][0]
# i2 = np.where(P2.time == int(filetime))[0][0]

T_lims = [18,21] # Temperature
Td_lims = [17,19] # Dewpoint
datalims = T_lims

az_rot = np.array([308., 309., 310., 311., 312., 
          313., 314., 315., 316., 317., 
          318., 319., 320., 321., 322., 
          323., 324., 325., 326., 327., 
          328., 329., 330.])
r_rot = np.array([1.70, 1.70, 1.69, 1.68, 1.68,
          1.68, 1.68, 1.67, 1.68, 1.68,
          1.67, 1.68, 1.68, 1.68, 1.69,
          1.70, 1.70, 1.72, 1.72, 1.72,
          1.73, 1.74, 1.74])
z_rot = np.array([0.03, 0.04, 0.06, 0.08, 0.13,
          0.16, 0.17, 0.18, 0.21, 0.22,
          0.23, 0.26, 0.30, 0.33, 0.39,
          0.41, 0.44, 0.46, 0.49, 0.53,
          0.54, 0.56, 0.60])

x_rot = r_rot * np.sin(az_rot*np.pi/180)
y_rot = r_rot * np.cos(az_rot*np.pi/180)

# irot = np.where(np.isclose(az_rot, azimuth))[0][0]

if True:
    fig,(ax1,ax2) = plt.subplots(1,2,figsize=(14,5))
    
    plot_cfill(vol[vi].xx[eli,:,:], vol[vi].yy[eli,:,:], vol[vi].dbz[eli,:,:], 'dbz', ax1, datalims=[0,70], xlims=[-rlim/2,0], ylims=[0,rlim/2])
    # ax1.scatter(x_rot, y_rot, s=5, c='k', marker='.')
    ax1.set_title(f"{vol[vi].elev[eli].round(1)}\N{DEGREE SIGN} reflectivity", fontsize=14)
    ax1.set_xlabel('E-W distance from radar (km)', fontsize=12)
    ax1.set_ylabel('N-S distance from radar (km)', fontsize=12)
    # s1 = ax1.scatter(P1.xx[i1], P1.yy[i1], s=50, c=P1.Utube[i1], cmap=cmaps['temp']['cm'],
    #                 vmin=datalims[0], vmax=datalims[1], marker='s', edgecolors='k')
    # s2 = ax1.scatter(P2.xx[i2], P2.yy[i2], s=50, c=P2.Utube[i2], cmap=cmaps['temp']['cm'],
    #                 vmin=datalims[0], vmax=datalims[1], marker='^', edgecolors='k')
    # ax1.scatter(0, 0, s=50, c='k')
    # ax1.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')
    # plt.colorbar(s1,label='Sfc temperature (C)')
    # plt.legend([s1,s2], ['Probe 1', 'Probe 2'], loc='upper right')
    
    plot_cfill(vol[vi].xx[eli,:,:], vol[vi].yy[eli,:,:], vol[vi].vel[eli,:,:], 'vel', ax2, datalims=[-va,va], xlims=[-rlim/2,0], ylims=[0,rlim/2])
    # ax2.scatter(x_rot, y_rot, s=5, c='k', marker='.')
    ax2.set_title(f"{vol[vi].elev[eli].round(1)}\N{DEGREE SIGN} radial velocity", fontsize=14)
    ax2.set_xlabel('E-W distance from radar (km)', fontsize=12)
    ax2.set_ylabel('N-S distance from radar (km)', fontsize=12)
    # s1 = ax2.scatter(P1.xx[i1], P1.yy[i1], s=50, c=P1.Utube[i1], cmap=cmaps['temp']['cm'],
    #                 vmin=datalims[0], vmax=datalims[1], marker='s', edgecolors='k')
    # s2 = ax2.scatter(P2.xx[i2], P2.yy[i2], s=50, c=P2.Utube[i2], cmap=cmaps['temp']['cm'],
    #                 vmin=datalims[0], vmax=datalims[1], marker='^', edgecolors='k')
    # ax2.scatter(0, 0, s=50, c='k')
    # ax2.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')
    # plt.colorbar(s1,label='Sfc temperature (C)')
    # plt.legend([s1,s2], ['Probe 1', 'Probe 2'], loc='upper right')
    
    # plt.suptitle(f"{filetime} UTC, azimuth = {azimuth}\N{DEGREE SIGN}", fontsize=14)
    plt.suptitle(f"{filetime} UTC", fontsize=14)
    if figsave:
        plt.savefig(f"Documents/perils2023/iop2/figs/rhis/vol{vi}_{filetime}_az{azimuth}_PPI.png", dpi=400)
        # plt.savefig(f"Documents/conferences/perilsmeeting/iop2_vol{vi}_{filetime}_az{azimuth}_PPI.png", dpi=400)




if True:
    
    fig,((ax1),(ax2)) = plt.subplots(2,1,figsize=(12,9))
    
    plot_cfill(rr, vol[vi].zz[:,azi,:], vol[vi].dbz[:,azi,:], 'dbz', ax1, datalims=[0,70], xlims=[0,rlim], ylims=[0,zlim])
    # ax1.scatter(r_rot[irot], z_rot[irot], s=10, c='k', marker='x')
    ax1.set_ylabel('Height ARL (km)', fontsize=14)
    ax1.set_title(f"Azimuth = {azimuth}\N{DEGREE SIGN}\n\n Reflectivity - {filetime} UTC", fontsize=14, fontweight='bold')
    ax1.invert_xaxis()
    ax1_ppi = inset_axes(ax1, '22%', '50%', loc=1)
    c,cb1 = plot_cfill(vol[vi].xx[eli,:,:], vol[vi].yy[eli,:,:], vol[vi].dbz[0,:,:], 'dbz', ax1_ppi, datalims=[0,70], xlims=[-rlim,rlim], ylims=[-rlim,rlim])
    ax1_ppi.plot(vol[vi].xx[eli,azi,:], vol[vi].yy[eli,azi,:], '--k', linewidth=1.25)
    cb1.remove()
    ax1_ppi.set_xticks([])
    ax1_ppi.set_yticks([])
    
    plot_cfill(rr, vol[vi].zz[:,azi,:], vol[vi].vel[:,azi,:], 'vel', ax2, datalims=[-va,va], xlims=[0,rlim], ylims=[0,zlim])
    # ax2.scatter(r_rot[irot], z_rot[irot], s=10, c='k', marker='x')
    ax2.set_xlabel('Range from radar (km)', fontsize=14)
    ax2.set_ylabel('Height ARL (km)', fontsize=14)
    ax2.set_title(f"Radial velocity - {filetime} UTC", fontsize=14, fontweight='bold')
    ax2.invert_xaxis()
    ax2_ppi = inset_axes(ax2, '22%', '50%', loc=1)
    c,cb2 = plot_cfill(vol[vi].xx[eli,:,:], vol[vi].yy[eli,:,:], vol[vi].vel[eli,:,:], 'vel', ax2_ppi, datalims=[-va,va], xlims=[-rlim,rlim], ylims=[-rlim,rlim])
    ax2_ppi.plot(vol[vi].xx[eli,azi,:], vol[vi].yy[eli,azi,:], '--k', linewidth=1.25)
    cb2.remove()
    ax2_ppi.set_xticks([])
    ax2_ppi.set_yticks([])
    
    plt.draw()
    
    # plt.suptitle(f"{filetime} UTC, azimuth = {azimuth}\N{DEGREE SIGN}")
    if figsave:
        plt.savefig(f"Documents/perils2023/iop2/figs/rhis/vol{vi}_{filetime}_az{azimuth}_RHI.png", dpi=400)
        # plt.savefig(f"Documents/conferences/perilsmeeting/vol{vi}_{filetime}_az{azimuth}_RHI.png", dpi=400)


if True:
    azs = slice(285,321)
    ri = 82
    
    fig,(ax1,ax2) = plt.subplots(1,2,figsize=(14,5))
    
    plot_cfill(vol[vi].az[0,azs], vol[vi].zz[:,315,ri], vol[vi].dbz[:,azs,ri], 'dbz', ax1, datalims=[0,70], xlims=[285,320], ylims=[0,0.8])
    ax1.set_title(f"Reflectivity at r = {d.r[ri]:.3} km", fontsize=14)
    ax1.set_xlabel('Azimuth (deg)', fontsize=12)
    ax1.set_ylabel('Height ARL (km)', fontsize=12)
    
    plot_cfill(vol[vi].az[0,azs], vol[vi].zz[:,315,ri], vol[vi].vel[:,azs,ri], 'vel', ax2, datalims=[-va,va], xlims=[285,320], ylims=[0,0.8])
    ax2.set_title(f"Radial velocity at r = {d.r[ri]:.3} km", fontsize=14)
    ax2.set_xlabel('Azimuth (deg)', fontsize=12)
    ax2.set_ylabel('Height ARL (km)', fontsize=12)
    
    # plt.suptitle(f"{filetime} UTC, azimuth = {azimuth}\N{DEGREE SIGN}", fontsize=14)
    plt.suptitle(f"{filetime} UTC", fontsize=14)
    if figsave:
        plt.savefig(f"Documents/perils2023/iop2/figs/rhis/vol{vi}_{filetime}_az{azimuth}_RHI.png", dpi=400)
        # plt.savefig(f"Documents/conferences/perilsmeeting/vol{vi}_{filetime}_az{azimuth}_RHI.png", dpi=400)


# fig = plt.figure(figsize=(8,6))
# ax = fig.add_subplot(111)
# plot_cfill(vol[vi].xx[eli,:,:], vol[vi].yy[eli,:,:], vol[vi].dbz[0,:,:], 'dbz', ax, datalims=[0,70], xlims=[-rlim,rlim], ylims=[-rlim,rlim])
# ax.plot(vol[vi].xx[eli,azi,:], vol[vi].yy[eli,azi,:], '--k', linewidth=1.5)
# ax.set_xlabel('E-W distance from radar (km)')
# ax.set_ylabel('N-S distance from radar (km)')
# ax.set_title(f"{vol[vi].elev[eli].round(1)}\N{DEGREE SIGN} reflectivity - {filetime} UTC", fontsize=14)
# # plt.suptitle(f"{filetime} UTC PPI at {vol[vi].elev[eli].round(1)}\N{DEGREE SIGN}, azimuth = {azimuth}\N{DEGREE SIGN}")
# if figsave:
#     plt.savefig(f"Documents/conferences/perilsmeeting/vol{vi}_{filetime}_az{azimuth}_PPI_Z.png", dpi=400)


# fig = plt.figure(figsize=(8,6))
# ax = fig.add_subplot(111)
# plot_cfill(vol[vi].xx[eli,:,:], vol[vi].yy[eli,:,:], vol[vi].vel[eli,:,:], 'vel', ax, datalims=[-va,va], xlims=[-rlim,rlim], ylims=[-rlim,rlim])
# ax.plot(vol[vi].xx[eli,azi,:], vol[vi].yy[eli,azi,:], '--k', linewidth=1.5)
# ax.set_xlabel('E-W distance from radar (km)')
# ax.set_ylabel('N-S distance from radar (km)')
# ax.set_title(f"{vol[vi].elev[eli].round(1)}\N{DEGREE SIGN} radial velocity - {filetime} UTC", fontsize=14)
# ax.text(6, 7, f"{azimuth}\N{DEGREE SIGN}", fontsize=16, fontweight='bold')
# # plt.suptitle(f"{filetime} UTC PPI at {vol[vi].elev[eli].round(1)}\N{DEGREE SIGN}, azimuth = {azimuth}\N{DEGREE SIGN}")
# if figsave:
#     plt.savefig(f"Documents/conferences/perilsmeeting/vol{vi}_{filetime}_az{azimuth}_PPI_V.png", dpi=400)

if False:
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)
    plot_cfill(vol[vi].xx[eli,:,:], vol[vi].yy[eli,:,:], vol[vi].vel[eli,:,:], 'vel', ax, datalims=[-va,va], xlims=[-rlim,rlim], ylims=[-rlim,rlim])
    # s1 = ax.scatter(P1.xx[i1], P1.yy[i1], s=30, c='k', marker='s')
    # s2 = ax.scatter(P2.xx[i2], P2.yy[i2], s=30, c='k', marker='^')
    # ax.text(P1.xx[i1]+0.3, P1.yy[i1], 'P1', fontsize=12, fontweight='bold')
    # ax.text(P2.xx[i2]+0.3, P2.yy[i2], 'P2', fontsize=12, fontweight='bold')
    # ax.scatter(0, 0, s=30, c='k')
    # ax.text(-1, 0.3, 'RaXPol', fontsize=12, fontweight='bold')
    # plt.legend([s1,s2], ['Probe 1', 'Probe 2'], loc='upper right')
    ax.set_xlabel('E-W distance from radar (km)', fontsize=12)
    ax.set_ylabel('N-S distance from radar (km)', fontsize=12)
    ax.set_title(f"{vol[vi].elev[eli].round(1)}\N{DEGREE SIGN} radial velocity - {filetime} UTC", fontsize=14, fontweight='bold')
    # plt.suptitle(f"{vol[vi].elev[eli].round(1)}\N{DEGREE SIGN} radial velocity - {filetime} UTC")
    if figsave:
        plt.savefig(f"Documents/conferences/perilsmeeting/vol{vi}_{filetime}_az{azimuth}_PPI_V_noAz.png", dpi=400)



#%% PERiLS 2023 Z, ZDR calibration constants

# IOP2
# KLZK beam height 1.2 km -- Rax range 3.7 km
fn_rax = 'Documents/perils2023/iop2/raxpol/CFradial/cfrad.20230303_080850_RaXPol_v38_s13.nc' # 17.5 deg
fn_88d = 'Documents/perils2023/iop2/klzk/KLZK20230303_080235_V06.ar2v'
# fn_rax = 'cfrad.20230303_081922_RaXPol_v97_s1.nc' # 19 deg
# fn_88d = 'Documents/perils2023/iop2/klzk/KLZK20230303_081934_V06.ar2v'

raxpol = pyart.io.read(fn_rax)
klzk = pyart.io.read(fn_88d)

rax = read_raxpol(fn_rax)
r_lzk = klzk.range['data']/1000
az_lzk = klzk.azimuth['data']
el_lzk = np.median(klzk.elevation['data'][0:719])

r_mat = np.tile(r_lzk, (len(az_lzk), 1))
az_mat = np.transpose(np.tile(az_lzk, (len(r_lzk), 1)))

xx_lzk = r_mat * np.sin(az_mat*np.pi/180) * np.cos(el_lzk*np.pi/180)
yy_lzk = r_mat * np.cos(az_mat*np.pi/180) * np.cos(el_lzk*np.pi/180)
zz_lzk = r_mat * np.sin(el_lzk*np.pi/180)



xr,yr = latlon2xy(raxpol.latitude['data'][0], raxpol.longitude['data'][0], klzk.latitude['data'][0], klzk.longitude['data'][0])


c1 = pyart.graph.RadarDisplay(raxpol)
c2 = pyart.graph.RadarDisplay(klzk)



#%%

fig = plt.figure(figsize=(14,5))

ax1 = fig.add_subplot(121)
c2.plot('reflectivity', 0, vmin=0, vmax=70, cmap='pyart_NWSRef')
c2.plot_range_rings([20,40,60,80,100])
ax2.set_xlim([-100,100])
ax2.set_ylim([-100,100])
# ax1.scatter(0, 0, s=50, c='k')
# ax1.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')

ax2 = fig.add_subplot(122)
c2.plot('velocity', vmin=-34, vmax=34, cmap='pyart_Carbone42')
c2.plot_range_rings([20,40,60,80,100])
ax2.set_xlim([-100,100])
ax2.set_ylim([-100,100])
ax2.scatter(xr, yr, s=25, c='k')

plt.show()

#%%

fig = plt.figure(figsize=(14,5))

ax1 = fig.add_subplot(121)
c1.plot('DBZ', 0, vmin=0, vmax=70, cmap='pyart_NWSRef')
c1.plot_range_rings([5, 10])
ax1.set_xlim([-10,10])
ax1.set_ylim([-10,10])
# ax1.scatter(0, 0, s=50, c='k')
# ax1.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')

ax2 = fig.add_subplot(122)
c2.plot('reflectivity', vmin=0, vmax=70, cmap='pyart_NWSRef')
c2.plot_range_rings([20,40,60,80,100])
ax2.set_xlim([-100,100])
ax2.set_ylim([-100,100])
ax2.scatter(xr, yr, s=25, c='k')

plt.show()



fig = plt.figure(figsize=(14,5))

ax1 = fig.add_subplot(121)
c1.plot('DBZ', 0, vmin=0, vmax=70, cmap='pyart_NWSRef')
c1.plot_range_rings([5, 10])
ax1.set_xlim([-10,10])
ax1.set_ylim([-10,10])
# ax1.scatter(0, 0, s=50, c='k')
# ax1.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')

ax2 = fig.add_subplot(122)
c2.plot('reflectivity', vmin=0, vmax=70, cmap='pyart_NWSRef')
c2.plot_range_rings([20,40,60,80,100])
ax2.set_xlim([xr-10,xr+10])
ax2.set_ylim([yr-10,yr+10])
ax2.scatter(xr, yr, s=25, c='k')

plt.show()



fig = plt.figure(figsize=(14,5))

ax1 = fig.add_subplot(121)
c1.plot('ZDR', 0, vmin=-5, vmax=5, cmap='pyart_NWSRef')
c1.plot_range_rings([5, 10])
ax1.set_xlim([-10,10])
ax1.set_ylim([-10,10])
# ax1.scatter(0, 0, s=50, c='k')
# ax1.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')

ax2 = fig.add_subplot(122)
c2.plot('differential_reflectivity', vmin=-5, vmax=5, cmap='pyart_NWSRef')
c2.plot_range_rings([20,40,60,80,100])
ax2.set_xlim([-100,100])
ax2.set_ylim([-100,100])
ax2.scatter(xr, yr, s=25, c='k')

plt.show()



fig = plt.figure(figsize=(14,5))

ax1 = fig.add_subplot(121)
c1.plot('ZDR', 0, vmin=-5, vmax=5, cmap='pyart_NWSRef')
c1.plot_range_rings([5, 10])
ax1.set_xlim([-10,10])
ax1.set_ylim([-10,10])
# ax1.scatter(0, 0, s=50, c='k')
# ax1.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')

ax2 = fig.add_subplot(122)
c2.plot('differential_reflectivity', vmin=-5, vmax=5, cmap='pyart_NWSRef')
c2.plot_range_rings([20,40,60,80,100])
ax2.set_xlim([xr-10,xr+10])
ax2.set_ylim([yr-10,yr+10])
# ax2.set_xlim([-100,100])
# ax2.set_ylim([-100,100])
ax2.scatter(xr, yr, s=25, c='k')

plt.show()



#%%

r_rax = (xr**2 + yr**2)**0.5
az_rax = np.arctan2(xr,yr)*180/np.pi
ir = np.where(r_lzk >= r_rax-4)[0][0]
# bh_lzk = np.median(klzk.gate_z['data'][0:359,ir])

# j1 = np.where(r_lzk >= r_rax-10)[0][-1]
# j2 = np.where(r_lzk <= r_rax+10)[0][0]

# i1 = np.where(az_lzk[0:719] >= 120)[0]
# i2 = np.where(az_lzk[0:719] <= 140)[0]

ii = slice(359,399)
jj = slice(318,358)

jj1 = slice(66,166)

# Z offset
zmean_lzk = np.min(klzk.fields['reflectivity']['data'][ii,jj])
zmean_rax = np.min(raxpol.fields['DBZ']['data'][:,jj1])

z_offset = zmean_lzk - zmean_rax

# ZDR offset
zdrmean_lzk = np.min(klzk.fields['differential_reflectivity']['data'][ii,jj])
zdrmean_rax = np.min(raxpol.fields['ZDR']['data'][:,jj1])

zdr_offset = zdrmean_lzk - zdrmean_rax




#%%

raxpol2 = raxpol
raxpol2.fields['DBZ']['data'] = raxpol2.fields['DBZ']['data'] + z_offset
raxpol2.fields['ZDR']['data'] = raxpol2.fields['ZDR']['data'] + zdr_offset

c3 = pyart.graph.RadarDisplay(raxpol2)



# fig = plt.figure(figsize=(14,5))

# ax1 = fig.add_subplot(121)
# c1.plot('DBZ', 0, vmin=0, vmax=70, cmap='pyart_NWSRef')
# c1.plot_range_rings([5, 10])
# ax1.set_xlim([-10,10])
# ax1.set_ylim([-10,10])
# # ax1.scatter(0, 0, s=50, c='k')
# # ax1.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')

# ax2 = fig.add_subplot(122)
# c2.plot('reflectivity', vmin=0, vmax=70, cmap='pyart_NWSRef')
# c2.plot_range_rings([20,40,60,80,100])
# ax2.set_xlim([-100,100])
# ax2.set_ylim([-100,100])
# ax2.scatter(xr, yr, s=25, c='k')

# plt.show()



fig = plt.figure(figsize=(14,5))

ax1 = fig.add_subplot(121)
c3.plot('DBZ', vmin=0, vmax=70, cmap='pyart_NWSRef')
c3.plot_range_rings([5, 10])
ax1.set_xlim([-10,10])
ax1.set_ylim([-10,10])
# ax1.scatter(0, 0, s=50, c='k')
# ax1.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')

ax2 = fig.add_subplot(122)
c2.plot('reflectivity', vmin=0, vmax=70, cmap='pyart_NWSRef')
c2.plot_range_rings([20,40,60,80,100])
ax2.set_xlim([xr-10,xr+10])
ax2.set_ylim([yr-10,yr+10])
ax2.scatter(xr, yr, s=25, c='k')

plt.show()



# fig = plt.figure(figsize=(14,5))

# ax1 = fig.add_subplot(121)
# c1.plot('ZDR', 0, vmin=-5, vmax=5, cmap='pyart_NWSRef')
# c1.plot_range_rings([5, 10])
# ax1.set_xlim([-10,10])
# ax1.set_ylim([-10,10])
# # ax1.scatter(0, 0, s=50, c='k')
# # ax1.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')

# ax2 = fig.add_subplot(122)
# c2.plot('differential_reflectivity', vmin=-5, vmax=5, cmap='pyart_NWSRef')
# c2.plot_range_rings([20,40,60,80,100])
# ax2.set_xlim([-100,100])
# ax2.set_ylim([-100,100])
# ax2.scatter(xr, yr, s=25, c='k')

# plt.show()



fig = plt.figure(figsize=(14,5))

ax1 = fig.add_subplot(121)
c3.plot('ZDR', vmin=-5, vmax=5, cmap='pyart_NWSRef')
c3.plot_range_rings([5, 10])
ax1.set_xlim([-10,10])
ax1.set_ylim([-10,10])
# ax1.scatter(0, 0, s=50, c='k')
# ax1.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')

ax2 = fig.add_subplot(122)
c2.plot('differential_reflectivity', vmin=-5, vmax=5, cmap='pyart_NWSRef')
c2.plot_range_rings([20,40,60,80,100])
ax2.set_xlim([xr-10,xr+10])
ax2.set_ylim([yr-10,yr+10])
# ax2.set_xlim([-100,100])
# ax2.set_ylim([-100,100])
ax2.scatter(xr, yr, s=25, c='k')

plt.show()


#%%
# IOP5
# KNQA beam height 1.2-1.5 km -- Rax range 2.5-3.1 km (29 deg), 3.1-4 km (23 deg), 4.6-5.8 km (15 deg)
fp = 'Documents/perils2023/iop5/raxpol/CFradial/'
fn_rax = fp+'cfrad.20230405_163236_RaXPol_v136_s3.nc'
fn_88d = 'Documents/perils2023/iop5/knqa/KNQA20230405_163216_V06.ar2v'
# knqa 162222, rax 162212 v123 s3
# knqa 162702, rax 162700 v129 s3
# knqa 163216, rax 163236 v136 s3/163148 v135 s3
raxpol = pyart.io.read(fn_rax)
knqa = pyart.io.read(fn_88d)
c1 = pyart.graph.RadarDisplay(raxpol)
c2 = pyart.graph.RadarDisplay(knqa)


xr,yr = latlon2xy(raxpol.latitude['data'][0], raxpol.longitude['data'][0], knqa.latitude['data'][0], knqa.longitude['data'][0])


fig = plt.figure(figsize=(14,5))

ax1 = fig.add_subplot(121)
c2.plot('reflectivity', sweep=0, vmin=0, vmax=70, cmap='pyart_NWSRef')
c2.plot_range_rings([25,50,75,100])
# ax1.set_xlim([-150,150])
# ax1.set_ylim([-150,150])
ax1.set_xlim([xr-15,xr+15])
ax1.set_ylim([yr-15,yr+15])
ax1.scatter(xr, yr, s=25, c='k')
ax1.scatter(xr-15, yr, s=10, c='k')
ax1.scatter(xr, yr+15, s=10, c='k')
# ax1.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')

ax2 = fig.add_subplot(122)
c2.plot('differential_reflectivity', sweep=0, vmin=-6, vmax=6, cmap='pyart_NWSRef')
c2.plot_range_rings([25,50,75,100])
# ax2.set_xlim([-150,150])
# ax2.set_ylim([-150,150])
ax2.set_xlim([xr-15,xr+15])
ax2.set_ylim([yr-15,yr+15])
ax2.scatter(xr, yr, s=25, c='k')

plt.show()



fig = plt.figure(figsize=(14,5))

ax1 = fig.add_subplot(121)
c1.plot('DBZ', vmin=0, vmax=70, cmap='pyart_NWSRef')
c1.plot_range_rings([5,10,15])
ax1.set_xlim([-15,15])
ax1.set_ylim([-15,15])
# ax1.scatter(xr, yr, s=25, c='k')
# ax1.scatter(xr-10, yr, s=10, c='k')
# ax1.scatter(xr, yr+10, s=10, c='k')
# ax1.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')

ax2 = fig.add_subplot(122)
c1.plot('ZDR', vmin=-6, vmax=6, cmap='pyart_NWSRef')
c1.plot_range_rings([5,10,15])
ax2.set_xlim([-15,15])
ax2.set_ylim([-15,15])
# ax2.scatter(xr, yr, s=25, c='k')

plt.show()

#%%

r_nqa = knqa.range['data']/1000
az_nqa = knqa.azimuth['data']
el_nqa = np.median(knqa.elevation['data'][0:719])

r_mat = np.tile(r_nqa, (len(az_nqa), 1))
az_mat = np.transpose(np.tile(az_nqa, (len(r_nqa), 1)))

xx_nqa = r_mat * np.sin(az_mat*np.pi/180) * np.cos(el_nqa*np.pi/180)
yy_nqa = r_mat * np.cos(az_mat*np.pi/180) * np.cos(el_nqa*np.pi/180)
zz_nqa = r_mat * np.sin(el_nqa*np.pi/180)



r_rax = (xr**2 + yr**2)**0.5
az_rax = np.arctan2(xr,yr)*180/np.pi + 360
ir = np.where(r_nqa >= r_rax-4)[0][0]
# bh_lzk = np.median(klzk.gate_z['data'][0:359,ir])

# j1 = np.where(r_lzk >= r_rax-10)[0][-1]
# j2 = np.where(r_lzk <= r_rax+10)[0][0]

# i1 = np.where(az_lzk[0:719] >= 120)[0]
# i2 = np.where(az_lzk[0:719] <= 140)[0]

ii = slice(52,88)
jj = slice(370,430)

ii1 = slice(88,268)
ii2 = slice(448,628)
jj1 = slice(333,400)

# Z offset
zmean_nqa = np.mean(knqa.fields['reflectivity']['data'][ii,jj])
zmean_rax = np.mean(raxpol.fields['DBZ']['data'][ii1,jj1])

z_offset = zmean_nqa - zmean_rax

# ZDR offset
zdrmean_nqa = np.mean(knqa.fields['differential_reflectivity']['data'][ii,jj])
zdrmean_rax = np.mean(raxpol.fields['ZDR']['data'][ii1,jj1])

zdr_offset = zdrmean_nqa - zdrmean_rax

#%%

raxpol2 = pyart.io.read(fn_rax)
raxpol2.fields['DBZ']['data'] = raxpol2.fields['DBZ']['data'] + z_offset
raxpol2.fields['ZDR']['data'] = raxpol2.fields['ZDR']['data'] + zdr_offset
c3 = pyart.graph.RadarDisplay(raxpol2)


fig = plt.figure(figsize=(14,5))

ax1 = fig.add_subplot(121)
c2.plot('reflectivity', sweep=0, vmin=0, vmax=70, cmap='pyart_NWSRef')
c2.plot_range_rings([25,50,75,100])
# ax1.set_xlim([-150,150])
# ax1.set_ylim([-150,150])
ax1.set_xlim([xr-15,xr+15])
ax1.set_ylim([yr-15,yr+15])
ax1.scatter(xr, yr, s=25, c='k')
ax1.scatter(xr-15, yr, s=10, c='k')
ax1.scatter(xr, yr+15, s=10, c='k')
# ax1.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')

ax2 = fig.add_subplot(122)
c2.plot('differential_reflectivity', sweep=0, vmin=-6, vmax=6, cmap='pyart_NWSRef')
c2.plot_range_rings([25,50,75,100])
# ax2.set_xlim([-150,150])
# ax2.set_ylim([-150,150])
ax2.set_xlim([xr-15,xr+15])
ax2.set_ylim([yr-15,yr+15])
ax2.scatter(xr, yr, s=25, c='k')

plt.show()



fig = plt.figure(figsize=(14,5))

ax1 = fig.add_subplot(121)
c3.plot('DBZ', vmin=0, vmax=70, cmap='pyart_NWSRef')
c3.plot_range_rings([5,10,15])
ax1.set_xlim([-15,15])
ax1.set_ylim([-15,15])
# ax1.scatter(xr, yr, s=25, c='k')
# ax1.scatter(xr-10, yr, s=10, c='k')
# ax1.scatter(xr, yr+10, s=10, c='k')
# ax1.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')

ax2 = fig.add_subplot(122)
c3.plot('ZDR', vmin=-6, vmax=6, cmap='pyart_NWSRef')
c3.plot_range_rings([5,10,15])
ax2.set_xlim([-15,15])
ax2.set_ylim([-15,15])
# ax2.scatter(xr, yr, s=25, c='k')

plt.show()







#%% PERiLS 2023 IOP5 heading correction

# if False:
#     fp = 'Documents/perils2023/iop5/raxpol/CFradial-corrected/'
#     files = sorted(glob(fp+'*.nc'))
    
#     for file in files:
#         ds = nc.Dataset(file, 'r+')
        
#         az = ds.variables['azimuth'][:] + 90
#         az_new = np.where(az > 360, az - 360, az)
#         ds.variables['azimuth'][:] = az_new
        
#         # az = ds.variables['azimuth'][:] - 90
#         # az_new = np.where(az < 0, az + 360, az)
#         # ds.variables['azimuth'][:] = az_new
        
#         ds.close()


#%% PERiLS 2023 IOP2 and IOP5 Z/ZDR calibration


z_offset = -6
zdr_offset = -3.5


fp = 'Documents/perils2023/iop2/raxpol/CFradial_cal/'
files = sorted(glob(fp+'*.nc'))

for file in files:
    ds = nc.Dataset(file, 'r+')
    
    z_new = ds.variables['DBZ'][:] + z_offset
    zdr_new = ds.variables['ZDR'][:] + zdr_offset
    ds.variables['DBZ'][:] = z_new
    ds.variables['ZDR'][:] = zdr_new
    
    ds.close()


fp = 'Documents/perils2023/iop5/raxpol/CFradial_cal/'
files = sorted(glob(fp+'*.nc'))

for file in files:
    ds = nc.Dataset(file, 'r+')
    
    z_new = ds.variables['DBZ'][:] + z_offset
    zdr_new = ds.variables['ZDR'][:] + zdr_offset
    ds.variables['DBZ'][:] = z_new
    ds.variables['ZDR'][:] = zdr_new
    
    ds.close()















