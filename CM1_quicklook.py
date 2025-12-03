#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 11:44:52 2023

@author: morgan.schneider

Quick look at Matt's supercell simulations for the general exam.
"""

from CM1utils import *
import os
from skimage.feature import peak_local_max, blob_doh, blob_dog, blob_log
from skimage.color import rgb2gray
from PIL import Image
from scipy.ndimage import gaussian_filter

#%% TLV times identified by criteria
# Matt's TG times:
# 53(R2), 79(R5), 62(R10), 83(R11), 61(R14), 76(R16), 55(R20), 55(R21), 55(R22), 53(R23),
# 49(R24), 53(R25), 52(R26), 56(R27), 56(R28), 55(R29), 52(R30), 53(R31), 53(R32), 52(R33)

# () indicates criteria met only at the surface

# R2: 53-83, (87), 119-120
# R5: 79-85, 87-91, 93, 110-112, 114
# R10: 62-64, 67-68, 70-87, 104-106, (109), (113), 117-120
# R11: 83-93, 98-101, 103-106, 115-117, 120
# R14: 61-71, 76-77, 87, (91), 103, (112)
# R16: 76-81, 83-91, (105), 117-120
# R20: 55-57, 59-67, (68), 74-76, 78-80, 82-86, 89, 91-95, 98, 101, 115-120
# R21: 55-73, 75-77, 93, (101), 107-108, 112-115, 117-118
# R22: 55-67, 69-74, 79, 81, (83), 112-113
# R23: 53-67, 69-76, 79-81, 97-100, 102-111, 119-120
# R24: 49-68, 79-80, 82-98, 101-102, 104-105, 107-109, 111-112, 115-116, 118-120
# R25: 53-72, 105-120
# R26: 52-80, 82-83, 87-92, 100-103, 105-111, 113-117, 120
# R27: 56-86, 106, 108-119
# R28: 56-57, 73-80, 82, 84
# R29: 55-64, 68-78, 92-93, 95-98, 109-115, 118-120
# R30: 52-77, 80-90, 98-101, 103, (107), 120
# R31: 53-73, (74), (84), 85, (87), 99, (100), 105, (107), 109-120
# R32: 55-56, 59, 63-81, 83-84, 86, 88-90, 92, 120
# R33: 52-62, 64-66, 68-72, 74-78, (81), 93-95, (102), 104, 106-107, 116-117, 119-120

### Near-field z_lfc = 1.835 km ###
### Far-field z_lfc = 2.115 km ###

#%% TLV criteria plot

r5 = 21 * np.ones(shape=(len(tlv_times.R5.tor),))
r10 = 20 * np.ones(shape=(len(tlv_times.R10.tor),))
r11 = 19 * np.ones(shape=(len(tlv_times.R11.tor),))
r14 = 18 * np.ones(shape=(len(tlv_times.R14.tor),))
r16 = 17 * np.ones(shape=(len(tlv_times.R16.tor),))
r2 = 15 * np.ones(shape=(len(tlv_times.R2.tor),))
r20 = 14 * np.ones(shape=(len(tlv_times.R20.tor),))
r21 = 13 * np.ones(shape=(len(tlv_times.R21.tor),))
r22 = 12 * np.ones(shape=(len(tlv_times.R22.tor),))
r23 = 11 * np.ones(shape=(len(tlv_times.R23.tor),))
r24 = 10 * np.ones(shape=(len(tlv_times.R24.tor),))
r25 = 9 * np.ones(shape=(len(tlv_times.R25.tor),))
r26 = 8 * np.ones(shape=(len(tlv_times.R26.tor),))
r27 = 7 * np.ones(shape=(len(tlv_times.R27.tor),))
r28 = 6 * np.ones(shape=(len(tlv_times.R28.tor),))
r29 = 5 * np.ones(shape=(len(tlv_times.R29.tor),))
r30 = 4 * np.ones(shape=(len(tlv_times.R30.tor),))
r31 = 3 * np.ones(shape=(len(tlv_times.R31.tor),))
r32 = 2 * np.ones(shape=(len(tlv_times.R32.tor),))
r33 = np.ones(shape=(len(tlv_times.R33.tor),))


fig = plt.figure(figsize=(16,10))

plt.grid(visible=True, which='both')
plt.scatter(tlv_times.R5.tor, r5, marker='o', c='k'); plt.plot(np.arange(40,121), 21*np.ones(shape=(81,)), 'k', linewidth=1)
plt.scatter(tlv_times.R10.tor, r10, marker='o', c='k'); plt.plot(np.arange(40,121), 20*np.ones(shape=(81,)), 'k', linewidth=1)
plt.scatter(tlv_times.R11.tor, r11, marker='o', c='k'); plt.plot(np.arange(40,121), 19*np.ones(shape=(81,)), 'k', linewidth=1)
plt.scatter(tlv_times.R14.tor, r14, marker='o', c='k'); plt.plot(np.arange(40,121), 18*np.ones(shape=(81,)), 'k', linewidth=1)
plt.scatter(tlv_times.R16.tor, r16, marker='o', c='k'); plt.plot(np.arange(40,121), 17*np.ones(shape=(81,)), 'k', linewidth=1)
plt.scatter(tlv_times.R2.tor, r2, marker='o', c='k'); plt.plot(np.arange(40,121), 15*np.ones(shape=(81,)), 'k', linewidth=1)
plt.scatter(tlv_times.R20.tor, r20, marker='o', c='k'); plt.plot(np.arange(40,121), 14*np.ones(shape=(81,)), 'k', linewidth=1)
plt.scatter(tlv_times.R21.tor, r21, marker='o', c='k'); plt.plot(np.arange(40,121), 13*np.ones(shape=(81,)), 'k', linewidth=1)
plt.scatter(tlv_times.R22.tor, r22, marker='o', c='k'); plt.plot(np.arange(40,121), 12*np.ones(shape=(81,)), 'k', linewidth=1)
plt.scatter(tlv_times.R23.tor, r23, marker='o', c='k'); plt.plot(np.arange(40,121), 11*np.ones(shape=(81,)), 'k', linewidth=1)
plt.scatter(tlv_times.R24.tor, r24, marker='o', c='k'); plt.plot(np.arange(40,121), 10*np.ones(shape=(81,)), 'k', linewidth=1)
plt.scatter(tlv_times.R25.tor, r25, marker='o', c='k'); plt.plot(np.arange(40,121), 9*np.ones(shape=(81,)), 'k', linewidth=1)
plt.scatter(tlv_times.R26.tor, r26, marker='o', c='k'); plt.plot(np.arange(40,121), 8*np.ones(shape=(81,)), 'k', linewidth=1)
plt.scatter(tlv_times.R27.tor, r27, marker='o', c='k'); plt.plot(np.arange(40,121), 7*np.ones(shape=(81,)), 'k', linewidth=1)
plt.scatter(tlv_times.R28.tor, r28, marker='o', c='k'); plt.plot(np.arange(40,121), 6*np.ones(shape=(81,)), 'k', linewidth=1)
plt.scatter(tlv_times.R29.tor, r29, marker='o', c='k'); plt.plot(np.arange(40,121), 5*np.ones(shape=(81,)), 'k', linewidth=1)
plt.scatter(tlv_times.R30.tor, r30, marker='o', c='k'); plt.plot(np.arange(40,121), 4*np.ones(shape=(81,)), 'k', linewidth=1)
plt.scatter(tlv_times.R31.tor, r31, marker='o', c='k'); plt.plot(np.arange(40,121), 3*np.ones(shape=(81,)), 'k', linewidth=1)
plt.scatter(tlv_times.R32.tor, r32, marker='o', c='k'); plt.plot(np.arange(40,121), 2*np.ones(shape=(81,)), 'k', linewidth=1)
plt.scatter(tlv_times.R33.tor, r33, marker='o', c='k'); plt.plot(np.arange(40,121), np.ones(shape=(81,)), 'k', linewidth=1)

ylab = ['R33','R32','R31','R30','R29','R28','R27','R26','R25','R24','R23','R22','R21','R20','R2','','R16','R14','R11','R10','R5']
plt.xlim([40,121])
plt.ylim([0,22])
plt.xticks(ticks=np.arange(40,121,5), fontsize=14)
plt.yticks(ticks=np.arange(1,22), labels=ylab, fontsize=14)
plt.xlabel('Time (min)', fontsize=16)
plt.ylabel('Simulation', fontsize=16)
plt.title('TLV criteria', fontsize=18)
plt.savefig('/Volumes/Promise_Pegasus_70TB/general/TLV_criteria.png', dpi=400)




#%% Intra-ensemble TLV statistics

# Total TLV time
# FF_time = [17, 30, 23, 15, 19] # mean=20.8, median=19.0, max=30, min=15
FF_time = [15, 30, 22, 13, 19] # mean=19.8, median=19.0, max=30, min=13
# NF_time = [33, 37, 31, 23, 42, 53, 36, 54, 44, 12, 37, 43, 36, 30, 34] # mean=36.3, median=36.0, max=54, min=12
NF_time = [33, 34, 30, 21, 42, 53, 36, 53, 43, 10, 37, 41, 33, 26, 33] # mean=35.0, median=34.0, max=53, min=10

# Max TLV duration
FF_duration = [7, 18, 11, 11, 9] # mean=11.2, median=11.0, max=18, min=7
NF_duration = [31, 9, 19, 13, 15, 20, 20, 29, 31, 8, 11, 26, 21, 19, 11] # mean=18.9, median=19.0, max=31, min=8




#%% Load data from one file

fp = '/Volumes/Promise_Pegasus_70TB/2021-SupSims-final/'
fd = 'R10_rerun'
sim_time = 85
# File number is 1 more than sim time

# Simulation output
ds = nc.Dataset(f"/Volumes/Promise_Pegasus_70TB/2021-SupSims-final/{fd}/cm1out_000001.nc")
prs0 = ds.variables['prs'][:].data # Pa
th0 = ds.variables['th'][:].data
qv0 = ds.variables['qv'][:].data
thr0 = th0 * (1 + 0.61*qv0) # K
# thv0 = thr0
# T0 = ds.variables['th'][:].data * (prs0/100000.)**0.286
# the0 = (T0 + (2.5e6*ds.variables['qv'][:].data)/(1004.5*T0)) * (100000./prs0)**0.286
u0 = ds.variables['uinterp'][:].data
v0 = ds.variables['vinterp'][:].data
ds.close()


ds = nc.Dataset(f"/Volumes/Promise_Pegasus_70TB/2021-SupSims-final/{fd}/cm1out_{sim_time+1:06d}.nc")

time = ds.variables['time'][:].data/60 # s -> min
xh = ds.variables['xh'][:].data # km
xf = ds.variables['xf'][:].data
yh = ds.variables['yh'][:].data
yf = ds.variables['yf'][:].data
z = ds.variables['z'][:].data
zf = ds.variables['zf'][:].data

prs = ds.variables['prs'][:].data # Pa
# pi = ds.variables['pi'][:].data
# pt1 = ds.variables['pt1'][:].data # Passive tracer kg/kg
th = ds.variables['th'][:].data
qv = ds.variables['qv'][:].data # kg/kg
dbz = ds.variables['dbz'][:].data # dBZ
qr = ds.variables['qr'][:].data
# qg = ds.variables['qg'][:].data
# qhl = ds.variables['qhl'][:].data
# crw = ds.variables['crw'][:].data
# chw = ds.variables['chw'][:].data
# chl = ds.variables['chl'][:].data
uinterp = ds.variables['uinterp'][:].data # m/s
vinterp = ds.variables['vinterp'][:].data
winterp = ds.variables['winterp'][:].data
# wspd = (uinterp**2 + vinterp**2)**0.5
# xvort = ds.variables['xvort'][:].data # 1/s
# yvort = ds.variables['yvort'][:].data # 1/s
zvort = ds.variables['zvort'][:].data # 1/s

thr = th * (1 + 0.61*qv - (ds.variables['qc'][:].data + ds.variables['qr'][:].data + 
                            ds.variables['qi'][:].data + ds.variables['qs'][:].data + 
                            ds.variables['qg'][:].data + ds.variables['qhl'][:].data)) # K
thrpert = thr - thr0 # K
prspert = prs - prs0 # Pa
thpert = th - th0
qvpert = qv - qv0

du_dx = np.gradient(uinterp, xh*1000, axis=3)
dv_dy = np.gradient(vinterp, yh*1000, axis=2)
divh = du_dx + dv_dy
B = 9.8 * (thpert/th0 + 0.61*qvpert - (ds.variables['qc'][:].data + ds.variables['qr'][:].data + 
                            ds.variables['qi'][:].data + ds.variables['qs'][:].data + 
                            ds.variables['qg'][:].data + ds.variables['qhl'][:].data))
ds.close()

F_tot,F_nl,F_l,F_b = pdcomp_forcing(xh, yh, z, uinterp, vinterp, winterp, u0, v0, B)




#%% Random plots of supercell and zoomed-in TLV

plt.close('all')
figsave = True

xlims = [-10,10]
ylims = [-10,10]
# xlims2 = [-2,8]
# ylims2 = [-4,6]
xlims2 = [1,5]
ylims2 = [0,4]
zlims = [0,2]

n = 'zoom_2'
ix = np.where(xh>=3.5)[0][0]
iy = np.where(yh>=1.8)[0][0]
iz = np.where(z<=0.75)[0][-1]

iz1 = np.where(z>=1)[0][0]

qix = 2
qiz = 3
qsc = 450

wlevs1 = [10,20]
wlevs2 = [-20,-10]
qrlevs = [1.0]
thlims = [-6,6]
fnlims = [-0.02,0.02]
fblims = [-0.002,0.002]
cm = cmocean.cm.balance

makefig = [1,1,1,1]

if False:
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    # c = plot_cfill(xf, yf, dbz[0,0,:,:], 'dbz', ax, datalims=[0,70], xlims=xlims, ylims=ylims)
    c = plot_cfill(xf, yf, zvort[0,0,:,:], 'zvort', ax, datalims=[-0.3,0.3], xlims=xlims2, ylims=ylims2)
    # ax.contour(xh, yh, dbz[0,0,:,:], levels=[15], colors='k', linewidths=1.5)
    ax.contour(xh, yh, thrpert[0,0,:,:], levels=[thr_lev], colors='blue', linewidths=1, linestyles='-')
    ax.contour(xh, yh, winterp[0,iz,:,:], levels=[10], colors='orange', linewidths=1.5)
    ax.quiver(xh[::4], yh[::4], uinterp[0,0,::4,::4], vinterp[0,0,::4,::4], scale=450, width=0.003, pivot='middle')
    ax.plot(xlims2, [yh[iy],yh[iy]], '--k', linewidth=1.5)
    ax.plot([xh[ix],xh[ix]], ylims2, '--k', linewidth=1.5)
    ax.set_xlabel('x distance (km)')
    ax.set_ylabel('y distance (km)')
    # ax.set_title(f"\u03B6$_{{10m}}$, $Z_{{10m}}$=15 dBZ, $w_{{1km}}$=10 m/s", fontsize=14)
    ax.set_title(f"\u03B6$_{{10m}}$, \u03B8\u1D68'$_{{10m}}$ = -2 K, $w_{{500m}}$ = 10 m/s", fontsize=14)
    
    plt.suptitle(f"$t={time[0]:.0f}$ min", fontsize=14)
    plt.show()

if False:
    iz3 = np.where(z>=3)[0][0]
    inds = [ np.unravel_index(np.argmax(zvort[0,i,:,:]), zvort[0,0,:,:].shape) for i in np.arange(0,iz3+1) ]
    iiy,iix = zip(*inds)
    iy1 = list(iiy)
    ix1 = list(iix)
    
    # fig,ax = plt.subplots(1,1,figsize=(8,6))
    # plot_cfill(xf, yf, zvort[0,33,:,:], 'zvort', ax, datalims=[-0.3,0.3])
    # ax.set_xlim(xlims2)
    # ax.set_ylim(ylims2)
    
    fig,ax = plt.subplots(1,1,figsize=(9,6))
    plot_cfill(xf, yf, thrpert[0,0,:,:], 'thrpert', ax, datalims=[-8,8], cmap=cmocean.cm.curl)
    ax.contour(xh, yh, dbz[0,0,:,:], levels=[10], colors='k', linewidths=1)
    # ax.contourf(xh, yh, winterp[0,iz,:,:], levels=[10,100], colors='orange', alpha=0.2)
    # ax.contourf(xh, yh, winterp[0,iz,:,:], levels=[-100,-10], colors='deepskyblue', alpha=0.2)
    ax.contour(xh, yh, winterp[0,iz,:,:], levels=[-10,10], colors='deepskyblue', linewidths=1)
    ax.quiver(xh[::qix], yh[::qix], uinterp[0,0,::qix,::qix], vinterp[0,0,::qix,::qix], scale=600, width=0.003, pivot='middle', color='gray')
    ax.plot(xh[ix1], yh[iy1], 'k', linewidth=1.5, zorder=10)
    c = ax.scatter(xh[ix1], yh[iy1], s=5, c=z[0:iz3+1], cmap=cmocean.cm.matter_r, zorder=20)
    ax.set_xlim(xlims2)
    ax.set_ylim(ylims2)
    cb = plt.colorbar(c)
    cb.set_label('Height (km)')
    ax.set_title("10 m \u03B8\u1D68', 1 km $w$, and trace of vorticity max")
    ax.set_xlabel('x distance (km)')
    ax.set_ylabel('y distance (km)')
    
    plt.suptitle(f"$t={time[0]:.0f}$ min", fontsize=14)
    if figsave:
        plt.savefig(f"/Volumes/Promise_Pegasus_70TB/general/{fd}/TLV_plots/vortex_trace_{time[0]:03.0f}min.png", dpi=400)


if makefig[0]:
    
    fig,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(12,9))
    
    plot_cfill(xf, yf, dbz[0,iz,:,:], 'dbz', ax1, datalims=[0,70], xlims=xlims2, ylims=ylims2)
    ax1.plot(xlims2, [yh[iy],yh[iy]], '-k', linewidth=1)
    ax1.plot([xh[ix],xh[ix]], ylims2, '-k', linewidth=1)
    ax1.contour(xh, yh, winterp[0,iz1,:,:], levels=wlevs1, colors='k', linewidths=1)
    ax1.contour(xh, yh, winterp[0,iz1,:,:], levels=wlevs2, colors='k', linewidths=1, linestyles='-.')
    ax1.contour(xh, yh, zvort[0,iz,:,:], levels=[0.2,0.3], colors='yellow', linewidths=1.2, linestyles='-')
    ax1.quiver(xh[::qix], yh[::qix], uinterp[0,iz,::qix,::qix], vinterp[0,iz,::qix,::qix], scale=500, width=0.003, pivot='middle')
    ax1.set_xlabel('x distance (km)')
    ax1.set_ylabel('y distance (km)')
    
    plot_cfill(xf, yf, thrpert[0,iz,:,:], 'thrpert', ax2, datalims=thlims, xlims=xlims2, ylims=ylims2, cmap=cmocean.cm.curl)
    ax2.plot(xlims2, [yh[iy],yh[iy]], '-k', linewidth=1)
    ax2.plot([xh[ix],xh[ix]], ylims2, '-k', linewidth=1)
    ax2.contour(xh, yh, winterp[0,iz1,:,:], levels=wlevs1, colors='k', linewidths=1, linestyles='-')
    ax2.contour(xh, yh, winterp[0,iz1,:,:], levels=wlevs2, colors='k', linewidths=1, linestyles='-.')
    # ax2.contour(xh, yh, prspert[0,iz,:,:]/100, levels=[-15,-10,-5,-2], colors='cyan', linewidths=1.2, linestyles='-.')
    # ax2.contour(xh, yh, prspert[0,iz,:,:]/100, levels=[1,5,10], colors='cyan', linewidths=1.2, linestyles='-')
    ax2.contour(xh, yh, zvort[0,iz,:,:], levels=[0.2,0.3], colors='yellow', linewidths=1.2, linestyles='-')
    ax2.contour(xh, yh, qr[0,iz,:,:]*1000, levels=qrlevs, colors='deepskyblue', linewidths=1)
    ax2.quiver(xh[::qix], yh[::qix], uinterp[0,iz,::qix,::qix], vinterp[0,iz,::qix,::qix], scale=500, width=0.003, pivot='middle')
    ax2.set_xlabel('x distance (km)')
    ax2.set_ylabel('y distance (km)')
    # ax2.set_title(f"\u03B6$_{{10m}}$, \u03B8\u1D68'$_{{10m}}$, $w_{{500m}}$", fontsize=14)
    
    plot_cfill(xf, yf, prspert[0,iz,:,:]/100, 'prspert', ax3, datalims=[-20,20], xlims=xlims2, ylims=ylims2)
    ax3.plot(xlims2, [yh[iy],yh[iy]], '-k', linewidth=1)
    ax3.plot([xh[ix],xh[ix]], ylims2, '-k', linewidth=1)
    ax3.contour(xh, yh, winterp[0,iz1,:,:], levels=wlevs1, colors='k', linewidths=1)
    ax3.contour(xh, yh, winterp[0,iz1,:,:], levels=wlevs2, colors='k', linewidths=1, linestyles='-.')
    ax3.contour(xh, yh, zvort[0,iz,:,:], levels=[0.2,0.3], colors='yellow', linewidths=1.2, linestyles='-')
    ax3.quiver(xh[::qix], yh[::qix], uinterp[0,iz,::qix,::qix], vinterp[0,iz,::qix,::qix], scale=500, width=0.003, pivot='middle')
    ax3.set_xlabel('x distance (km)')
    ax3.set_ylabel('y distance (km)')
    
    plot_cfill(xf, yf, divh[0,iz,:,:], 'divh', ax4, datalims=[-0.3,0.3], xlims=xlims2, ylims=ylims2)
    ax4.plot(xlims2, [yh[iy],yh[iy]], '-k', linewidth=1)
    ax4.plot([xh[ix],xh[ix]], ylims2, '-k', linewidth=1)
    ax4.contour(xh, yh, winterp[0,iz1,:,:], levels=wlevs1, colors='k', linewidths=1)
    ax4.contour(xh, yh, winterp[0,iz1,:,:], levels=wlevs2, colors='k', linewidths=1, linestyles='-.')
    ax4.contour(xh, yh, zvort[0,iz,:,:], levels=[0.2,0.3], colors='yellow', linewidths=1.2, linestyles='-')
    ax4.quiver(xh[::qix], yh[::qix], uinterp[0,iz,::qix,::qix], vinterp[0,iz,::qix,::qix], scale=500, width=0.003, pivot='middle')
    ax4.set_xlabel('x distance (km)')
    ax4.set_ylabel('y distance (km)')
    
    # plt.suptitle(f"\u03B6$_{{10m}}$, \u03B8\u1D68'$_{{10m}}$ = {thr_lev:.0f} K (blue), $w_{{500m}}$ = 10 m/s (orange) | $t={time[0]:.0f}$ min", fontsize=14)
    plt.suptitle(f"$z={z[iz]:.2f}$ km, $t={time[0]:.0f}$ min", fontsize=14)
    
    if figsave:
        # plt.savefig(f"/Volumes/Promise_Pegasus_70TB/general/{fd}/TLV_plots/sfc_{time[0]:03.0f}min_{n}.png", dpi=400)
        plt.savefig(f"/Volumes/Promise_Pegasus_70TB/general/{fd}/TLV_plots/birdseye_{z[iz]*1000:.0f}m_{time[0]:03.0f}min_{n}.png", dpi=400)
    
    
    
    
    
    
if makefig[1]:
    fig,((ax1,ax2),(ax3,ax4),(ax5,ax6)) = plt.subplots(3,2,figsize=(12,11))
    
    c1 = plot_cfill(xf, zf, thrpert[0,:,iy,:], 'thrpert', ax1, datalims=thlims, xlims=xlims2, ylims=zlims, cmap=cmocean.cm.curl)
    ax1.contour(xh, z, winterp[0,:,iy,:], levels=wlevs1, colors='k', linewidths=1, linestyles='-')
    ax1.contour(xh, z, winterp[0,:,iy,:], levels=wlevs2, colors='k', linewidths=1, linestyles='-.')
    ax1.contour(xh, z, qr[0,:,iy,:]*1000, levels=qrlevs, colors='deepskyblue', linewidths=1, linestyles='-.')
    # ax1.contour(xh, z, prspert[0,:,iy,:]/100, levels=[-15,-10,-5,-2], colors='aqua', linewidths=1.2, linestyles='-.')
    # ax1.contour(xh, z, prspert[0,:,iy,:]/100, levels=[1,5,10], colors='aqua', linewidths=1.2, linestyles='-')
    ax1.contour(xh, z, zvort[0,:,iy,:], levels=[0.1,0.2,0.3], colors='yellow', linewidths=1, linestyles='-')
    ax1.quiver(xh[::qix], z[::qiz], uinterp[0,::qiz,iy,::qix], winterp[0,::qiz,iy,::qix], scale=qsc, width=0.003, pivot='middle')
    ax1.set_xlabel('x distance (km)')
    ax1.set_ylabel('Height (km)')
    ax1.set_title(f"W-E cross-section, y={yh[iy]:.3f} km", fontsize=14)
    
    c2 = plot_cfill(yf, zf, thrpert[0,:,:,ix], 'thrpert', ax2, datalims=thlims, xlims=ylims2, ylims=zlims, cmap=cmocean.cm.curl)
    ax2.contour(yh, z, winterp[0,:,:,ix], levels=wlevs1, colors='k', linewidths=1, linestyles='-')
    ax2.contour(yh, z, winterp[0,:,:,ix], levels=wlevs2, colors='k', linewidths=1, linestyles='-.')
    ax2.contour(yh, z, qr[0,:,:,ix]*1000, levels=qrlevs, colors='deepskyblue', linewidths=1, linestyles='-.')
    # ax2.contour(yh, z, prspert[0,:,:,ix]/100, levels=[-15,-10,-5,-2], colors='aqua', linewidths=1.2, linestyles='-.')
    # ax2.contour(yh, z, prspert[0,:,:,ix]/100, levels=[1,5,10], colors='aqua', linewidths=1.2, linestyles='-')
    ax2.contour(yh, z, zvort[0,:,:,ix], levels=[0.1,0.2,0.3], colors='yellow', linewidths=1, linestyles='-')
    ax2.quiver(yh[::qix], z[::qiz], vinterp[0,::qiz,::qix,ix], winterp[0,::qiz,::qix,ix], scale=qsc, width=0.003, pivot='middle')
    ax2.set_xlabel('y distance (km)')
    ax2.set_ylabel('Height (km)')
    ax2.set_title(f"S-N cross-section, x={xh[ix]:.3f} km", fontsize=14)
    
    # plt.suptitle(f"\u03B6, \u03B8\u1D68' = {thr_lev:.0f} K (blue) | $t={time[0]:.0f}$ min", fontsize=14)
    # plt.suptitle(f"\u03B8\u1D68' (shaded), \u03B6 (yellow), $w$ (black) | $t={time[0]:.0f}$ min", fontsize=14)
    
    # plt.show()
    # if figsave:
    #     plt.savefig(f"/Volumes/Promise_Pegasus_70TB/general/{fd}/TLV_plots/vslice_thr+vort_{time[0]:03.0f}min_1.png", dpi=400)
    
    
    # fig,(ax1,ax2) = plt.subplots(1,2,figsize=(16,6))
    
    c3 = plot_cfill(xf, zf, prspert[0,:,iy,:]/100, 'prspert', ax3, datalims=[-20,20], xlims=xlims2, ylims=zlims, cmap=cmocean.cm.balance)
    ax3.contour(xh, z, winterp[0,:,iy,:], levels=wlevs1, colors='k', linewidths=1, linestyles='-')
    ax3.contour(xh, z, winterp[0,:,iy,:], levels=wlevs2, colors='k', linewidths=1, linestyles='-.')
    ax3.contour(xh, z, prspert[0,:,iy,:]/100, levels=[-10], colors='aqua', linewidths=1, linestyles='-.')
    # ax1.contour(xh, z, prspert[0,:,iy,:]/100, levels=[1,5,10], colors='aqua', linewidths=1.2, linestyles='-')
    ax3.contour(xh, z, zvort[0,:,iy,:], levels=[0.1,0.2,0.3], colors='yellow', linewidths=1, linestyles='-')
    ax3.quiver(xh[::qix], z[::qiz], uinterp[0,::qiz,iy,::qix], winterp[0,::qiz,iy,::qix], scale=qsc, width=0.003, pivot='middle')
    ax3.set_xlabel('x distance (km)')
    ax3.set_ylabel('Height (km)')
    # ax3.set_title(f"W-E cross-section, y={yh[iy]:.3f} km", fontsize=14)
    
    c4 = plot_cfill(yf, zf, prspert[0,:,:,ix]/100, 'prspert', ax4, datalims=[-20,20], xlims=ylims2, ylims=zlims, cmap=cmocean.cm.balance)
    ax4.contour(yh, z, winterp[0,:,:,ix], levels=wlevs1, colors='k', linewidths=1, linestyles='-')
    ax4.contour(yh, z, winterp[0,:,:,ix], levels=wlevs2, colors='k', linewidths=1, linestyles='-.')
    ax4.contour(yh, z, prspert[0,:,:,ix]/100, levels=[-10], colors='aqua', linewidths=1, linestyles='-.')
    # ax2.contour(yh, z, prspert[0,:,:,ix]/100, levels=[1,5,10], colors='aqua', linewidths=1.2, linestyles='-')
    ax4.contour(yh, z, zvort[0,:,:,ix], levels=[0.1,0.2,0.3], colors='yellow', linewidths=1, linestyles='-')
    ax4.quiver(yh[::qix], z[::qiz], vinterp[0,::qiz,::qix,ix], winterp[0,::qiz,::qix,ix], scale=qsc, width=0.003, pivot='middle')
    ax4.set_xlabel('y distance (km)')
    ax4.set_ylabel('Height (km)')
    # ax4.set_title(f"S-N cross-section, x={xh[ix]:.3f} km", fontsize=14)
    
    # plt.suptitle(f"\u03B6, \u03B8\u1D68' = {thr_lev:.0f} K (blue) | $t={time[0]:.0f}$ min", fontsize=14)
    # plt.suptitle(f"$p'$ (shaded), \u03B6 (yellow), $w$ (black) | $t={time[0]:.0f}$ min", fontsize=14)
    
    # plt.show()
    # if figsave:
    #     plt.savefig(f"/Volumes/Promise_Pegasus_70TB/general/{fd}/TLV_plots/vslice_pp+vort_{time[0]:03.0f}min_1.png", dpi=400)
    
    
    # fig,(ax1,ax2) = plt.subplots(1,2,figsize=(16,6))
    
    c5 = plot_cfill(xf, zf, winterp[0,:,iy,:], 'w', ax5, datalims=[-30,30], xlims=xlims2, ylims=zlims, cmap=cmocean.cm.balance)
    # ax1.contour(xh, z, prspert[0,:,iy,:]/100, levels=[-15,-10,-5,-2], colors='aqua', linewidths=1.2, linestyles='-.')
    # ax1.contour(xh, z, prspert[0,:,iy,:]/100, levels=[1,5,10], colors='aqua', linewidths=1.2, linestyles='-')
    ax5.contour(xh, z, zvort[0,:,iy,:], levels=[0.1,0.2,0.3], colors='yellow', linewidths=1, linestyles='-')
    ax5.quiver(xh[::qix], z[::qiz], uinterp[0,::qiz,iy,::qix], winterp[0,::qiz,iy,::qix], scale=qsc, width=0.003, pivot='middle')
    ax5.set_xlabel('x distance (km)')
    ax5.set_ylabel('Height (km)')
    # ax5.set_title(f"W-E cross-section, y={yh[iy]:.3f} km", fontsize=14)
    
    c6 = plot_cfill(yf, zf, winterp[0,:,:,ix], 'w', ax6, datalims=[-30,30], xlims=ylims2, ylims=zlims, cmap=cmocean.cm.balance)
    # ax2.contour(yh, z, prspert[0,:,:,ix]/100, levels=[-15,-10,-5,-2], colors='aqua', linewidths=1.2, linestyles='-.')
    # ax2.contour(yh, z, prspert[0,:,:,ix]/100, levels=[1,5,10], colors='aqua', linewidths=1.2, linestyles='-')
    ax6.contour(yh, z, zvort[0,:,:,ix], levels=[0.1,0.2,0.3], colors='yellow', linewidths=1, linestyles='-')
    ax6.quiver(yh[::qix], z[::qiz], vinterp[0,::qiz,::qix,ix], winterp[0,::qiz,::qix,ix], scale=qsc, width=0.003, pivot='middle')
    ax6.set_xlabel('y distance (km)')
    ax6.set_ylabel('Height (km)')
    # ax6.set_title(f"S-N cross-section, x={xh[ix]:.3f} km", fontsize=14)
    
    # plt.suptitle(f"\u03B6, \u03B8\u1D68' = {thr_lev:.0f} K (blue) | $t={time[0]:.0f}$ min", fontsize=14)
    plt.suptitle(f"Contours: \u03B6 (yellow), $w$ (black) | $t={time[0]:.0f}$ min", fontsize=14)
    
    plt.show()
    if figsave:
        plt.savefig(f"/Volumes/Promise_Pegasus_70TB/general/{fd}/TLV_plots/vslice_{time[0]:03.0f}min_{n}.png", dpi=400)
    


if makefig[2]:
    fig,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(12,9))
    # fig,(ax1,ax2) = plt.subplots(1,2,figsize=(14,5))
    
    c1 = plot_cfill(xf, zf, F_nl[0,:,iy,:], 'pipert', ax1, datalims=fnlims, 
                    xlims=xlims2, ylims=zlims, cmap=cm)
    ax1.contour(xh, z, winterp[0,:,iy,:], levels=wlevs1, colors='k', linewidths=1, linestyles='-')
    ax1.contour(xh, z, winterp[0,:,iy,:], levels=wlevs2, colors='k', linewidths=1, linestyles='-.')
    ax1.contour(xh, z, prspert[0,:,iy,:]/100, levels=[-10], colors='aqua', linewidths=1, linestyles='-.')
    ax1.contour(xh, z, zvort[0,:,iy,:], levels=[0.1,0.2,0.3], colors='yellow', linewidths=1, linestyles='-')
    ax1.quiver(xh[::qix], z[::qiz], uinterp[0,::qiz,iy,::qix], winterp[0,::qiz,iy,::qix], scale=qsc, width=0.003, pivot='middle')
    ax1.set_xlabel('x distance (km)')
    ax1.set_ylabel('Height (km)')
    ax1.set_title(f"W-E cross-section, y={yh[iy]:.3f} km", fontsize=14)
    
    c2 = plot_cfill(yf, zf, F_nl[0,:,:,ix], 'pipert', ax2, datalims=fnlims, 
                    xlims=ylims2, ylims=zlims, cmap=cm)
    ax2.contour(yh, z, winterp[0,:,:,ix], levels=wlevs1, colors='k', linewidths=1, linestyles='-')
    ax2.contour(yh, z, winterp[0,:,:,ix], levels=wlevs2, colors='k', linewidths=1, linestyles='-.')
    ax2.contour(yh, z, prspert[0,:,:,ix]/100, levels=[-10], colors='aqua', linewidths=1, linestyles='-.')
    ax2.contour(yh, z, zvort[0,:,:,ix], levels=[0.1,0.2,0.3], colors='yellow', linewidths=1, linestyles='-')
    ax2.quiver(yh[::qix], z[::qiz], vinterp[0,::qiz,::qix,ix], winterp[0,::qiz,::qix,ix], scale=qsc, width=0.003, pivot='middle')
    ax2.set_xlabel('y distance (km)')
    ax2.set_ylabel('Height (km)')
    ax2.set_title(f"S-N cross-section, x={xh[ix]:.3f} km", fontsize=14)
    
    
    c3 = plot_cfill(xf, zf, F_b[0,:,iy,:], 'pipert', ax3, datalims=fblims, 
                    xlims=xlims2, ylims=zlims, cmap=cm)
    ax3.contour(xh, z, winterp[0,:,iy,:], levels=wlevs1, colors='k', linewidths=1, linestyles='-')
    ax3.contour(xh, z, winterp[0,:,iy,:], levels=wlevs2, colors='k', linewidths=1, linestyles='-.')
    ax3.contour(xh, z, prspert[0,:,iy,:]/100, levels=[-10], colors='aqua', linewidths=1, linestyles='-.')
    ax3.contour(xh, z, zvort[0,:,iy,:], levels=[0.1,0.2,0.3], colors='yellow', linewidths=1, linestyles='-')
    # ax3.contour(xh, z, gaussian_filter(thrpert[0,:,iy,:],1), levels=[-5,-4,-3,-2,-1], colors='dodgerblue', linewidths=1, linestyles='-')
    # ax3.contour(xh, z, gaussian_filter(thrpert[0,:,iy,:],1), levels=[1,2,3,4,5], colors='orange', linewidths=1)
    # ax3.contour(xh, z, gaussian_filter(B[0,:,iy,:],1), levels=[-0.2,-0.15,-0.1,-0.05,-0.01], colors='dodgerblue', linewidths=1, linestyles='-')
    # ax3.contour(xh, z, gaussian_filter(B[0,:,iy,:],1), levels=[0.01,0.05,0.1,0.15,0.2], colors='orange', linewidths=1)
    # ax3.contour(xh, z, qr[0,:,iy,:]*1000, levels=qrlevs, colors='deepskyblue', linewidths=1)
    # ax3.contour(xh, z, crw[0,:,iy,:], levels=[100,1000,10000,100000], colors='r', linewidths=1)
    ax3.quiver(xh[::qix], z[::qiz], uinterp[0,::qiz,iy,::qix], winterp[0,::qiz,iy,::qix], scale=qsc, width=0.003, pivot='middle')
    ax3.set_xlabel('x distance (km)')
    ax3.set_ylabel('Height (km)')
    # ax3.set_title(f"W-E cross-section, y={yh[iy]:.3f} km", fontsize=14)
    
    c4 = plot_cfill(yf, zf, F_b[0,:,:,ix], 'pipert', ax4, datalims=fblims, 
                    xlims=ylims2, ylims=zlims, cmap=cm)
    ax4.contour(yh, z, winterp[0,:,:,ix], levels=wlevs1, colors='k', linewidths=1, linestyles='-')
    ax4.contour(yh, z, winterp[0,:,:,ix], levels=wlevs2, colors='k', linewidths=1, linestyles='-.')
    ax4.contour(yh, z, prspert[0,:,:,ix]/100, levels=[-10], colors='aqua', linewidths=1, linestyles='-.')
    ax4.contour(yh, z, zvort[0,:,:,ix], levels=[0.1,0.2,0.3], colors='yellow', linewidths=1, linestyles='-')
    # ax4.contour(yh, z, gaussian_filter(thrpert[0,:,:,ix],1), levels=[-5,-4,-3,-2,-1], colors='dodgerblue', linewidths=1, linestyles='-')
    # ax4.contour(yh, z, gaussian_filter(thrpert[0,:,:,ix],1), levels=[1,2,3,4,5], colors='orange', linewidths=1)
    # ax4.contour(yh, z, gaussian_filter(B[0,:,:,ix],1), levels=[-0.2,-0.15,-0.1,-0.05,-0.01], colors='dodgerblue', linewidths=1, linestyles='-')
    # ax4.contour(yh, z, gaussian_filter(B[0,:,:,ix],1), levels=[0.01,0.05,0.1,0.15,0.2], colors='orange', linewidths=1)
    # ax4.contour(yh, z, qr[0,:,:,ix]*1000, levels=qrlevs, colors='deepskyblue', linewidths=1)
    # ax4.contour(yh, z, crw[0,:,:,ix], levels=[100,1000,10000,100000], colors='r', linewidths=1)
    ax4.quiver(yh[::qix], z[::qiz], vinterp[0,::qiz,::qix,ix], winterp[0,::qiz,::qix,ix], scale=qsc, width=0.003, pivot='middle')
    ax4.set_xlabel('y distance (km)')
    ax4.set_ylabel('Height (km)')
    # ax4.set_title(f"S-N cross-section, x={xh[ix]:.3f} km", fontsize=14)
    
    # plt.suptitle(f"\u03B6, \u03B8\u1D68' = {thr_lev:.0f} K (blue) | $t={time[0]:.0f}$ min", fontsize=14)
    plt.suptitle(f"Nonlinear dynamic & buoyancy forcing, \u03B6 (yellow), $w$ (black) | $t={time[0]:.0f}$ min", fontsize=14)
    
    plt.show()
    if figsave:
        plt.savefig(f"/Volumes/Promise_Pegasus_70TB/general/{fd}/TLV_plots/ppnl_vslice_{time[0]:03.0f}min_{n}.png", dpi=400)
    
    
if makefig[3]:
    fig,(ax1,ax2) = plt.subplots(1,2,figsize=(14,5))
    
    plot_cfill(xf, yf, F_nl[0,iz,:,:], 'pipert', ax1, datalims=fnlims, 
               xlims=xlims2, ylims=ylims2, cmap=cm)
    ax1.contour(xh, yh, winterp[0,iz1,:,:], levels=wlevs2, colors='k', linewidths=1, linestyles='-.')
    ax1.contour(xh, yh, winterp[0,iz1,:,:], levels=wlevs1, colors='k', linewidths=1, linestyles='-')
    ax1.contour(xh, yh, prspert[0,iz,:,:]/100, levels=[-10], colors='aqua', linewidths=1, linestyles='-.')
    ax1.contour(xh, yh, zvort[0,iz,:,:], levels=[0.1,0.2,0.3], colors='yellow', linewidths=1)
    ax1.quiver(xh[::qix], yh[::qix], uinterp[0,iz,::qix,::qix], vinterp[0,iz,::qix,::qix], scale=500, width=0.003, pivot='middle')
    ax1.set_title('Nonlinear dynamic forcing')
    ax1.set_xlabel('x distance (km)')
    ax1.set_ylabel('y distance (km)')
    
    plot_cfill(xf, yf, F_b[0,iz,:,:], 'pipert', ax2, datalims=fblims, 
               xlims=xlims2, ylims=ylims2, cmap=cm)
    ax2.contour(xh, yh, winterp[0,iz1,:,:], levels=wlevs2, colors='k', linewidths=1, linestyles='-.')
    ax2.contour(xh, yh, winterp[0,iz1,:,:], levels=wlevs1, colors='k', linewidths=1, linestyles='-')
    ax2.contour(xh, yh, prspert[0,iz,:,:]/100, levels=[-10], colors='aqua', linewidths=1, linestyles='-.')
    ax2.contour(xh, yh, zvort[0,iz,:,:], levels=[0.1,0.2,0.3], colors='yellow', linewidths=1)
    # ax2.contourf(xh, yh, qr[0,iz,:,:]*1000, levels=[0.1,0.5,1.0,1.5,2.0], cmap='viridis', alpha=0.2)
    ax2.quiver(xh[::qix], yh[::qix], uinterp[0,iz,::qix,::qix], vinterp[0,iz,::qix,::qix], scale=500, width=0.003, pivot='middle')
    ax2.set_title('Buoyancy forcing')
    ax2.set_xlabel('x distance (km)')
    ax2.set_ylabel('y distance (km)')
    
    plt.suptitle(f"\u25BD$^2$p' forcing, $w$ (black) | $z={z[iz]:.2f}$ km, $t={time[0]:.0f}$ min", fontsize=14)
    
    plt.show()
    if figsave:
        # plt.savefig(f"/Volumes/Promise_Pegasus_70TB/general/{fd}/TLV_plots/ppsfc_{time[0]:03.0f}min_{n}.png", dpi=400)
        plt.savefig(f"/Volumes/Promise_Pegasus_70TB/general/{fd}/TLV_plots/pp{z[iz]*1000:.0f}m_{time[0]:03.0f}min_{n}.png", dpi=400)


#%% Tracking local maxima in vorticity and TLV flags (actually correct this time)

base_dir = '/Volumes/Promise_Pegasus_70TB/2021-SupSims-final/'
tor_members = [2,5,10,11,14,16,20,21,22,23,24,25,26,27,28,29,30,31,32,33]
sim_dirs = [ f"R{i}_rerun" for i in tor_members ]

tlv_flags = struct()
tlv_times = struct()
for fd in sim_dirs:
    print(f"--- Simulation {fd} ---")
    ds = nc.Dataset(base_dir+fd+'/cm1out_000001.nc')
    prs0 = ds.variables['prs'][:].data
    ds.close()
    
    files = glob(base_dir+fd+'/cm1out_00*.nc')
    tor_time = int(glob(base_dir+fd+'/pp/*.nc')[-1][-9:-3])
    files = files[tor_time::]
    
    flags = np.zeros(shape=(len(files),), dtype=float)
    sfc_flags = np.zeros(shape=(len(files),), dtype=float)
    vort_flags = np.zeros(shape=(len(files),), dtype=float)
    wspd_flags = np.zeros(shape=(len(files),), dtype=float)
    pp1_flags = np.zeros(shape=(len(files),), dtype=float)
    pps_flags = np.zeros(shape=(len(files),), dtype=float)
    locs = []; times = []; sfc_times = []; t = []
    for n in range(len(files)):
        print(f"File {files[n][-6:-3]}...")
        ds = nc.Dataset(files[n])
        time = ds.variables['time'][:].data[0]/60
        xh = ds.variables['xh'][:].data
        yh = ds.variables['yh'][:].data
        z = ds.variables['z'][:].data
        iz1 = np.where(z>=1)[0][0]
        zvort_sfc = ds.variables['zvort'][:].data[0,0,:,:]
        prspert = ds.variables['prs'][:].data[0,0:iz1,:,:] - prs0[0,0:iz1,:,:]
        wspd_sfc = (ds.variables['uinterp'][:].data[0,0,:,:]**2 + ds.variables['vinterp'][:].data[0,0,:,:]**2)**0.5
        ds.close()
        
        zvort_pos = np.where(zvort_sfc>0, zvort_sfc, 0)
        im = plt.imshow(zvort_pos, cmap='gray', origin='lower', vmin=0, vmax=0.3)
        img = Image.fromarray(np.uint8(im.get_cmap()(im.get_array())*255)).convert('RGB')
        zvort_img = rgb2gray(img)
        
        vort_blobs = blob_dog(zvort_img, min_sigma=0.1, max_sigma=30, threshold=0.05)
        vort_blobs[:,2] = vort_blobs[:,2] * np.sqrt(2)
        # xy = peak_local_max(zvort_img, min_distance=10, threshold_abs=0.1)
        
        if len(vort_blobs) > 0:
            flg = np.zeros(shape=(len(vort_blobs),), dtype=float)
            sfc_flg = np.zeros(shape=(len(vort_blobs),), dtype=float)
            zv_flg = np.zeros(shape=(len(vort_blobs),), dtype=float)
            ws_flg = np.zeros(shape=(len(vort_blobs),), dtype=float)
            pp1_flg = np.zeros(shape=(len(vort_blobs),), dtype=float)
            pps_flg = np.zeros(shape=(len(vort_blobs),), dtype=float)
            ncrit = np.zeros(shape=(len(vort_blobs),), dtype=float)
            xc = np.zeros(shape=(len(vort_blobs),), dtype=float)
            yc = np.zeros(shape=(len(vort_blobs),), dtype=float)
            for i in range(len(vort_blobs)):
                yi,xi,r = vort_blobs[i]
                idx = slice(int(xi-16), int(xi+17))
                jdy = slice(int(yi-16), int(yi+17))
                xc[i] = xh[int(xi)]
                yc[i] = yh[int(yi)]
                
                zv_flg[i] = 1 if np.any(zvort_sfc[jdy,idx]>=0.3, axis=(0,1)) else 0
                ws_flg[i] = 1 if np.any(wspd_sfc[jdy,idx]>=35, axis=(0,1)) else 0
                pp1_flg[i] = 1 if np.all(np.any(prspert[:,jdy,idx]/100<=-10, axis=(1,2)),axis=0) else 0
                pps_flg[i] = 1 if np.any(prspert[0,jdy,idx]/100<=-10, axis=(0,1)) else 0
                flg[i] = 1 if np.all([zv_flg[i], ws_flg[i], pp1_flg[i]]) else 0
                sfc_flg[i] = 1 if np.all([zv_flg[i], ws_flg[i], pps_flg[i]]) else 0
                ncrit[i] = np.count_nonzero([zv_flg[i], ws_flg[i], pps_flg[i]])
            
            flags[n] = 1 if np.any(flg) else 0
            sfc_flags[n] = 1 if np.any(sfc_flg) else 0
            t.append(int(time))
            if np.any(flg):
                i = np.nonzero(flg)[0][0]
                times.append(int(time))
                sfc_times.append(int(time))
                vort_flags[n] = 1
                wspd_flags[n] = 1
                pp1_flags[n] = 1
                pps_flags[n] = 1
                locs.append([xc[i],yc[i]])
            elif np.any(sfc_flg):
                i = np.nonzero(sfc_flg)[0][0]
                sfc_times.append(int(time))
                vort_flags[n] = 1
                wspd_flags[n] = 1
                pp1_flags[n] = pp1_flg[i]
                pps_flags[n] = 1
            else:
                i = np.where(ncrit == np.max(ncrit))[0][0]
                vort_flags[n] = zv_flg[i]
                wspd_flags[n] = ws_flg[i]
                pp1_flags[n] = pp1_flg[i]
                pps_flags[n] = pps_flg[i]
        else:
            flags[n] = 0; sfc_flags[n] = 0
            vort_flags[n] = 0; wspd_flags[n] = 0; pp1_flags[n] = 0; pps_flags[n] = 0
            t.append(int(time))
    
    
    tmp = struct()
    tmp.SetAttr('tor', flags)
    tmp.SetAttr('sfc', sfc_flags)
    tmp.SetAttr('vort', vort_flags)
    tmp.SetAttr('wspd', wspd_flags)
    tmp.SetAttr('pp_1km', pp1_flags)
    tmp.SetAttr('pp_sfc', pps_flags)
    tlv_flags.SetAttr(fd[0:-6], tmp)
    
    tmp2 = struct()
    tmp2.SetAttr('tor', np.array(times))
    tmp2.SetAttr('sfc', np.array(sfc_times))
    tmp2.SetAttr('all', np.array(t))
    tmp2.SetAttr('locs', np.array(locs))
    tlv_times.SetAttr(fd[0:-6], tmp2)
    
    del tmp,tmp2


#%% Make overview surface plots of all files in all simulations

base_dir = '/Volumes/Promise_Pegasus_70TB/2021-SupSims-final/'
#sim_dirs = os.listdir(base_dir)[1:-1]
# sim_dirs = [f"R{n}_rerun" for n in [2,5,10,11,14,16,21,24,25,26,27,29,32]]
# tor_members = [2,5,10,11,14,16,20,21,22,23,24,25,26,27,28,29,30,31,32,33]
# sim_dirs = [f"R{n}_rerun" for n in tor_members]
sim_dirs = ['R3_rerun']
figsave = False

for fd in sim_dirs:
    print(f"--- Simulation {fd} ---")
    # Base state
    fn_bs = base_dir+fd+'/cm1out_000001.nc'
    
    ds = nc.Dataset(fn_bs)
    th0 = ds.variables['th'][:].data
    qv0 = ds.variables['qv'][:].data
    prs0 = ds.variables['prs'][:].data
    # pi0 = ds.variables['pi'][:].data
    thr0 = th0 * (1 + 0.61*qv0)
    ds.close()
    
    # Simulation output
    files = glob(base_dir+fd+'/cm1out_00*.nc')
    # ppfiles = glob(base_dir+fd+'/pp/*.nc')
    # fnum = int(ppfiles[0][-9:-3])
    # files = files[fnum-1::]
    zvort_max = []
    wspd_max = []
    wmax = []
    ppmin = []
    times = []
    zvort_sfc = []
    wspd_sfc = []
    wmax_1km = []
    pp_1km = []
    for fn in files:
        print(f"--> File {fn[-6:-3]}...")
        
        ds = nc.Dataset(fn)
        
        time = ds.variables['time'][:].data # s
        xh = ds.variables['xh'][:].data # km
        xf = ds.variables['xf'][:].data
        yh = ds.variables['yh'][:].data
        yf = ds.variables['yf'][:].data
        z = ds.variables['z'][:].data
        zf = ds.variables['zf'][:].data
        
        th = ds.variables['th'][:].data # K
        prs = ds.variables['prs'][:].data # Pa
        # pi = ds.variables['pi'][:].data
        # pt1 = ds.variables['pt1'][:].data # kg/kg - passive tracer mixing ratio
        qv = ds.variables['qv'][:].data # kg/kg
        qc = ds.variables['qc'][:].data
        qr = ds.variables['qr'][:].data
        qi = ds.variables['qi'][:].data
        qs = ds.variables['qs'][:].data
        qg = ds.variables['qg'][:].data
        qhl = ds.variables['qhl'][:].data
        dbz = ds.variables['dbz'][:].data # dBZ
        
        uinterp = ds.variables['uinterp'][:].data # m/s
        vinterp = ds.variables['vinterp'][:].data
        wspd = ((uinterp**2) + (vinterp**2))**0.5
        winterp = ds.variables['winterp'][:].data
        # xvort = ds.variables['xvort'][:].data # 1/s
        # yvort = ds.variables['yvort'][:].data # 1/s
        zvort = ds.variables['zvort'][:].data # 1/s
        
        qx = qc + qr + qi + qs + qg + qhl
        thr = th * (1 + 0.61*qv - qx)
        thrpert = thr - thr0
        prspert = prs - prs0
        # qvpert = qv - qv0
        # thpert = th - th0
        # pipert = pi - pi0
        # B = 9.8 * (thpert/th0 + 0.61*qvpert - qx)
        
        ds.close()
        
        iz = np.where(z>=1)[0][0]
        
        times.append(time[0])
        zvort_max.append(np.max(zvort[0,:,:,:], axis=(1,2)))
        wspd_max.append(np.max(wspd[0,:,:,:], axis=(1,2)))
        wmax.append(np.max(winterp[0,:,:,:], axis=(1,2)))
        ppmin.append(np.min(prspert[0,:,:,:], axis=(1,2)))
        
        zvort_sfc.append(np.max(zvort[0,0,:,:], axis=(0,1)))
        wspd_sfc.append(np.max(wspd[0,0,:,:], axis=(0,1)))
        wmax_1km.append(np.max(winterp[0,iz,:,:], axis=(0,1)))
        pp_1km.append(np.mean(np.min(prspert[0,::iz,:,:], axis=(1,2)), axis=0))
        
        del qc,qr,qi,qs,qg,qhl
        
        if False:
            fig,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(12,8))
            
            c1 = plot_cfill(xf, yf, dbz[0,0,:,:], 'dbz', ax1, datalims=[0,70])
            # q1 = ax1.quiver(xh[::20], yh[::20], uinterp[0,0,::20,::20], vinterp[0,0,::20,::20], scale=400, width=0.003, pivot='middle')
            ax1.set_ylabel('y distance (km)')
            ax1.set_title(f"$Z_H$ at {1000*z[0]:.0f} m AGL", fontsize=12)
            ax1.set_xlim([-30,30])
            ax1.set_ylim([-20,40])
            
            c2 = plot_cfill(xf, yf, prspert[0,0,:,:]/100, 'prspert', ax2, datalims=[-10,10])
            # q2 = ax2.quiver(xh[::25], yh[::25], uinterp[0,0,::25,::25], vinterp[0,0,::25,::25], scale=500, width=0.003, pivot='middle')
            ax2.set_title(f"$p$' at {1000*z[0]:.0f} m AGL", fontsize=12)
            ax2.set_xlim([-30,30])
            ax2.set_ylim([-20,40])
            
            c3 = plot_cfill(xf, yf, thrpert[0,0,:,:], 'thrpert', ax3, datalims=[-10,10])
            # q3 = ax3.quiver(xh[::20], yh[::20], uinterp[0,0,::20,::20], vinterp[0,0,::20,::20], scale=500, width=0.003, pivot='middle')
            ax3.set_xlabel('x distance (km)')
            ax3.set_ylabel('y distance (km)')
            ax3.set_title(f"\u03B8\u1D68' at {1000*z[0]:.0f} m AGL", fontsize=12)
            ax3.set_xlim([-30,30])
            ax3.set_ylim([-20,40])
            
            c4 = plot_cfill(xf, yf, zvort[0,0,:,:], 'zvort', ax4, datalims=[-0.3,0.3])
            # q4 = ax4.quiver(xh[::20], yh[::20], uinterp[0,0,::20,::20], vinterp[0,0,::20,::20], scale=500, width=0.003, pivot='middle')
            ax4.set_xlabel('x distance (km)')
            ax4.set_title(f"\u03B6 at {1000*z[0]:.0f} m AGL", fontsize=12)
            ax4.set_xlim([-30,30])
            ax4.set_ylim([-20,40])
            
            plt.suptitle(f"$t={time[0]/60:.0f}$ min", fontsize=14)
            
            if figsave:
                plt.savefig('/Volumes/Promise_Pegasus_70TB/general/'+fd+f"/full_sim_plots/sfc_{fn[-6:-3]}_{time[0]/60:.0f}min.png")
            # plt.show()
            
            plt.close()
        
        # fig,(ax1,ax2) = plt.subplots(1,2,figsize=(14,5))
        
        # c1 = plot_cfill(xf, yf, thrpert[0,0,:,:], 'thrpert', ax1, datalims=[-10,10])
        # # q1 = ax1.quiver(xh[::20], yh[::20], uinterp[0,0,::20,::20], vinterp[0,0,::20,::20], scale=500, width=0.003, pivot='middle')
        # ax1.contour(xh, yh, zvort[0,0,:,:], levels=[0.1, 0.2, 0.3, 0.4, 0.5], colors=['k','k','r','r','r'], linewidths=1.25)
        # ax1.contour(xh, yh, winterp[0,iz,:,:], levels=[5,15,25,35], colors='orange', linewidths=1.5)
        # ax1.set_xlabel('x distance (km)')
        # ax1.set_ylabel('y distance (km)')
        # ax1.set_xlim([-20,40])
        # ax1.set_ylim([-20,40])
        
        # c2 = plot_cfill(xf, yf, thrpert[0,0,:,:], 'thrpert', ax2, datalims[-10,10])
        # ax2.contour(xh, yh, prspert[0,0,:,:]/100, levels=[-10,-8,-6,-4,-2], colors=['r','k','k','k','k'], linewidths=1.25)
        # ax2.contour(xh, yh, prspert[0,0,:,:]/100, levels=[2,4,6,8,10], colors=['k','k','k','k','r'], linewidths=1.25)
        # ax2.set_xlabel('x distance (km)')
        # ax2.set_ylabel('y distance (km)')
        # ax2.set_xlim([-20,40])
        # ax2.set_ylim([-20,40])
    
    times = np.array(times)
    zvort_max = np.array(zvort_max)
    wspd_max = np.array(wspd_max)
    wmax = np.array(wmax)
    ppmin = np.array(ppmin)
    
    zvort_sfc = np.array(zvort_sfc)
    wspd_sfc = np.array(wspd_sfc)
    wmax_1km = np.array(wmax_1km)
    pp_1km = np.array(pp_1km)
    
    
    
    fig,((ax1),(ax2),(ax3)) = plt.subplots(3,1,figsize=(12,12))
    
    c1 = plot_cfill(times/60, z, zvort_max.transpose(), 'zvort', ax1, datalims=[0,0.5], cmap='Reds')
    ax1.contour(times/60, z, zvort_max.transpose(), levels=[0.3], colors='w', linewidths=1.5)
    # ax1.set_xlabel('Time (min)')
    ax1.set_ylabel('Height (km)')
    ax1.set_ylim([0,5])
    ax1.set_title('Maximum vertical vorticity (1/s)')
    
    c2 = plot_cfill(times/60, z, wmax.transpose(), 'w', ax2, datalims=[0,50], cmap='Reds')
    # ax2.contour(times/60, z, wmax.transpose(), levels=[20], colors='w', linewidths=1.5)
    # ax2.set_xlabel('Time (min)')
    ax2.set_ylabel('Height (km)')
    ax2.set_ylim([0,5])
    ax2.set_title('Maximum updraft speed (m/s)')
    
    c3 = plot_cfill(times/60, z, ppmin.transpose()/100, 'prspert', ax3, datalims=[-30,0], cmap='Blues_r')
    ax3.contour(times/60, z, ppmin.transpose()/100, levels=[-10], colors='k', linewidths=1.5, linestyles='-')
    ax3.set_xlabel('Time (min)')
    ax3.set_ylabel('Height (km)')
    ax3.set_ylim([0,5])
    ax3.set_title('Minimum pressure perturbation (hPa)')
    
    plt.suptitle(f"{fd}")
    plt.savefig('/Volumes/Promise_Pegasus_70TB/general/'+fd+'/timeheight.png')
    # plt.close()
    
    
    # zvort_flag = np.where(zvort_sfc>=0.3, 1, 0)
    # wspd_flag = np.where(wspd_sfc>=35, 1, 0)
    # pp_flag = np.all(np.where(ppmin[:,0:iz]/100<=-10, 1, 0), axis=1)
    
    
    # fig,ax = plt.subplots(1,1,figsize=(14,4))
    
    # ax.plot(times/60, 3*np.ones(shape=times.shape), 'b', linewidth=0.5)
    # ax.plot(times/60, 2*np.ones(shape=times.shape), 'r', linewidth=0.5)
    # ax.plot(times/60, np.ones(shape=times.shape), 'g', linewidth=0.5)
    # ax.scatter(times[zvort_flag==1]/60, 3*np.ones(shape=times[zvort_flag==1].shape), s=20, c='b', marker='.')
    # ax.scatter(times[zvort_flag==0]/60, 3*np.ones(shape=times[zvort_flag==0].shape), s=20, c='b', marker='x', linewidth=1)
    # ax.scatter(times[wspd_flag==1]/60, 2*np.ones(shape=times[wspd_flag==1].shape), s=20, c='r', marker='.')
    # ax.scatter(times[wspd_flag==0]/60, 2*np.ones(shape=times[wspd_flag==0].shape), s=20, c='r', marker='x', linewidth=1)
    # ax.scatter(times[pp_flag==1]/60, np.ones(shape=times[pp_flag==1].shape), s=20, c='g', marker='.')
    # ax.scatter(times[pp_flag==0]/60, np.ones(shape=times[pp_flag==0].shape), s=20, c='g', marker='x', linewidth=1)
    
    # ax.set_xlabel('Time (min)')
    # ax.set_xlim([-1,121])
    # ax.set_ylim([0,4])
    # ax.set_yticks([1,2,3])
    # ax.set_yticklabels(["0-1 km p' < -10 hPa", "10 m wspd > 35 m/s", "10 m \u03B6 > 0.3 s$^{-1}$"])
    
    



#%% Base state LFC

### Near-field z_lfc = 1.835 km ###
### Far-field z_lfc = 2.115 km ###


# sim_dirs = os.listdir('/Volumes/Promise_Pegasus_70TB/general/')[1:-1]

# Near field environmental LFC
ds = nc.Dataset('/Volumes/Promise_Pegasus_70TB/2021-SupSims-final/R2_rerun/cm1out_000001.nc')
z = ds.variables['z'][:].data
th0 = ds.variables['th'][:].data[0,:,430,430]
qv0 = ds.variables['qv'][:].data[0,:,430,430]
prs0 = ds.variables['prs'][:].data[0,:,430,430]
# u0 = ds.variables['uinterp'][:].data[0,:,430,430]
# v0 = ds.variables['vinterp'][:].data[0,:,430,430]
ds.close()

T0 = th0 * (prs0/100000.)**0.286
e0 = (qv0 * prs0/100) / (0.622+qv0)
Td0 = 243.5 / ((17.67/(np.log(e0/6.112)))-1) + 273.15

LFC = mc.lfc(prs0*units.Pa, T0*units.K, Td0*units.K)
pLFC = LFC[0].magnitude
iz1 = np.where(prs0>=pLFC)[0][-1]
iz2 = np.where(prs0<=pLFC)[0][0]
zLFC_near = z[iz2] - (z[iz2]-z[iz1]) * (pLFC-prs0[iz2]) / (prs0[iz1]-prs0[iz2])


# Far field environmental LFC
ds = nc.Dataset('/Volumes/Promise_Pegasus_70TB/2021-SupSims-final/R3_rerun/cm1out_000001.nc')
z = ds.variables['z'][:].data
th0 = ds.variables['th'][:].data[0,:,430,430]
qv0 = ds.variables['qv'][:].data[0,:,430,430]
prs0 = ds.variables['prs'][:].data[0,:,430,430]
# u0 = ds.variables['uinterp'][:].data[0,:,430,430]
# v0 = ds.variables['vinterp'][:].data[0,:,430,430]
ds.close()

T0 = th0 * (prs0/100000.)**0.286
e0 = (qv0 * prs0/100) / (0.622+qv0)
Td0 = 243.5 / ((17.67/(np.log(e0/6.112)))-1) + 273.15

LFC = mc.lfc(prs0*units.Pa, T0*units.K, Td0*units.K)
pLFC = LFC[0].magnitude
iz1 = np.where(prs0>=pLFC)[0][-1]
iz2 = np.where(prs0<=pLFC)[0][0]
zLFC_far = z[iz2] - (z[iz2]-z[iz1]) * (pLFC-prs0[iz2]) / (prs0[iz1]-prs0[iz2])

# Check LFC of environmental inflow, primary RFGF, RFDIS at time of dissipation?
# This will have to be done manually -- maybe only do for a couple of cases

#%% More LFC calculations

simtimes = [77,78,79,80,81,82,83,84,85,86,87,88]

# Base state
ds = nc.Dataset('/Volumes/Promise_Pegasus_70TB/2021-SupSims-final/R10_rerun/cm1out_000001.nc')
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
z = ds.variables['z'][:].data
th0 = ds.variables['th'][:].data
qv0 = ds.variables['qv'][:].data
prs0 = ds.variables['prs'][:].data
thr0 = th0 * (1 + 0.61*qv0)
ds.close()

zLFC = np.zeros(shape=(len(simtimes),), dtype=float)
zLCL = np.zeros(shape=(len(simtimes),), dtype=float)
for t in simtimes:
    ds = nc.Dataset(f"/Volumes/Promise_Pegasus_70TB/2021-SupSims-final/R10_rerun/cm1out_{t+1:06d}.nc")
    zvort = ds.variables['zvort'][:].data
    uinterp = ds.variables['uinterp'][:].data
    vinterp = ds.variables['vinterp'][:].data
    th = ds.variables['th'][:].data
    qv = ds.variables['qv'][:].data
    prs = ds.variables['prs'][:].data
    thr = th * (1 + 0.61*qv - (ds.variables['qc'][:].data + ds.variables['qr'][:].data + 
                               ds.variables['qi'][:].data + ds.variables['qs'][:].data + 
                               ds.variables['qg'][:].data + ds.variables['qhl'][:].data)) # K
    thrpert = thr - thr0 # K
    B = 9.8 * (thrpert / thr0)
    wspd = (uinterp**2 + vinterp**2)**0.5
    ds.close()
    
    inds = np.where(zvort[0,0,:,:] == np.max(zvort[0,0,:,:], axis=(0,1)))
    iy = inds[0][0]
    ix = inds[1][0]
    
    T = th[0,:,iy,ix] * (prs[0,:,iy,ix]/100000.)**0.286
    e = (qv[0,:,iy,ix] * prs[0,:,iy,ix]/100) / (0.622+qv[0,:,iy,ix])
    Td = 243.5 / ((17.67/(np.log(e/6.112)))-1) + 273.15
    
    LFC = mc.lfc(prs[0,:,iy,ix]*units.Pa, T*units.K, Td*units.K)
    pLFC = LFC[0].magnitude
    iz1 = np.where(prs[0,:,iy,ix]>=pLFC)[0][-1]
    iz2 = np.where(prs[0,:,iy,ix]<=pLFC)[0][0]
    zLFC[t-77] = z[iz2] - (z[iz2]-z[iz1]) * (pLFC-prs[0,iz2,iy,ix]) / (prs[0,iz1,iy,ix]-prs[0,iz2,iy,ix])
    
    LCL = mc.lcl(prs[0,0,iy,ix]*units.Pa, T[0]*units.K, Td[0]*units.K)
    pLCL = LCL[0].magnitude
    iz1 = np.where(prs[0,:,iy,ix]>=pLCL)[0][-1]
    iz2 = np.where(prs[0,:,iy,ix]<=pLCL)[0][0]
    zLCL[t-77] = z[iz2] - (z[iz2]-z[iz1]) * (pLCL-prs[0,iz2,iy,ix]) / (prs[0,iz1,iy,ix]-prs[0,iz2,iy,ix])
    








#%%

fp = '/Volumes/Promise_Pegasus_70TB/2021-SupSims-final/'
fd = 'R27_rerun'
sim_time = 87
# File number is 1 more than sim time

# Simulation output
ds = nc.Dataset(f"/Volumes/Promise_Pegasus_70TB/2021-SupSims-final/{fd}/cm1out_000001.nc")
u0 = ds.variables['uinterp'][:].data
v0 = ds.variables['vinterp'][:].data
ds.close()


ds = nc.Dataset(f"/Volumes/Promise_Pegasus_70TB/2021-SupSims-final/{fd}/cm1out_{sim_time+1:06d}.nc")
time = ds.variables['time'][:].data/60 # s -> min
xh = ds.variables['xh'][:].data # km
xf = ds.variables['xf'][:].data
yh = ds.variables['yh'][:].data
yf = ds.variables['yf'][:].data
z = ds.variables['z'][:].data
zf = ds.variables['zf'][:].data

qr = ds.variables['qr'][:].data # kg/kg
dbz = ds.variables['dbz'][:].data # dBZ
uinterp = ds.variables['uinterp'][:].data # m/s
vinterp = ds.variables['vinterp'][:].data
winterp = ds.variables['winterp'][:].data
zvort = ds.variables['zvort'][:].data # 1/s

du_dx = np.gradient(uinterp, xh*1000, axis=3)
dv_dy = np.gradient(vinterp, yh*1000, axis=2)
divh = du_dx + dv_dy
ds.close()

# pgfx = struct(); pgfy = struct(); pgfz = struct()

ds = nc.Dataset(f"/Volumes/Promise_Pegasus_70TB/2021-SupSims-final/{fd}/pp/dyn_{sim_time+1:06d}.nc")
thrpert = ds.variables['thrpert'][:].data
p_dl = ds.variables['p_dl'][:].data # linear dynamic prspert
p_dn = ds.variables['p_dn'][:].data # nonlinear dynamic prspert
p_b = ds.variables['p_b'][:].data # buoyancy prspert
prspert = p_dl + p_dn + p_b
# pgfx.SetAttr('dl', ds.variables['pgfx_dl'][:].data) # linear dynamic ppgf
# pgfy.SetAttr('dl', ds.variables['pgfy_dl'][:].data)
# pgfz.SetAttr('dl', ds.variables['pgfz_dl'][:].data)
# pgfx.SetAttr('dn', ds.variables['pgfx_dn'][:].data) # nonlinear dynamic ppgf
# pgfy.SetAttr('dn', ds.variables['pgfy_dn'][:].data)
# pgfz.SetAttr('dn', ds.variables['pgfz_dn'][:].data)
# pgfx.SetAttr('b', ds.variables['pgfx_b'][:].data) # buoyancy ppgf
# pgfy.SetAttr('b', ds.variables['pgfy_b'][:].data)
# pgfz.SetAttr('b', ds.variables['pgfz_b'][:].data)
# pi_dl = ds.variables['pi_dl'][:].data # linear dynamic pipert
# pi_dn = ds.variables['pi_dn'][:].data # nonlinear dynamic pipert
# pi_b = ds.variables['pi_b'][:].data # buoyancy pipert
ds.close()



#%% Random plots of supercell and zoomed-in TLV

# plt.close('all')
figsave = False

xlims = [-10,10]
ylims = [-10,10]
xlims2 = [-4,2]
ylims2 = [2,8]
# xlims2 = [1,5]
# ylims2 = [0,4]
zlims = [0,2]

n = 'zoom_2'
ix = np.where(xh>=-1.9)[0][0]
iy = np.where(yh>=5.3)[0][0]
# iz = np.where(z<=0.75)[0][-1]
iz = np.where(z>=0.0)[0][0]

iz1 = np.where(z>=1)[0][0]

qix = 2
qiz = 3
qsc = 450

wlevs1 = [10,20]
wlevs2 = [-20,-10]
qrlevs = [1.0]
thlims = [-6,6]
pplims = [-10,10]

cm = cmocean.cm.balance

makefig = [0,0,1,1]

if False:
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    # c = plot_cfill(xf, yf, dbz[0,0,:,:], 'dbz', ax, datalims=[0,70], xlims=xlims, ylims=ylims)
    c = plot_cfill(xf, yf, zvort[0,0,:,:], 'zvort', ax, datalims=[-0.3,0.3], xlims=xlims2, ylims=ylims2)
    ax.contour(xh, yh, thrpert[0,0,:,:], levels=[thr_lev], colors='blue', linewidths=1, linestyles='-')
    ax.contour(xh, yh, winterp[0,iz,:,:], levels=[10], colors='orange', linewidths=1.5)
    ax.quiver(xh[::4], yh[::4], uinterp[0,0,::4,::4], vinterp[0,0,::4,::4], scale=450, width=0.003, pivot='middle')
    ax.plot(xlims2, [yh[iy],yh[iy]], '--k', linewidth=1.5)
    ax.plot([xh[ix],xh[ix]], ylims2, '--k', linewidth=1.5)
    ax.set_xlabel('x distance (km)')
    ax.set_ylabel('y distance (km)')
    ax.set_title(f"\u03B6$_{{10m}}$, \u03B8\u1D68'$_{{10m}}$ = -2 K, $w_{{500m}}$ = 10 m/s", fontsize=14)
    
    plt.suptitle(f"$t={time[0]:.0f}$ min", fontsize=14)
    plt.show()



if makefig[0]:
    
    fig,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(12,9))
    
    plot_cfill(xf, yf, dbz[0,iz,:,:], 'dbz', ax1, datalims=[0,70], xlims=xlims2, ylims=ylims2)
    ax1.plot(xlims2, [yh[iy],yh[iy]], '-k', linewidth=1)
    ax1.plot([xh[ix],xh[ix]], ylims2, '-k', linewidth=1)
    ax1.contour(xh, yh, winterp[0,iz1,:,:], levels=wlevs1, colors='k', linewidths=1)
    ax1.contour(xh, yh, winterp[0,iz1,:,:], levels=wlevs2, colors='k', linewidths=1, linestyles='-.')
    ax1.contour(xh, yh, zvort[0,iz,:,:], levels=[0.2,0.3], colors='yellow', linewidths=1.2, linestyles='-')
    ax1.quiver(xh[::qix], yh[::qix], uinterp[0,iz,::qix,::qix], vinterp[0,iz,::qix,::qix], scale=500, width=0.003, pivot='middle')
    ax1.set_xlabel('x distance (km)')
    ax1.set_ylabel('y distance (km)')
    
    plot_cfill(xf, yf, thrpert[0,iz,:,:], 'thrpert', ax2, datalims=thlims, xlims=xlims2, ylims=ylims2, cmap=cmocean.cm.curl)
    ax2.plot(xlims2, [yh[iy],yh[iy]], '-k', linewidth=1)
    ax2.plot([xh[ix],xh[ix]], ylims2, '-k', linewidth=1)
    ax2.contour(xh, yh, winterp[0,iz1,:,:], levels=wlevs1, colors='k', linewidths=1, linestyles='-')
    ax2.contour(xh, yh, winterp[0,iz1,:,:], levels=wlevs2, colors='k', linewidths=1, linestyles='-.')
    ax2.contour(xh, yh, zvort[0,iz,:,:], levels=[0.2,0.3], colors='yellow', linewidths=1.2, linestyles='-')
    ax2.contour(xh, yh, qr[0,iz,:,:]*1000, levels=qrlevs, colors='deepskyblue', linewidths=1)
    ax2.quiver(xh[::qix], yh[::qix], uinterp[0,iz,::qix,::qix], vinterp[0,iz,::qix,::qix], scale=500, width=0.003, pivot='middle')
    ax2.set_xlabel('x distance (km)')
    ax2.set_ylabel('y distance (km)')
    
    plot_cfill(xf, yf, prspert[iz,:,:]/100, 'prspert', ax3, datalims=[-20,20], xlims=xlims2, ylims=ylims2)
    ax3.plot(xlims2, [yh[iy],yh[iy]], '-k', linewidth=1)
    ax3.plot([xh[ix],xh[ix]], ylims2, '-k', linewidth=1)
    ax3.contour(xh, yh, winterp[0,iz1,:,:], levels=wlevs1, colors='k', linewidths=1)
    ax3.contour(xh, yh, winterp[0,iz1,:,:], levels=wlevs2, colors='k', linewidths=1, linestyles='-.')
    ax3.contour(xh, yh, zvort[0,iz,:,:], levels=[0.2,0.3], colors='yellow', linewidths=1.2, linestyles='-')
    ax3.quiver(xh[::qix], yh[::qix], uinterp[0,iz,::qix,::qix], vinterp[0,iz,::qix,::qix], scale=500, width=0.003, pivot='middle')
    ax3.set_xlabel('x distance (km)')
    ax3.set_ylabel('y distance (km)')
    
    plot_cfill(xf, yf, divh[0,iz,:,:], 'divh', ax4, datalims=[-0.3,0.3], xlims=xlims2, ylims=ylims2)
    ax4.plot(xlims2, [yh[iy],yh[iy]], '-k', linewidth=1)
    ax4.plot([xh[ix],xh[ix]], ylims2, '-k', linewidth=1)
    ax4.contour(xh, yh, winterp[0,iz1,:,:], levels=wlevs1, colors='k', linewidths=1)
    ax4.contour(xh, yh, winterp[0,iz1,:,:], levels=wlevs2, colors='k', linewidths=1, linestyles='-.')
    ax4.contour(xh, yh, zvort[0,iz,:,:], levels=[0.2,0.3], colors='yellow', linewidths=1.2, linestyles='-')
    ax4.quiver(xh[::qix], yh[::qix], uinterp[0,iz,::qix,::qix], vinterp[0,iz,::qix,::qix], scale=500, width=0.003, pivot='middle')
    ax4.set_xlabel('x distance (km)')
    ax4.set_ylabel('y distance (km)')
    
    plt.suptitle(f"$z={z[iz]:.2f}$ km, $t={time[0]:.0f}$ min", fontsize=14)
    
    if figsave:
        plt.savefig(f"/Volumes/Promise_Pegasus_70TB/general/{fd}/TLV_plots/birdseye_{z[iz]*1000:.0f}m_{time[0]:03.0f}min_{n}.png", dpi=400)
    
    
    
    
    
    
if makefig[1]:
    fig,((ax1,ax2),(ax3,ax4),(ax5,ax6)) = plt.subplots(3,2,figsize=(12,11))
    
    c1 = plot_cfill(xf, zf, thrpert[0,:,iy,:], 'thrpert', ax1, datalims=thlims, xlims=xlims2, ylims=zlims, cmap=cmocean.cm.curl)
    ax1.contour(xh, z, winterp[0,:,iy,:], levels=wlevs1, colors='k', linewidths=1, linestyles='-')
    ax1.contour(xh, z, winterp[0,:,iy,:], levels=wlevs2, colors='k', linewidths=1, linestyles='-.')
    # ax1.contour(xh, z, qr[0,:,iy,:]*1000, levels=qrlevs, colors='deepskyblue', linewidths=1, linestyles='-.')
    ax1.contour(xh, z, zvort[0,:,iy,:], levels=[0.1,0.2,0.3], colors='yellow', linewidths=1, linestyles='-')
    ax1.quiver(xh[::qix], z[::qiz], uinterp[0,::qiz,iy,::qix], winterp[0,::qiz,iy,::qix], scale=qsc, width=0.003, pivot='middle')
    ax1.set_xlabel('x distance (km)')
    ax1.set_ylabel('Height (km)')
    ax1.set_title(f"W-E cross-section, y={yh[iy]:.3f} km", fontsize=14)
    
    c2 = plot_cfill(yf, zf, thrpert[0,:,:,ix], 'thrpert', ax2, datalims=thlims, xlims=ylims2, ylims=zlims, cmap=cmocean.cm.curl)
    ax2.contour(yh, z, winterp[0,:,:,ix], levels=wlevs1, colors='k', linewidths=1, linestyles='-')
    ax2.contour(yh, z, winterp[0,:,:,ix], levels=wlevs2, colors='k', linewidths=1, linestyles='-.')
    # ax2.contour(yh, z, qr[0,:,:,ix]*1000, levels=qrlevs, colors='deepskyblue', linewidths=1, linestyles='-.')
    ax2.contour(yh, z, zvort[0,:,:,ix], levels=[0.1,0.2,0.3], colors='yellow', linewidths=1, linestyles='-')
    ax2.quiver(yh[::qix], z[::qiz], vinterp[0,::qiz,::qix,ix], winterp[0,::qiz,::qix,ix], scale=qsc, width=0.003, pivot='middle')
    ax2.set_xlabel('y distance (km)')
    ax2.set_ylabel('Height (km)')
    ax2.set_title(f"S-N cross-section, x={xh[ix]:.3f} km", fontsize=14)
    
    
    c3 = plot_cfill(xf, zf, prspert[:,iy,:]/100, 'prspert', ax3, datalims=[-20,20], xlims=xlims2, ylims=zlims, cmap=cmocean.cm.balance)
    ax3.contour(xh, z, winterp[0,:,iy,:], levels=wlevs1, colors='k', linewidths=1, linestyles='-')
    ax3.contour(xh, z, winterp[0,:,iy,:], levels=wlevs2, colors='k', linewidths=1, linestyles='-.')
    ax3.contour(xh, z, prspert[:,iy,:]/100, levels=[-10], colors='aqua', linewidths=1, linestyles='-.')
    ax3.contour(xh, z, zvort[0,:,iy,:], levels=[0.1,0.2,0.3], colors='yellow', linewidths=1, linestyles='-')
    ax3.quiver(xh[::qix], z[::qiz], uinterp[0,::qiz,iy,::qix], winterp[0,::qiz,iy,::qix], scale=qsc, width=0.003, pivot='middle')
    ax3.set_xlabel('x distance (km)')
    ax3.set_ylabel('Height (km)')
    
    c4 = plot_cfill(yf, zf, prspert[:,:,ix]/100, 'prspert', ax4, datalims=[-20,20], xlims=ylims2, ylims=zlims, cmap=cmocean.cm.balance)
    ax4.contour(yh, z, winterp[0,:,:,ix], levels=wlevs1, colors='k', linewidths=1, linestyles='-')
    ax4.contour(yh, z, winterp[0,:,:,ix], levels=wlevs2, colors='k', linewidths=1, linestyles='-.')
    ax4.contour(yh, z, prspert[:,:,ix]/100, levels=[-10], colors='aqua', linewidths=1, linestyles='-.')
    ax4.contour(yh, z, zvort[0,:,:,ix], levels=[0.1,0.2,0.3], colors='yellow', linewidths=1, linestyles='-')
    ax4.quiver(yh[::qix], z[::qiz], vinterp[0,::qiz,::qix,ix], winterp[0,::qiz,::qix,ix], scale=qsc, width=0.003, pivot='middle')
    ax4.set_xlabel('y distance (km)')
    ax4.set_ylabel('Height (km)')
    
    
    c5 = plot_cfill(xf, zf, winterp[0,:,iy,:], 'w', ax5, datalims=[-30,30], xlims=xlims2, ylims=zlims, cmap=cmocean.cm.balance)
    ax5.contour(xh, z, zvort[0,:,iy,:], levels=[0.1,0.2,0.3], colors='yellow', linewidths=1, linestyles='-')
    ax5.quiver(xh[::qix], z[::qiz], uinterp[0,::qiz,iy,::qix], winterp[0,::qiz,iy,::qix], scale=qsc, width=0.003, pivot='middle')
    ax5.set_xlabel('x distance (km)')
    ax5.set_ylabel('Height (km)')
    
    c6 = plot_cfill(yf, zf, winterp[0,:,:,ix], 'w', ax6, datalims=[-30,30], xlims=ylims2, ylims=zlims, cmap=cmocean.cm.balance)
    ax6.contour(yh, z, zvort[0,:,:,ix], levels=[0.1,0.2,0.3], colors='yellow', linewidths=1, linestyles='-')
    ax6.quiver(yh[::qix], z[::qiz], vinterp[0,::qiz,::qix,ix], winterp[0,::qiz,::qix,ix], scale=qsc, width=0.003, pivot='middle')
    ax6.set_xlabel('y distance (km)')
    ax6.set_ylabel('Height (km)')
    
    plt.suptitle(f"Contours: \u03B6 (yellow), $w$ (black) | $t={time[0]:.0f}$ min", fontsize=14)
    
    plt.show()
    if figsave:
        plt.savefig(f"/Volumes/Promise_Pegasus_70TB/general/{fd}/TLV_plots/vslice_{time[0]:03.0f}min_{n}.png", dpi=400)
    


if makefig[2]:
    fig,((ax1,ax2),(ax3,ax4),(ax5,ax6)) = plt.subplots(3,2,figsize=(12,11))
    
    c1 = plot_cfill(xf, zf, p_dn[:,iy,:]/100, 'prspert', ax1, datalims=pplims, 
                    xlims=xlims2, ylims=zlims, cmap=cm)
    ax1.contour(xh, z, winterp[0,:,iy,:], levels=wlevs1, colors='k', linewidths=1, linestyles='-')
    ax1.contour(xh, z, winterp[0,:,iy,:], levels=wlevs2, colors='k', linewidths=1, linestyles='-.')
    ax1.contour(xh, z, zvort[0,:,iy,:], levels=[0.1,0.2,0.3], colors='yellow', linewidths=1, linestyles='-')
    ax1.quiver(xh[::qix], z[::qiz], uinterp[0,::qiz,iy,::qix], winterp[0,::qiz,iy,::qix], scale=qsc, width=0.003, pivot='middle')
    ax1.set_xlabel('x distance (km)')
    ax1.set_ylabel('Height (km)')
    ax1.set_title(f"W-E cross-section, y={yh[iy]:.3f} km", fontsize=14)
    
    c2 = plot_cfill(yf, zf, p_dn[:,:,ix]/100, 'prspert', ax2, datalims=pplims, 
                    xlims=ylims2, ylims=zlims, cmap=cm)
    ax2.contour(yh, z, winterp[0,:,:,ix], levels=wlevs1, colors='k', linewidths=1, linestyles='-')
    ax2.contour(yh, z, winterp[0,:,:,ix], levels=wlevs2, colors='k', linewidths=1, linestyles='-.')
    ax2.contour(yh, z, zvort[0,:,:,ix], levels=[0.1,0.2,0.3], colors='yellow', linewidths=1, linestyles='-')
    ax2.quiver(yh[::qix], z[::qiz], vinterp[0,::qiz,::qix,ix], winterp[0,::qiz,::qix,ix], scale=qsc, width=0.003, pivot='middle')
    ax2.set_xlabel('y distance (km)')
    ax2.set_ylabel('Height (km)')
    ax2.set_title(f"S-N cross-section, x={xh[ix]:.3f} km", fontsize=14)
    
    
    c3 = plot_cfill(xf, zf, p_b[:,iy,:]/100, 'prspert', ax3, datalims=pplims, 
                    xlims=xlims2, ylims=zlims, cmap=cm)
    ax3.contour(xh, z, winterp[0,:,iy,:], levels=wlevs1, colors='k', linewidths=1, linestyles='-')
    ax3.contour(xh, z, winterp[0,:,iy,:], levels=wlevs2, colors='k', linewidths=1, linestyles='-.')
    ax3.contour(xh, z, zvort[0,:,iy,:], levels=[0.1,0.2,0.3], colors='yellow', linewidths=1, linestyles='-')
    # ax3.contour(xh, z, gaussian_filter(thrpert[0,:,iy,:],1), levels=[-5,-4,-3,-2,-1], colors='dodgerblue', linewidths=1, linestyles='-')
    # ax3.contour(xh, z, gaussian_filter(thrpert[0,:,iy,:],1), levels=[1,2,3,4,5], colors='orange', linewidths=1)
    ax3.quiver(xh[::qix], z[::qiz], uinterp[0,::qiz,iy,::qix], winterp[0,::qiz,iy,::qix], scale=qsc, width=0.003, pivot='middle')
    ax3.set_xlabel('x distance (km)')
    ax3.set_ylabel('Height (km)')
    
    c4 = plot_cfill(yf, zf, p_b[:,:,ix]/100, 'prspert', ax4, datalims=pplims, 
                    xlims=ylims2, ylims=zlims, cmap=cm)
    ax4.contour(yh, z, winterp[0,:,:,ix], levels=wlevs1, colors='k', linewidths=1, linestyles='-')
    ax4.contour(yh, z, winterp[0,:,:,ix], levels=wlevs2, colors='k', linewidths=1, linestyles='-.')
    ax4.contour(yh, z, zvort[0,:,:,ix], levels=[0.1,0.2,0.3], colors='yellow', linewidths=1, linestyles='-')
    # ax4.contour(yh, z, gaussian_filter(thrpert[0,:,:,ix],1), levels=[-5,-4,-3,-2,-1], colors='dodgerblue', linewidths=1, linestyles='-')
    # ax4.contour(yh, z, gaussian_filter(thrpert[0,:,:,ix],1), levels=[1,2,3,4,5], colors='orange', linewidths=1)
    ax4.quiver(yh[::qix], z[::qiz], vinterp[0,::qiz,::qix,ix], winterp[0,::qiz,::qix,ix], scale=qsc, width=0.003, pivot='middle')
    ax4.set_xlabel('y distance (km)')
    ax4.set_ylabel('Height (km)')
    
    c5 = plot_cfill(xf, zf, p_dl[:,iy,:]/100, 'prspert', ax5, datalims=pplims, 
                    xlims=xlims2, ylims=zlims, cmap=cm)
    ax5.contour(xh, z, winterp[0,:,iy,:], levels=wlevs1, colors='k', linewidths=1, linestyles='-')
    ax5.contour(xh, z, winterp[0,:,iy,:], levels=wlevs2, colors='k', linewidths=1, linestyles='-.')
    ax5.contour(xh, z, zvort[0,:,iy,:], levels=[0.1,0.2,0.3], colors='yellow', linewidths=1, linestyles='-')
    ax5.quiver(xh[::qix], z[::qiz], uinterp[0,::qiz,iy,::qix], winterp[0,::qiz,iy,::qix], scale=qsc, width=0.003, pivot='middle')
    ax5.set_xlabel('x distance (km)')
    ax5.set_ylabel('Height (km)')
    
    c6 = plot_cfill(yf, zf, p_dl[:,:,ix]/100, 'prspert', ax6, datalims=pplims, 
                    xlims=ylims2, ylims=zlims, cmap=cm)
    ax6.contour(yh, z, winterp[0,:,:,ix], levels=wlevs1, colors='k', linewidths=1, linestyles='-')
    ax6.contour(yh, z, winterp[0,:,:,ix], levels=wlevs2, colors='k', linewidths=1, linestyles='-.')
    ax6.contour(yh, z, zvort[0,:,:,ix], levels=[0.1,0.2,0.3], colors='yellow', linewidths=1, linestyles='-')
    ax6.quiver(yh[::qix], z[::qiz], vinterp[0,::qiz,::qix,ix], winterp[0,::qiz,::qix,ix], scale=qsc, width=0.003, pivot='middle')
    ax6.set_xlabel('y distance (km)')
    ax6.set_ylabel('Height (km)')
    
    plt.suptitle(f"Nonlinear dyn, buoyancy, & linear dyn pp, \u03B6 (yellow), $w$ (black) | $t={time[0]:.0f}$ min", fontsize=14)
    
    plt.show()
    if figsave:
        plt.savefig(f"/Volumes/Promise_Pegasus_70TB/general/{fd}/TLV_plots/ppnl_vslice_{time[0]:03.0f}min_{n}.png", dpi=400)
    
    
if makefig[3]:
    fig,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(12,9))
    
    plot_cfill(xf, yf, p_dn[iz,:,:]/100, 'prspert', ax1, datalims=pplims, 
               xlims=xlims2, ylims=ylims2, cmap=cm)
    ax1.contour(xh, yh, winterp[0,iz1,:,:], levels=wlevs2, colors='k', linewidths=1, linestyles='-.')
    ax1.contour(xh, yh, winterp[0,iz1,:,:], levels=wlevs1, colors='k', linewidths=1, linestyles='-')
    ax1.contour(xh, yh, zvort[0,iz,:,:], levels=[0.1,0.2,0.3], colors='yellow', linewidths=1)
    ax1.quiver(xh[::qix], yh[::qix], uinterp[0,iz,::qix,::qix], vinterp[0,iz,::qix,::qix], scale=500, width=0.003, pivot='middle')
    ax1.set_title('Nonlinear dynamic')
    ax1.set_ylabel('y distance (km)')
    
    plot_cfill(xf, yf, p_b[iz,:,:]/100, 'prspert', ax2, datalims=pplims, 
               xlims=xlims2, ylims=ylims2, cmap=cm)
    ax2.contour(xh, yh, winterp[0,iz1,:,:], levels=wlevs2, colors='k', linewidths=1, linestyles='-.')
    ax2.contour(xh, yh, winterp[0,iz1,:,:], levels=wlevs1, colors='k', linewidths=1, linestyles='-')
    ax2.contour(xh, yh, zvort[0,iz,:,:], levels=[0.1,0.2,0.3], colors='yellow', linewidths=1)
    ax2.quiver(xh[::qix], yh[::qix], uinterp[0,iz,::qix,::qix], vinterp[0,iz,::qix,::qix], scale=500, width=0.003, pivot='middle')
    ax2.set_title('Buoyancy')
    ax2.set_ylabel('y distance (km)')
    
    plot_cfill(xf, yf, p_dl[iz,:,:]/100, 'prspert', ax3, datalims=pplims, 
               xlims=xlims2, ylims=ylims2, cmap=cm)
    ax3.contour(xh, yh, winterp[0,iz1,:,:], levels=wlevs2, colors='k', linewidths=1, linestyles='-.')
    ax3.contour(xh, yh, winterp[0,iz1,:,:], levels=wlevs1, colors='k', linewidths=1, linestyles='-')
    ax3.contour(xh, yh, zvort[0,iz,:,:], levels=[0.1,0.2,0.3], colors='yellow', linewidths=1)
    ax3.quiver(xh[::qix], yh[::qix], uinterp[0,iz,::qix,::qix], vinterp[0,iz,::qix,::qix], scale=500, width=0.003, pivot='middle')
    ax3.set_title('Linear dynamic')
    ax3.set_xlabel('x distance (km)')
    ax3.set_ylabel('y distance (km)')
    
    plot_cfill(xf, yf, prspert[iz,:,:]/100, 'prspert', ax4, datalims=pplims, 
               xlims=xlims2, ylims=ylims2, cmap=cm)
    ax4.contour(xh, yh, winterp[0,iz1,:,:], levels=wlevs2, colors='k', linewidths=1, linestyles='-.')
    ax4.contour(xh, yh, winterp[0,iz1,:,:], levels=wlevs1, colors='k', linewidths=1, linestyles='-')
    ax4.contour(xh, yh, zvort[0,iz,:,:], levels=[0.1,0.2,0.3], colors='yellow', linewidths=1)
    ax4.quiver(xh[::qix], yh[::qix], uinterp[0,iz,::qix,::qix], vinterp[0,iz,::qix,::qix], scale=500, width=0.003, pivot='middle')
    ax4.set_title('Total')
    ax4.set_xlabel('x distance (km)')
    ax4.set_ylabel('y distance (km)')
    
    plt.suptitle(f"\u25BD$^2$p' forcing, $w$ (black) | $z={z[iz]:.2f}$ km, $t={time[0]:.0f}$ min", fontsize=14)
    
    plt.show()
    if figsave:
        # plt.savefig(f"/Volumes/Promise_Pegasus_70TB/general/{fd}/TLV_plots/ppsfc_{time[0]:03.0f}min_{n}.png", dpi=400)
        plt.savefig(f"/Volumes/Promise_Pegasus_70TB/general/{fd}/TLV_plots/pp{z[iz]*1000:.0f}m_{time[0]:03.0f}min_{n}.png", dpi=400)











