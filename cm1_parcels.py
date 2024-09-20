#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 14:34:11 2023

@author: morgan.schneider
"""

####################
### Load modules ###
####################

from CM1utils import *
from matplotlib import patches

#%% Load parcel file

# Parcel data dimensions (time, pid)

fp = '/Volumes/Promise_Pegasus_70TB/merger/qlcs-125m/'
ip = '/Users/morgan.schneider/Documents/merger/qlcs-125m/'

# Read parcel data
ds = nc.Dataset(fp+'cm1out_pdata.nc')
ptime = ds.variables['time'][:].data
pid = ds.variables['xh'][:].data
x = ds.variables['x'][:].data
y = ds.variables['y'][:].data
z = ds.variables['z'][:].data
# u = ds.variables['u'][:].data
# v = ds.variables['v'][:].data
w = ds.variables['w'][:].data
# th = ds.variables['th'][:].data
# prs = ds.variables['prs'][:].data
# qv = ds.variables['qv'][:].data
# thr = ds.variables['th'][:].data * (1 + 0.61*ds.variables['qv'][:].data - 
#             (ds.variables['qc'][:].data + ds.variables['qr'][:].data + ds.variables['qi'][:].data
#             + ds.variables['qs'][:].data + ds.variables['qg'][:].data + ds.variables['qhl'][:].data))
# qc = ds.variables['qc'][:].data
# qr = ds.variables['qr'][:].data
# qi = ds.variables['qi'][:].data
# qs = ds.variables['qs'][:].data
# qg = ds.variables['qg'][:].data
# qhl = ds.variables['qhl'][:].data
# dbz = ds.variables['dbz'][:].data
b = ds.variables['b'][:].data
# vpg = ds.variables['vpg'][:].data
zvort = ds.variables['zvort'][:].data
ds.close()


#%% Load data and filter trajectories

# Read gridded output
fnum = 72
ds = nc.Dataset(fp+f"cm1out_{fnum:06d}.nc")
stime = ds.variables['time'][:].data[0]
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
zh = ds.variables['z'][:].data
xf = ds.variables['xf'][:].data
yf = ds.variables['yf'][:].data
zf = ds.variables['zf'][:].data

iz1 = np.where(zh >= 1)[0][0]
dbz = ds.variables['dbz'][:].data[0,0,:,:]
# winterp = ds.variables['winterp'][:].data[0,iz1,:,:]
ds.close()

ds = nc.Dataset(fp+"base/cm1out_000013.nc")
dbz0 = ds.variables['dbz2'][:].data[0,0,:,:]
ds.close()

# Filter trajectories

# Simple DBZ with box for the parcels
# Rework this to plot parcels in the 210 min box over time
# Can probably save this stuff to a pickle file too

ti = np.where(ptime == stime)[0][0]

box_name = 'QLCS'

dbfile = open(ip+'boxes_q.pkl', 'rb')
box = pickle.load(dbfile)
x1 = box['x1_pp'][fnum-14]
x2 = box['x2_pp'][fnum-14]
y1 = box['y1_pp'][fnum-14]
y2 = box['y2_pp'][fnum-14]
dbfile.close()

# z_cond = (z[ti,:]<=1000) & (z[ti,:]>=100) & (x[ti,:]>=x1*1000) & (x[ti,:]<=x2*1000) & (y[ti,:]>=y1*1000) & (y[ti,:]<=y2*1000)
# w_cond = (w[ti,:]>=1) & (z[ti,:]<=1000) & (z[ti,:]>=100) & (x[ti,:]>=x1*1000) & (x[ti,:]<=x2*1000) & (y[ti,:]>=y1*1000) & (y[ti,:]<=y2*1000)
# vort_cond = (zvort[ti,:]>=0.01) & (z[ti,:]<=1000) & (z[ti,:]>=100) & (x[ti,:]>=x1*1000) & (x[ti,:]<=x2*1000) & (y[ti,:]>=y1*1000) & (y[ti,:]<=y2*1000)
wvort_cond = (zvort[ti,:]>=0.01) & (w[ti,:]>=5) & (z[ti,:]<=1000) & (z[ti,:]>=100) & (x[ti,:]>=x1*1000) & (x[ti,:]<=x2*1000) & (y[ti,:]>=y1*1000) & (y[ti,:]<=y2*1000)
# z_cond_str = "z<1000"
# w_cond_str = "w>1, z<1000"
# vort_cond_str = "\u03B6>0.01, z<1000"
wvort_cond_str = "\u03B6>0.01, w>5, z<1000"

pids = pid[wvort_cond]
x_mv = x[0:ti+1, wvort_cond]
y_mv = y[0:ti+1, wvort_cond]
z_mv = z[0:ti+1, wvort_cond]
w_mv = w[0:ti+1, wvort_cond]
zvort_mv = zvort[0:ti+1, wvort_cond]
b_mv = b[0:ti+1, wvort_cond]

cond_str = wvort_cond_str



#% Set source regions (need to combine this with the other code later)

if np.equal(box_name, 'MERGER (MV1)'):
    dbfile = open(ip+f"traj_clusters_{stime/60:.0f}min.pkl", 'rb')
    ccs = pickle.load(dbfile)
    cc = ccs['mv1']
    dbfile.close()
    
    # cc[(ccs['mv1_w5']==1)] = 2 # can't tell if this source is more LL or ML?
    
    c0 = (cc == 0) # low-level inflow
    c1 = (cc == 1) # supercell cold pool
    c2 = (cc == 2) # mid-level inflow
    c3 = (cc == 3) # QLCS outflow
    c4 = (cc == 4) # parcels starting in mesocyclone
    
    
    # for i in range(max(cc['mv1_w5'])+1):
    #     exec(f"c{i+1} = (cc['mv1_w5'] == i)")
    
    # c1 = (labels1 == 0)
    # c2 = (labels1 == 1)
    # c3 = (labels1 == 2)
    # c4 = (labels1 == 3)
    
    print(f"{len(pids[c0])} parcels from near-sfc inflow ({len(pids[c0])/len(pids)*100:.02f}%)")
    print(f"{len(pids[c1])} parcels from supercell cold pool ({len(pids[c1])/len(pids)*100:.02f}%)")
    print(f"{len(pids[c2])} parcels from mid-level inflow ({len(pids[c2])/len(pids)*100:.02f}%)")
    print(f"{len(pids[c3])} parcels from QLCS outflow ({len(pids[c3])/len(pids)*100:.02f}%)")
    print(f"{len(pids[c4])} parcels starting in updraft ({len(pids[c4])/len(pids)*100:.02f}%)")
    print(f"{len(pids)} total parcels")
    
    img_str = 'MV1'
    xl = [-50,30]
    yl = [-120,-40]


if np.equal(box_name, 'MERGER (MV2)'):
    dbfile = open(ip+f"traj_clusters_{stime/60:.0f}min.pkl", 'rb')
    ccs = pickle.load(dbfile)
    cc = ccs['mv2']
    dbfile.close()
    
    c0 = (cc == 0) # low-level inflow
    c1 = (cc == 1) # supercell cold pool
    c2 = (cc == 2) # mid-level inflow
    c3 = (cc == 3) # QLCS low-level outflow
    c4 = (cc == 4) # QLCS outflow aloft
    
    print(f"{len(pids[c0])} parcels from near-sfc inflow ({len(pids[c0])/len(pids)*100:.02f}%)")
    print(f"{len(pids[c1])} parcels from supercell cold pool ({len(pids[c1])/len(pids)*100:.02f}%)")
    print(f"{len(pids[c2])} parcels from mid-level inflow ({len(pids[c2])/len(pids)*100:.02f}%)")
    print(f"{len(pids[c3])} parcels from QLCS low-level outflow ({len(pids[c3])/len(pids)*100:.02f}%)")
    print(f"{len(pids[c4])} parcels from QLCS outflow aloft ({len(pids[c4])/len(pids)*100:.02f}%)")
    print(f"{len(pids)} total parcels")
    
    img_str = 'MV2'
    xl = [-50,30]
    yl = [-100,-20]


if np.equal(box_name, 'QLCS'):
    dbfile = open(ip+f"traj_clusters_{stime/60:.0f}min.pkl", 'rb')
    ccs = pickle.load(dbfile)
    cc = ccs['q']
    dbfile.close()
    
    c0 = (cc == 0) # low-level inflow
    c1 = (cc == 1) # QLCS outflow
    
    print(f"{len(pids[c0])} parcels from low-level inflow ({len(pids[c0])/len(pids)*100:.02f}%)")
    print(f"{len(pids[c1])} parcels from QLCS outflow ({len(pids[c1])/len(pids)*100:.02f}%)")
    
    img_str = 'Q'
    xl = [-50,30]
    yl = [-115,-35]
    

#%% Save source regions to pkl

if False:
    dbfile = open(ip+f"traj_clusters_{stime/60:.0f}min.pkl", 'wb')
    ccs = {'q':cc_new}
    pickle.dump(ccs, dbfile)
    dbfile.close()

if False:
    dbfile = open(ip+f"traj_clusters_{stime/60:.0f}min.pkl", 'rb')
    ccs = pickle.load(dbfile)
    dbfile.close()
    
    new_vars = {'q':cc_new}
    ccs.update(new_vars)
    dbfile = open(ip+f"traj_clusters_{stime/60:.0f}min.pkl", 'wb')
    pickle.dump(ccs, dbfile)
    dbfile.close()


if True:
    cc = labels1
    cc_new = np.zeros(shape=cc.shape, dtype=cc.dtype)
    cc_new[(cc == 0)] = 1
    cc_new[(cc == 1)] = 0
    # cc_new[(cc == 2)] = 1
    # cc_new[(labels1 == 3)] = 1
    # cc_new[(cc == 4)] = 3
    
    # cc_new[(cc==2) & (x_mv[0,:]<-10000) & (y_mv[0,:]>-110000)] = 0
    # cc_new[(z_mv[0,:]>1100)] = 2
    # cc_new[(x_mv[0,:]<-30000)] = 3
    # cc_new[(cc_new == 1) & (b_mv[0,:] > -0.03)] = 0
    # cc_new[(cc_new == 4) & (b_mv[0,:] < -0.01) & (y_mv[0,:] > -105000)] = 1
    # cc_new[(labels1 == 1) & (x_mv[0,:] > -10000)] = 0
    
    cc_old = cc
    cc = cc_new
    
    


#%% Make plots

# from scipy.ndimage import gaussian_filter

qix = int(np.round(len(pids)/90))

figsave = 0

xl = [-50,30]
yl = [-110,-30]

# Plot trajectories over dbz
if False:
    qix = 20
    
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    plot_cfill(xf, yf, np.ma.masked_array(dbz, dbz<1), 'dbz', ax, datalims=[1,70], xlims=xl, ylims=yl)
    ax.set_xlabel('E-W distance (km)')
    ax.set_ylabel('N-S distance (km)')
    ax.set_title(f"{box_name} - 10 m $Z_H$ at {stime/60:.0f} min \n {cond_str}")
    ax.plot(x_mv[:,::qix]/1000, y_mv[:,::qix]/1000, '-k', linewidth=1)
    if figsave:
        plt.savefig(ip+f"imgs_trajectories/{img_str}_traj_{stime/60:.0f}min.png", dpi=300)

# Plot trajectories colored by z, w, zvort, buoyancy
if False:
    qix = 1
    
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    p = ax.scatter(x_mv[:,::qix]/1000, y_mv[:,::qix]/1000, s=1, c=z_mv[:,::qix]/1000, marker='.', cmap='pyart_HomeyerRainbow', vmin=0, vmax=3.5)
    ax.contour(xh, yh, dbz, levels=[30], colors='k', linewidths=1)
    ax.set_xlabel('E-W distance (km)')
    ax.set_ylabel('N-S distance (km)')
    ax.set_xlim(xl)
    ax.set_ylim(yl)
    ax.set_title(f"{box_name} - Trajectories by $z$ at {stime/60:.0f} min \n {cond_str}")
    plt.colorbar(p, ax=ax)
    if figsave:
        plt.savefig(ip+f"imgs_trajectories/{img_str}_traj-z_{stime/60:.0f}min.png", dpi=300)
    
    
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    p = ax.scatter(x_mv[:,::qix]/1000, y_mv[:,::qix]/1000, s=1, c=w_mv[:,::qix], marker='.', cmap='RdBu_r', vmin=-10, vmax=10)
    ax.contour(xh, yh, dbz, levels=[30], colors='k', linewidths=1)
    ax.set_xlabel('E-W distance (km)')
    ax.set_ylabel('N-S distance (km)')
    ax.set_xlim(xl)
    ax.set_ylim(yl)
    ax.set_title(f"{box_name} - Trajectories by $w$ at {stime/60:.0f} min \n {cond_str}")
    plt.colorbar(p, ax=ax)
    if figsave:
        plt.savefig(ip+f"imgs_trajectories/{img_str}_traj-w_{stime/60:.0f}min.png", dpi=300)
    
    
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    p = ax.scatter(x_mv[:,::qix]/1000, y_mv[:,::qix]/1000, s=1, c=zvort_mv[:,::qix], marker='.', cmap='RdBu_r', vmin=-0.04, vmax=0.04)
    ax.contour(xh, yh, dbz, levels=[30], colors='k', linewidths=1)
    ax.set_xlabel('E-W distance (km)')
    ax.set_ylabel('N-S distance (km)')
    ax.set_xlim(xl)
    ax.set_ylim(yl)
    ax.set_title(f"{box_name} - Trajectories by zvort at {stime/60:.0f} min \n {cond_str}")
    plt.colorbar(p, ax=ax)
    if figsave:
        plt.savefig(ip+f"imgs_trajectories/{img_str}_traj-zvort_{stime/60:.0f}min.png", dpi=300)
        
    
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    p = ax.scatter(x_mv[:,::qix]/1000, y_mv[:,::qix]/1000, s=1, c=b_mv[:,::qix], marker='.', cmap='pyart_HomeyerRainbow', vmin=-0.4, vmax=0.1)
    ax.contour(xh, yh, dbz, levels=[30], colors='k', linewidths=1)
    ax.set_xlim(xl)
    ax.set_ylim(yl)
    ax.set_title(f"{box_name} - Buoyancy, 30 dBZ at {stime/60:.0f} min \n {cond_str}")
    plt.colorbar(p, ax=ax)
    ax.set_xlabel('E-W distance (km)')
    ax.set_ylabel('N-S distance (km)')
    # ax.set_title(f"{box_name} - Source A - 10 m $Z_H$ at {stime/60:.0f} min \n {cond_str}")
    # ax.plot(x_src[:,::qix]/1000, y_src[:,::qix]/1000, '-k', linewidth=1)
    # plt.show()
    if figsave:
        plt.savefig(ip+f"imgs_trajectories/{img_str}_traj-buoy_{stime/60:.0f}min.png", dpi=300)
    
    plt.show()

# Plot trajectory starts by path buoyancy, path height, initial buoyancy, initial height, initial w
if True:
    qix = 1
    
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    p = ax.scatter(x_mv[:,::qix]/1000, y_mv[:,::qix]/1000, s=1, c=b_mv[:,::qix], marker='.', cmap='pyart_HomeyerRainbow', vmin=-0.2, vmax=0.1)
    ax.contour(xh, yh, dbz, levels=[30], colors='k', linewidths=1)
    # ax.contour(xh, yh, gaussian_filter(winterp,2), levels=[-5,5], colors='b', linewidths=1)
    ax.set_xlim(xl)
    ax.set_ylim(yl)
    ax.set_title(f"{box_name} - Buoyancy, 30 dBZ at {stime/60:.0f} min")
    plt.colorbar(p, ax=ax)
    ax.set_xlabel('E-W distance (km)')
    ax.set_ylabel('N-S distance (km)')
    # ax.set_title(f"{box_name} - Source A - 10 m $Z_H$ at {stime/60:.0f} min \n {cond_str}")
    # ax.plot(x_src[:,::qix]/1000, y_src[:,::qix]/1000, '-k', linewidth=1)
    # plt.show()
    if figsave:
        plt.savefig(ip+f"imgs_trajectories/{img_str}_traj-buoy_{stime/60:.0f}min.png", dpi=300)
        
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    p = ax.scatter(x_mv[:,::qix]/1000, y_mv[:,::qix]/1000, s=1, c=z_mv[:,::qix]/1000, marker='.', cmap='pyart_HomeyerRainbow', vmin=0, vmax=2)
    ax.contour(xh, yh, dbz, levels=[30], colors='k', linewidths=1)
    # ax.contour(xh, yh, gaussian_filter(winterp,2), levels=[-5,5], colors='b', linewidths=1)
    ax.set_xlim(xl)
    ax.set_ylim(yl)
    ax.set_title(f"{box_name} - Height, 30 dBZ at {stime/60:.0f} min")
    plt.colorbar(p, ax=ax)
    ax.set_xlabel('E-W distance (km)')
    ax.set_ylabel('N-S distance (km)')
    # ax.set_title(f"{box_name} - Source A - 10 m $Z_H$ at {stime/60:.0f} min \n {cond_str}")
    # ax.plot(x_src[:,::qix]/1000, y_src[:,::qix]/1000, '-k', linewidth=1)
    # plt.show()
    if figsave:
        plt.savefig(ip+f"imgs_trajectories/{img_str}_traj-z_{stime/60:.0f}min.png", dpi=300)
    
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    p = ax.scatter(x_mv[0,:]/1000, y_mv[0,:]/1000, s=5, c=b_mv[0,:], marker='.', cmap='pyart_HomeyerRainbow', vmin=-0.2, vmax=0.1)
    ax.contour(xh, yh, dbz0, levels=[30], colors='k', linewidths=1)
    ax.set_xlim(xl)
    ax.set_ylim(yl)
    ax.set_title(f"{box_name} - Starting buoyancy for parcels at {stime/60:.0f} min")
    plt.colorbar(p, ax=ax)
    ax.set_xlabel('E-W distance (km)')
    ax.set_ylabel('N-S distance (km)')
    if figsave:
        plt.savefig(ip+f"imgs_trajectories/{img_str}_traj-buoystart_{stime/60:.0f}min.png", dpi=300)
    
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    p = ax.scatter(x_mv[0,:]/1000, y_mv[0,:]/1000, s=5, c=z_mv[0,:]/1000, marker='.', cmap='pyart_HomeyerRainbow', vmin=0, vmax=2)
    ax.contour(xh, yh, dbz0, levels=[30], colors='k', linewidths=1)
    ax.set_xlim(xl)
    ax.set_ylim(yl)
    ax.set_title(f"{box_name} - Starting height for parcels at {stime/60:.0f} min")
    plt.colorbar(p, ax=ax)
    ax.set_xlabel('E-W distance (km)')
    ax.set_ylabel('N-S distance (km)')
    if figsave:
        plt.savefig(ip+f"imgs_trajectories/{img_str}_traj-zstart_{stime/60:.0f}min.png", dpi=300)
    
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    p = ax.scatter(x_mv[0,:]/1000, y_mv[0,:]/1000, s=5, c=w_mv[0,:], marker='.', cmap='pyart_HomeyerRainbow', vmin=-10, vmax=10)
    ax.contour(xh, yh, dbz0, levels=[30], colors='k', linewidths=1)
    ax.set_xlim(xl)
    ax.set_ylim(yl)
    ax.set_title(f"{box_name} - Starting vertical velocity for parcels at {stime/60:.0f} min")
    plt.colorbar(p, ax=ax)
    ax.set_xlabel('E-W distance (km)')
    ax.set_ylabel('N-S distance (km)')
    if figsave:
        plt.savefig(ip+f"imgs_trajectories/{img_str}_traj-wstart_{stime/60:.0f}min.png", dpi=300)
    
    plt.show()

    
# Plot trajectories colored by source
if False:
    fig,(ax1,ax2) = plt.subplots(1,2,figsize=(12,5))
    ax1.scatter(x_mv/1000, y_mv/1000, s=1, c=np.tile(cc, (len(x_mv[:,0]),1)), marker='.', cmap='pyart_HomeyerRainbow', vmin=0, vmax=4)
    ax1.contour(xh, yh, dbz, levels=[30], colors='k', linewidths=1)
    ax1.set_xlim(xl)
    ax1.set_ylim(yl)
    ax1.set_title(f"Trajectory clusters, 30 dBZ")
    p = ax2.scatter(x_mv[0,:]/1000, y_mv[0,:]/1000, s=1, c=cc, marker='.', cmap='pyart_HomeyerRainbow', vmin=0, vmax=4)
    ax2.contour(xh, yh, dbz0, levels=[30], colors='k', linewidths=1)
    # ax2.contour(xh, yh, winterp, levels=[5], colors='r')
    ax2.set_xlim(xl)
    ax2.set_ylim(yl)
    ax2.set_title(f"Trajectory starting points, initial 30 dBZ")
    plt.colorbar(p, ax=ax2)
    plt.suptitle(f"{box_name} - Trajectory clusters at {stime/60:.0f} min - {cond_str}")
    # if figsave:
    #     plt.savefig(ip+f"imgs_trajectories/{img_str}_traj-clusters_{stime/60:.0f}min.png", dpi=300)
    plt.show()


# Plot trajectories by source over dbz - single source
if False:
    c = (cc_new == 1)
    x_src = x_mv[:,c]
    y_src = y_mv[:,c]
    qix = 1
    
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    plot_cfill(xf, yf, np.ma.masked_array(dbz, dbz<1), 'dbz', ax, datalims=[1,70], xlims=xl, ylims=yl)
    ax.set_xlabel('E-W distance (km)')
    ax.set_ylabel('N-S distance (km)')
    ax.set_title(f"{box_name} - Source A - 10 m $Z_H$ at {stime/60:.0f} min \n {cond_str}")
    ax.plot(x_src[:,::qix]/1000, y_src[:,::qix]/1000, '-k', linewidth=1)
    plt.show()
    if figsave:
        plt.savefig(ip+f"imgs_trajectories/{img_str}_src1_traj_{stime/60:.0f}min.png", dpi=300)

# Plot trajectories by source colored by z, w, zvort, buoyancy
if False:
    c = (cc_new == 1)
    x_src = x_mv[:,c]
    y_src = y_mv[:,c]
    z_src = z_mv[:,c]
    w_src = w_mv[:,c]
    zvort_src = zvort_mv[:,c]
    b_src = b_mv[:,c]
    qix = 1
    
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    p = ax.scatter(x_src[:,::qix]/1000, y_src[:,::qix]/1000, s=1, c=z_src[:,::qix]/1000, marker='.', cmap='pyart_HomeyerRainbow', vmin=0, vmax=3.5)
    ax.contour(xh, yh, dbz0, levels=[30], colors='k', linewidths=1)
    ax.set_xlabel('E-W distance (km)')
    ax.set_ylabel('N-S distance (km)')
    ax.set_xlim(xl)
    ax.set_ylim(yl)
    ax.set_title(f"{box_name} - Trajectories by $z$ at {stime/60:.0f} min \n {cond_str}")
    plt.colorbar(p, ax=ax)
    if figsave:
        plt.savefig(ip+f"imgs_trajectories/{img_str}_src1_traj-z_{stime/60:.0f}min.png", dpi=300)
    
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    p = ax.scatter(x_src[:,::qix]/1000, y_src[:,::qix]/1000, s=1, c=w_src[:,::qix], marker='.', cmap='RdBu_r', vmin=-10, vmax=10)
    ax.contour(xh, yh, dbz, levels=[30], colors='k', linewidths=1)
    ax.set_xlabel('E-W distance (km)')
    ax.set_ylabel('N-S distance (km)')
    ax.set_xlim(xl)
    ax.set_ylim(yl)
    ax.set_title(f"{box_name} - Trajectories by $w$ at {stime/60:.0f} min \n {cond_str}")
    plt.colorbar(p, ax=ax)
    if figsave:
        plt.savefig(ip+f"imgs_trajectories/{img_str}_src1_traj-w_{stime/60:.0f}min.png", dpi=300)
    
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    p = ax.scatter(x_src[:,::qix]/1000, y_src[:,::qix]/1000, s=1, c=zvort_src[:,::qix], marker='.', cmap='RdBu_r', vmin=-0.04, vmax=0.04)
    ax.contour(xh, yh, dbz, levels=[30], colors='k', linewidths=1)
    ax.set_xlabel('E-W distance (km)')
    ax.set_ylabel('N-S distance (km)')
    ax.set_xlim(xl)
    ax.set_ylim(yl)
    ax.set_title(f"{box_name} - Trajectories by zvort at {stime/60:.0f} min \n {cond_str}")
    plt.colorbar(p, ax=ax)
    if figsave:
        plt.savefig(ip+f"imgs_trajectories/{img_str}_src1_traj-zvort_{stime/60:.0f}min.png", dpi=300)
    
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    p = ax.scatter(x_src[:,::qix]/1000, y_src[:,::qix]/1000, s=1, c=b_src[:,::qix], marker='.', cmap='pyart_HomeyerRainbow', vmin=-0.4, vmax=0.1)
    ax.contour(xh, yh, dbz, levels=[30], colors='k', linewidths=1)
    ax.set_xlabel('E-W distance (km)')
    ax.set_ylabel('N-S distance (km)')
    ax.set_xlim(xl)
    ax.set_ylim(yl)
    ax.set_title(f"{box_name} - Trajectories by buoyancy at {stime/60:.0f} min \n {cond_str}")
    plt.colorbar(p, ax=ax)
    if figsave:
        plt.savefig(ip+f"imgs_trajectories/{img_str}_src1_traj-buoy_{stime/60:.0f}min.png", dpi=300)



# Plot trajectories at various times over dbz and colored by z, w, zvort, buoyancy
if False:
    fn = 43
    ds = nc.Dataset(fp+f"cm1out_{fn:06d}.nc")
    stime2 = ds.variables['time'][:].data[0]
    dbz2 = ds.variables['dbz'][:].data[0,0,:,:]
    # w_field = ds.variables['winterp'][:].data[0,(zh>=1),:,:]
    # zvort_field = ds.variables['zvort'][:].data[0,0,:,:]
    ds.close()
    
    tj = np.where(ptime == stime2)[0][0]
    
    # Trajectories plotted over dbz
    if True:
        qix = 20
        
        fig,ax = plt.subplots(1,1,figsize=(8,6))
        plot_cfill(xf, yf, np.ma.masked_array(dbz2, dbz2<1), 'dbz', ax, datalims=[1,70], xlims=xl, ylims=yl)
        ax.set_xlabel('E-W distance (km)')
        ax.set_ylabel('N-S distance (km)')
        ax.set_title(f"{box_name} - 10 m $Z_H$ at {stime2/60:.0f} min \n {cond_str}")
        ax.plot(x_mv[0:tj+1,::qix]/1000, y_mv[0:tj+1,::qix]/1000, '-k', linewidth=1)
        plt.show()
        if figsave:
            plt.savefig(ip+f"imgs_trajectories/{img_str}_traj_{stime2/60:.0f}min.png", dpi=300)
    
    # Trajectories colored by z, w, zvort, buoyancy
    if False:
        qix = 1
        
        fig,ax = plt.subplots(1,1,figsize=(8,6))
        p = ax.scatter(x_mv[0:tj+1,::qix]/1000, y_mv[0:tj+1,::qix]/1000, s=1, c=z_mv[0:tj+1,::qix]/1000, marker='.', cmap='pyart_HomeyerRainbow', vmin=0, vmax=3.5)
        ax.set_xlabel('E-W distance (km)')
        ax.set_ylabel('N-S distance (km)')
        ax.set_xlim(xl)
        ax.set_ylim(yl)
        ax.set_title(f"{box_name} - Trajectories by $z$ at {stime2/60:.0f} min \n {cond_str}")
        plt.show()
        if figsave:
            plt.savefig(ip+f"imgs_trajectories/{img_str}_traj-z_{stime2/60:.0f}min.png", dpi=300)
        
        fig,ax = plt.subplots(1,1,figsize=(8,6))
        p = ax.scatter(x_mv[0:tj+1,::qix]/1000, y_mv[0:tj+1,::qix]/1000, s=1, c=w_mv[0:tj+1,::qix], marker='.', cmap='RdBu_r', vmin=-10, vmax=10)
        ax.set_xlabel('E-W distance (km)')
        ax.set_ylabel('N-S distance (km)')
        ax.set_xlim(xl)
        ax.set_ylim(yl)
        ax.set_title(f"{box_name} - Trajectories by $w$ at {stime2/60:.0f} min \n {cond_str}")
        plt.show()
        if figsave:
            plt.savefig(ip+f"imgs_trajectories/{img_str}_traj-w_{stime2/60:.0f}min.png", dpi=300)
        
        fig,ax = plt.subplots(1,1,figsize=(8,6))
        p = ax.scatter(x_mv[0:tj+1,::qix]/1000, y_mv[0:tj+1,::qix]/1000, s=1, c=zvort_mv[0:tj+1,::qix], marker='.', cmap='RdBu_r', vmin=-0.04, vmax=0.04)
        ax.set_xlabel('E-W distance (km)')
        ax.set_ylabel('N-S distance (km)')
        ax.set_xlim(xl)
        ax.set_ylim(yl)
        ax.set_title(f"{box_name} - Trajectories by zvort at {stime2/60:.0f} min \n {cond_str}")
        plt.show()
        if figsave:
            plt.savefig(ip+f"imgs_trajectories/{img_str}_traj-zvort_{stime2/60:.0f}min.png", dpi=300)
        
        fig,ax = plt.subplots(1,1,figsize=(8,6))
        p = ax.scatter(x_mv[0:tj+1,::qix]/1000, y_mv[0:tj+1,::qix]/1000, s=1, c=b_mv[0:tj+1,::qix], marker='.', cmap='pyart_HomeyerRainbow', vmin=-0.4, vmax=0.1)
        ax.set_xlabel('E-W distance (km)')
        ax.set_ylabel('N-S distance (km)')
        ax.set_xlim(xl)
        ax.set_ylim(yl)
        ax.set_title(f"{box_name} - Trajectories by buoyancy at {stime2/60:.0f} min \n {cond_str}")
        plt.show()
        if figsave:
            plt.savefig(ip+f"imgs_trajectories/{img_str}_traj-buoy_{stime2/60:.0f}min.png", dpi=300)
        

#%% Kmeans clustering

# from sklearn.cluster import KMeans

figsave = False

wcss1 = np.zeros(shape=(10,), dtype=float)
wcss2 = np.zeros(shape=(10,), dtype=float)
wcss3 = np.zeros(shape=(10,), dtype=float)
for i in range(1,11):
    km1 = KMeans(n_clusters=i, random_state=0, n_init='auto').fit(np.array([b_mv[0,:], z_mv[0,:]]).transpose())
    wcss1[i-1] = km1.inertia_
    km2 = KMeans(n_clusters=i, random_state=0, n_init='auto').fit(np.array([b_mv[0,:], z_mv[0,:], x_mv[0,:], y_mv[0,:]]).transpose())
    wcss2[i-1] = km2.inertia_
    # km1 = KMeans(n_clusters=i, random_state=0, n_init='auto').fit(np.array([np.mean(b_mv,axis=0), z_mv[0,:]]).transpose())
    # wcss1[i-1] = km1.inertia_
    # km2 = KMeans(n_clusters=i, random_state=0, n_init='auto').fit(np.array([np.mean(b_mv,axis=0), z_mv[0,:], x_mv[0,:], y_mv[0,:]]).transpose())
    # wcss2[i-1] = km2.inertia_
    
fig,((ax1),(ax2)) = plt.subplots(2,1,figsize=(9,6))
ax1.plot(range(1,11), wcss1)
ax1.set_title('starting B+z')
ax2.plot(range(1,11), wcss2)
ax2.set_title('starting B+z+x+y')
plt.suptitle(f"{box_name} - KMeans clustering WCSS \n {cond_str}")
plt.show()
if figsave:
    plt.savefig(ip+f"imgs_trajectories/kmeans_WCSS_{img_str}_{stime/60:.0f}.png", dpi=400)
    

n = 2
km1 = KMeans(n_clusters=n, random_state=0, n_init='auto').fit(np.array([b_mv[0,:], z_mv[0,:]]).transpose())
km2 = KMeans(n_clusters=n, random_state=0, n_init='auto').fit(np.array([b_mv[0,:], z_mv[0,:], x_mv[0,:], y_mv[0,:]]).transpose())


b1 = km1.cluster_centers_[:,0]
z1 = km1.cluster_centers_[:,1]
b2 = km2.cluster_centers_[:,0]
z2 = km2.cluster_centers_[:,1]
x2 = km2.cluster_centers_[:,2]
y2 = km2.cluster_centers_[:,3]


# Cluster centroids with trajectories colored by buoyancy, z
fig,(ax1,ax2) = plt.subplots(1,2,figsize=(14,5))

p = ax1.scatter(x_mv[:,:]/1000, y_mv[:,:]/1000, s=10, c=b_mv[:,:], marker='.', cmap='pyart_HomeyerRainbow', vmin=-0.4, vmax=0.1)
ax1.contour(xh, yh, dbz, levels=[30], colors='k', linewidths=1)
ax1.scatter(x2/1000, y2/1000, s=400, c='k', marker='.')
ax1.set_xlim(xl)
ax1.set_ylim(yl)
ax1.set_title(f"Buoyancy, 30 dBZ")
plt.colorbar(p, ax=ax1)
ax1.set_xlabel('E-W distance (km)')
ax1.set_ylabel('N-S distance (km)')

p = ax2.scatter(x_mv[:,:]/1000, y_mv[:,:]/1000, s=10, c=z_mv[0,:], marker='.', cmap='pyart_HomeyerRainbow', vmin=-0.4, vmax=0.1)
ax2.contour(xh, yh, dbz, levels=[30], colors='k', linewidths=1)
ax2.scatter(x2/1000, y2/1000, s=400, c='k', marker='.')
ax2.set_xlim(xl)
ax2.set_ylim(yl)
ax2.set_title(f"Height AGL, 30 dBZ")
plt.colorbar(p, ax=ax2)
ax2.set_xlabel('E-W distance (km)')
ax2.set_ylabel('N-S distance (km)')

plt.suptitle(f"{box_name} - {n} K-means cluster centroids \n {cond_str}")
plt.show()
if figsave:
    plt.savefig(ip+f"imgs_trajectories/kmeans-{n}clustercentroids_{img_str}_{stime/60:.0f}.png", dpi=400)


# Trajectories colored by cluster
fig,(ax1,ax2) = plt.subplots(1,2,figsize=(12,5))
ax1.scatter(x_mv[:,:]/1000, y_mv[:,:]/1000, s=10, c=np.tile(km1.labels_, (len(x_mv[:,0]),1)), marker='.', cmap='pyart_HomeyerRainbow')
ax1.contour(xh, yh, dbz, levels=[30], colors='k', linewidths=1)
ax1.set_xlim(xl)
ax1.set_ylim(yl)
ax1.set_title(f"Clusters using B+z, 30 dBZ")
ax2.scatter(x_mv[:,:]/1000, y_mv[:,:]/1000, s=10, c=np.tile(km2.labels_, (len(x_mv[:,0]),1)), marker='.', cmap='pyart_HomeyerRainbow')
ax2.contour(xh, yh, dbz, levels=[30], colors='k', linewidths=1)
ax2.set_xlim(xl)
ax2.set_ylim(yl)
ax2.set_title(f"Clusters using B+z+x+y, 30 dBZ")
plt.suptitle(f"{box_name} - {n} K-means trajectory clusters - {cond_str}")
plt.show()
if figsave:
    plt.savefig(ip+f"imgs_trajectories/kmeans-{n}clusters-full_{img_str}_{stime/60:.0f}.png", dpi=400)

# Trajectory starting points colored by cluster
fig,(ax1,ax2) = plt.subplots(1,2,figsize=(12,5))
ax1.scatter(x_mv[0,:]/1000, y_mv[0,:]/1000, s=10, c=km1.labels_, marker='.', cmap='pyart_HomeyerRainbow')
ax1.contour(xh, yh, dbz0, levels=[30], colors='k', linewidths=1)
ax1.set_xlim(xl)
ax1.set_ylim(yl)
ax1.set_title(f"Clusters using B+z, initial 30 dBZ")
ax2.scatter(x_mv[0,:]/1000, y_mv[0,:]/1000, s=10, c=km2.labels_, marker='.', cmap='pyart_HomeyerRainbow')
ax2.contour(xh, yh, dbz0, levels=[30], colors='k', linewidths=1)
ax2.set_xlim(xl)
ax2.set_ylim(yl)
ax2.set_title(f"Clusters using B+z+x+y, initial 30 dBZ")
plt.suptitle(f"{box_name} - {n} K-means trajectory cluster starting points - {cond_str}")
plt.show()
if figsave:
    plt.savefig(ip+f"imgs_trajectories/kmeans-{n}clusters-start_{img_str}_{stime/60:.0f}.png", dpi=400)

# Conclusion: B+z sucks. x+y and B+z+x+y are the same.
# MV1, w>1 - mid-level source is only distinct with b+z, but x+y looks overall better. maybe 4 or even 5 is the way to go here?
# MV1, w>5 - could be either 2 or 3 clusters. rear source not a distinct cluster, but mid-level source is and both forward sources are.
# MV2, w>1 - 3 clusters looks best
# MV2, w>5 - 3 clusters looks a little better than 2, but neither is great - the two forward sources aren't really distinct

#%% Gaussian mixture model clustering

# from sklearn.mixture import GaussianMixture
# from sklearn.metrics import silhouette_score

# figsave = False

aic1 = np.zeros(shape=(10,),dtype=float)
aic2 = np.zeros(shape=(10,),dtype=float)
bic1 = np.zeros(shape=(10,),dtype=float)
bic2 = np.zeros(shape=(10,),dtype=float)
score1 = np.zeros(shape=(10,),dtype=float)
score2 = np.zeros(shape=(10,),dtype=float)
sil1 = np.zeros(shape=(10,),dtype=float)
sil2 = np.zeros(shape=(10,),dtype=float)
for i in range(2,12):
    X1 = np.array([b_mv[0,:], z_mv[0,:]]).transpose()
    X2 = np.array([b_mv[0,:], z_mv[0,:], x_mv[0,:], y_mv[0,:]]).transpose()
    gm1 = GaussianMixture(n_components=i, init_params='k-means++', max_iter=1000).fit(X1)
    gm2 = GaussianMixture(n_components=i, init_params='k-means++', max_iter=1000).fit(X2)
    # X1 = np.array([np.mean(b_mv,axis=0), z_mv[0,:]]).transpose()
    # X2 = np.array([np.mean(b_mv,axis=0), z_mv[0,:], x_mv[0,:], y_mv[0,:]]).transpose()
    # gm1 = GaussianMixture(n_components=i, init_params='k-means++').fit(X1)
    # gm2 = GaussianMixture(n_components=i, init_params='k-means++').fit(X2)
    aic1[i-2] = gm1.aic(X1)
    aic2[i-2] = gm2.aic(X2)
    bic1[i-2] = gm1.bic(X1)
    bic2[i-2] = gm2.bic(X2)
    score1[i-2] = gm1.score(X1)
    score2[i-2] = gm2.score(X2)
    
    sil1[i-2] = silhouette_score(X1, gm1.predict(X1))
    sil2[i-2] = silhouette_score(X2, gm2.predict(X2))
    
    # labels1 = gm1.predict(X1)
    # labels2 = gm2.predict(X2)
    

fig,((ax1),(ax2)) = plt.subplots(2,1,figsize=(9,6))
ax1.plot(range(2,12), sil1)
ax1.set_title('starting B+z')
ax2.plot(range(2,12), sil2)
ax2.set_title('starting B+z+x+y')
plt.suptitle(f"{box_name} - GMM clustering silhouette score \n {cond_str}")
plt.show()
if figsave:
    plt.savefig(ip+f"imgs_trajectories/GMM-silscore_{img_str}_{stime/60:.0f}.png", dpi=400)
    
    

# for i in range(2,6):
#     X1 = np.array([b_mv[0,:], z_mv[0,:]]).transpose()
#     X2 = np.array([b_mv[0,:], z_mv[0,:], x_mv[0,:], y_mv[0,:]]).transpose()
#     gm1 = GaussianMixture(n_components=i, init_params='k-means++').fit(X1)
#     gm2 = GaussianMixture(n_components=i, init_params='k-means++').fit(X2)
#     # X1 = np.array([np.mean(b_mv,axis=0), z_mv[0,:]]).transpose()
#     # X2 = np.array([np.mean(b_mv,axis=0), z_mv[0,:], x_mv[0,:], y_mv[0,:]]).transpose()
#     # gm1 = GaussianMixture(n_components=i, init_params='k-means++')).fit(X1)
#     # gm2 = GaussianMixture(n_components=i, init_params='k-means++')).fit(X2))
#     labels1 = gm1.predict(X1)
#     labels2 = gm2.predict(X2)
#     sil1 = silhouette_score(X1, gm1.predict(X1))
#     sil2 = silhouette_score(X2, gm2.predict(X2))
    
    
#     # Trajectories colored by cluster
#     fig,(ax1,ax2) = plt.subplots(1,2,figsize=(12,5))
#     ax1.scatter(x_mv[:,:]/1000, y_mv[:,:]/1000, s=1, c=np.tile(labels1, (len(x_mv[:,0]),1)), marker='.', cmap='pyart_HomeyerRainbow')
#     ax1.contour(xh, yh, dbz, levels=[30], colors='k', lineiwdths=1)
#     ax1.set_xlim(xl)
#     ax1.set_ylim(yl)
#     ax1.set_title(f"Clusters using B+z, 30 dBZ (SIL = {sil1:.2f})")
#     ax2.scatter(x_mv[:,:]/1000, y_mv[:,:]/1000, s=1, c=np.tile(labels2, (len(x_mv[:,0]),1)), marker='.', cmap='pyart_HomeyerRainbow')
#     ax2.contour(xh, yh, dbz, levels=[30], colors='k', lineiwdths=1)
#     ax2.set_xlim(xl)
#     ax2.set_ylim(yl)
#     ax2.set_title(f"Clusters using B+z+x+y, 30 dBZ (SIL = {sil2:.2f})")
#     plt.suptitle(f"{box_name} - {i} GMM trajectory cluster starting points - {cond_str}")
    
#     # Trajectory starting points colored by cluster
#     fig,(ax1,ax2) = plt.subplots(1,2,figsize=(12,5))
#     ax1.scatter(x_mv[0,:]/1000, y_mv[0,:]/1000, s=1, c=labels1, marker='.', cmap='pyart_HomeyerRainbow')
#     ax1.contour(xh, yh, dbz0, levels=[30], colors='k', lineiwdths=1)
#     ax1.set_xlim(xl)
#     ax1.set_ylim(yl)
#     ax1.set_title(f"Clusters using B+z, initial 30 dBZ")
#     ax2.scatter(x_mv[0,:]/1000, y_mv[0,:]/1000, s=1, c=labels2, marker='.', cmap='pyart_HomeyerRainbow')
#     ax2.contour(xh, yh, dbz0, levels=[30], colors='k', lineiwdths=1)
#     ax2.set_xlim(xl)
#     ax2.set_ylim(yl)
#     ax2.set_title(f"Clusters using B+z+x+y, initial 30 dBZ")
#     plt.suptitle(f"{box_name} - {i} GMM trajectory cluster starting points - {cond_str}")


#%% Actually find GMM clusters

n = 3
X1 = np.array([b_mv[0,:], z_mv[0,:], w_mv[0,:], x_mv[0,:], y_mv[0,:]]).transpose()
# X2 = np.array([b_mv[0,:], z_mv[0,:], x_mv[0,:], y_mv[0,:]]).transpose()
gm1 = GaussianMixture(n_components=n, init_params='k-means++', max_iter=1000).fit(X1)
# gm2 = GaussianMixture(n_components=n, init_params='k-means++', max_iter=1000).fit(X2)
# X1 = np.array([np.mean(b_mv,axis=0), z_mv[0,:]]).transpose()
# X2 = np.array([np.mean(b_mv,axis=0), z_mv[0,:], x_mv[0,:], y_mv[0,:]]).transpose()
# gm1 = GaussianMixture(n_components=i, init_params='k-means++')).fit(X1)
# gm2 = GaussianMixture(n_components=i, init_params='k-means++')).fit(X2))
labels1 = gm1.predict(X1)
# labels2 = gm2.predict(X2)
sil1 = silhouette_score(X1, gm1.predict(X1))
# sil2 = silhouette_score(X2, gm2.predict(X2))

# # Trajectories colored by cluster
# fig,(ax1,ax2) = plt.subplots(1,2,figsize=(12,5))
# ax1.scatter(x_mv[:,:]/1000, y_mv[:,:]/1000, s=1, c=np.tile(labels1, (len(x_mv[:,0]),1)), marker='.', cmap='pyart_HomeyerRainbow')
# ax1.contour(xh, yh, dbz, levels=[30], colors='k', linewidths=1)
# ax1.set_xlim(xl)
# ax1.set_ylim(yl)
# ax1.set_title(f"Clusters using B+z, 30 dBZ (SIL = {sil1:.2f})")
# ax2.scatter(x_mv[:,:]/1000, y_mv[:,:]/1000, s=1, c=np.tile(labels2, (len(x_mv[:,0]),1)), marker='.', cmap='pyart_HomeyerRainbow')
# ax2.contour(xh, yh, dbz, levels=[30], colors='k', linewidths=1)
# ax2.set_xlim(xl)
# ax2.set_ylim(yl)
# ax2.set_title(f"Clusters using B+z+x+y, 30 dBZ (SIL = {sil2:.2f})")
# plt.suptitle(f"{box_name} - {n} GMM trajectory cluster starting points - {cond_str}")
# plt.show()

# # Trajectory starting points colored by cluster
# fig,(ax1,ax2) = plt.subplots(1,2,figsize=(12,5))
# ax1.scatter(x_mv[0,:]/1000, y_mv[0,:]/1000, s=1, c=labels1, marker='.', cmap='pyart_HomeyerRainbow')
# ax1.contour(xh, yh, dbz0, levels=[30], colors='k', linewidths=1)
# ax1.set_xlim(xl)
# ax1.set_ylim(yl)
# ax1.set_title(f"Clusters using B+z, initial 30 dBZ")
# p = ax2.scatter(x_mv[0,:]/1000, y_mv[0,:]/1000, s=1, c=labels2, marker='.', cmap='pyart_HomeyerRainbow')
# ax2.contour(xh, yh, dbz0, levels=[30], colors='k', linewidths=1)
# ax2.set_xlim(xl)
# ax2.set_ylim(yl)
# ax2.set_title(f"Clusters using B+z+x+y, initial 30 dBZ")
# plt.colorbar(p, ax=ax2)
# plt.suptitle(f"{box_name} - {n} GMM trajectory cluster starting points - {cond_str}")
# plt.show()


fig,(ax1,ax2) = plt.subplots(1,2,figsize=(12,5))
ax1.scatter(x_mv[:,:]/1000, y_mv[:,:]/1000, s=1, c=np.tile(labels1, (len(x_mv[:,0]),1)), marker='.', cmap='pyart_HomeyerRainbow')
ax1.contour(xh, yh, dbz, levels=[30], colors='k', linewidths=1)
ax1.set_xlim(xl)
ax1.set_ylim(yl)
ax1.set_title(f"Clusters, 30 dBZ (SIL = {sil1:.2f})")
p = ax2.scatter(x_mv[0,:]/1000, y_mv[0,:]/1000, s=1, c=labels1, marker='.', cmap='pyart_HomeyerRainbow')
ax2.contour(xh, yh, dbz0, levels=[30], colors='k', linewidths=1)
ax2.set_xlim(xl)
ax2.set_ylim(yl)
ax2.set_title(f"Cluster starting points, initial 30 dBZ")
plt.colorbar(p, ax=ax2)
plt.suptitle(f"{box_name} - {n} GMM trajectory clusters - {cond_str}")
plt.show()


#%% Time series

cc = cc
c0 = (cc == 0)
c1 = (cc == 1)
c2 = (cc == 2)
c3 = (cc == 3)
c4 = (cc == 4)



figsave = 0

#% Plot trajectory z, w, zvort time series by source
# mean trajectory is not ideal but neither is max or median - 90th percentile?
if False:
    tl = [stime/60-15,stime/60]
    
    fig,((ax1),(ax2),(ax3),(ax4)) = plt.subplots(4,1,figsize=(12,15))
    
    ax1.plot(ptime[0:ti+1]/60, z_mv[:,c0]/1000, '-b', linewidth=0.25)
    ax1.plot(ptime[0:ti+1]/60, z_mv[:,c1]/1000, '-g', linewidth=0.25)
    ax1.plot(ptime[0:ti+1]/60, z_mv[:,c2]/1000, 'orange', linewidth=0.25)
    ax1.plot(ptime[0:ti+1]/60, z_mv[:,c3]/1000, '-r', linewidth=0.25)
    ax1.plot(ptime[0:ti+1]/60, z_mv[:,c4]/1000, '-k', linewidth=0.25)
    ax1.set_ylabel('z (km)')
    ax1.set_title(f"{box_name} trajectory height")
    # ax1.legend(labels=['LL env', 'CP', 'ML env', 'QLCS'])
    ax1.set_xlim(tl)
    
    ax2.plot(ptime[0:ti+1]/60, w_mv[:,c0], '-b', linewidth=0.25)
    ax2.plot(ptime[0:ti+1]/60, w_mv[:,c1], '-g', linewidth=0.25)
    ax2.plot(ptime[0:ti+1]/60, w_mv[:,c2], 'orange', linewidth=0.25)
    ax2.plot(ptime[0:ti+1]/60, w_mv[:,c3], '-r', linewidth=0.25)
    ax2.plot(ptime[0:ti+1]/60, w_mv[:,c4], '-k', linewidth=0.25)
    ax2.set_ylabel('w (m/s)')
    ax2.set_title(f"{box_name} trajectory $w$")
    ax2.set_xlim(tl)
    
    ax3.plot(ptime[0:ti+1]/60, zvort_mv[:,c0], '-b', linewidth=0.25)
    ax3.plot(ptime[0:ti+1]/60, zvort_mv[:,c1], '-g', linewidth=0.25)
    ax3.plot(ptime[0:ti+1]/60, zvort_mv[:,c2], 'orange', linewidth=0.25)
    ax3.plot(ptime[0:ti+1]/60, zvort_mv[:,c3], '-r', linewidth=0.25)
    ax3.plot(ptime[0:ti+1]/60, zvort_mv[:,c4], '-k', linewidth=0.25)
    ax3.set_ylabel('zvort (1/s)')
    ax3.set_title(f"{box_name} trajectory zvort")
    ax3.set_xlim(tl)
    
    ax4.plot(ptime[0:ti+1]/60, b_mv[:,c0], '-b', linewidth=0.25)
    ax4.plot(ptime[0:ti+1]/60, b_mv[:,c1], '-g', linewidth=0.25)
    ax4.plot(ptime[0:ti+1]/60, b_mv[:,c2], 'orange', linewidth=0.25)
    ax4.plot(ptime[0:ti+1]/60, b_mv[:,c3], '-r', linewidth=0.25)
    ax4.plot(ptime[0:ti+1]/60, b_mv[:,c4], '-k', linewidth=0.25)
    ax4.set_ylabel('buoyancy')
    ax4.set_xlabel('Time (min)')
    ax4.set_xlim(tl)
    ax4.set_title(f"{box_name} trajectory buoyancy")
    
    plt.show()
    
if True:
    tl = [stime/60-15,stime/60]
    
    fig,((ax1),(ax2),(ax3),(ax4)) = plt.subplots(4,1,figsize=(12,15))
    
    ax1.plot(ptime[0:ti+1]/60, np.median(z_mv[:,c0],axis=1)/1000, '-b', linewidth=1)
    ax1.plot(ptime[0:ti+1]/60, np.median(z_mv[:,c1],axis=1)/1000, '-g', linewidth=1)
    ax1.plot(ptime[0:ti+1]/60, np.median(z_mv[:,c2],axis=1)/1000, 'orange', linewidth=1)
    ax1.plot(ptime[0:ti+1]/60, np.median(z_mv[:,c3],axis=1)/1000, '-r', linewidth=1)
    ax1.plot(ptime[0:ti+1]/60, np.median(z_mv[:,c4],axis=1)/1000, '-k', linewidth=1)
    if np.any(c0):
        ax1.fill_between(ptime[0:ti+1]/60, np.percentile(z_mv[:,c0],25,axis=1)/1000, np.percentile(z_mv[:,c0],75,axis=1)/1000,
                         alpha=0.2, color='b')
    if np.any(c1):
        ax1.fill_between(ptime[0:ti+1]/60, np.percentile(z_mv[:,c1],25,axis=1)/1000, np.percentile(z_mv[:,c1],75,axis=1)/1000,
                         alpha=0.2, color='g')
    if np.any(c2):
        ax1.fill_between(ptime[0:ti+1]/60, np.percentile(z_mv[:,c2],25,axis=1)/1000, np.percentile(z_mv[:,c2],75,axis=1)/1000,
                          alpha=0.2, color='orange')
    if np.any(c3):
        ax1.fill_between(ptime[0:ti+1]/60, np.percentile(z_mv[:,c3],25,axis=1)/1000, np.percentile(z_mv[:,c3],75,axis=1)/1000,
                          alpha=0.2, color='r')
    if np.any(c4):
        ax1.fill_between(ptime[0:ti+1]/60, np.percentile(z_mv[:,c4],25,axis=1)/1000, np.percentile(z_mv[:,c4],75,axis=1)/1000,
                          alpha=0.1, color='k')
    ax1.set_ylabel('z (km)')
    ax1.legend(labels=['Lower env.', 'Cold pool', 'Upper env.', 'QLCS', 'Upper QLCS'], loc=2)
    # ax1.legend(labels=['Lower env.', 'Upper env.', 'Cold pool'], loc=2)
    ax1.set_xlim(tl)
    ax1.set_ylim([0,3])
    ax1.set_title(f"{box_name} median trajectory height")
    
    
    ax2.plot(ptime[0:ti+1]/60, np.median(w_mv[:,c0],axis=1), '-b', linewidth=1)
    ax2.plot(ptime[0:ti+1]/60, np.median(w_mv[:,c1],axis=1), '-g', linewidth=1)
    ax2.plot(ptime[0:ti+1]/60, np.median(w_mv[:,c2],axis=1), 'orange', linewidth=1)
    ax2.plot(ptime[0:ti+1]/60, np.median(w_mv[:,c3],axis=1), '-r', linewidth=1)
    ax2.plot(ptime[0:ti+1]/60, np.median(w_mv[:,c4],axis=1), '-k', linewidth=1)
    if np.any(c0):
        ax2.fill_between(ptime[0:ti+1]/60, np.percentile(w_mv[:,c0],25,axis=1), np.percentile(w_mv[:,c0],75,axis=1),
                         alpha=0.2, color='b')
    if np.any(c1):
        ax2.fill_between(ptime[0:ti+1]/60, np.percentile(w_mv[:,c1],25,axis=1), np.percentile(w_mv[:,c1],75,axis=1),
                         alpha=0.2, color='g')
    if np.any(c2):
        ax2.fill_between(ptime[0:ti+1]/60, np.percentile(w_mv[:,c2],25,axis=1), np.percentile(w_mv[:,c2],75,axis=1),
                          alpha=0.2, color='orange')
    if np.any(c3):
        ax2.fill_between(ptime[0:ti+1]/60, np.percentile(w_mv[:,c3],25,axis=1), np.percentile(w_mv[:,c3],75,axis=1),
                          alpha=0.2, color='r')
    if np.any(c4):
        ax2.fill_between(ptime[0:ti+1]/60, np.percentile(w_mv[:,c4],25,axis=1), np.percentile(w_mv[:,c4],75,axis=1),
                          alpha=0.1, color='k')
    ax2.axhline(y=0, color='k', linestyle='--', linewidth=1)
    ax2.set_ylabel('w (m/s)')
    ax2.set_xlim(tl)
    ax2.set_ylim([-15,15])
    ax2.set_title(f"{box_name} median trajectory $w$")
    
    
    ax3.plot(ptime[0:ti+1]/60, np.median(zvort_mv[:,c0],axis=1), '-b', linewidth=1)
    ax3.plot(ptime[0:ti+1]/60, np.median(zvort_mv[:,c1],axis=1), '-g', linewidth=1)
    ax3.plot(ptime[0:ti+1]/60, np.median(zvort_mv[:,c2],axis=1), 'orange', linewidth=1)
    ax3.plot(ptime[0:ti+1]/60, np.median(zvort_mv[:,c3],axis=1), '-r', linewidth=1)
    ax3.plot(ptime[0:ti+1]/60, np.median(zvort_mv[:,c4],axis=1), '-k', linewidth=1)
    if np.any(c0):
        ax3.fill_between(ptime[0:ti+1]/60, np.percentile(zvort_mv[:,c0],25,axis=1), np.percentile(zvort_mv[:,c0],75,axis=1),
                         alpha=0.2, color='b')
    if np.any(c1):
        ax3.fill_between(ptime[0:ti+1]/60, np.percentile(zvort_mv[:,c1],25,axis=1), np.percentile(zvort_mv[:,c1],75,axis=1),
                         alpha=0.2, color='g')
    if np.any(c2):
        ax3.fill_between(ptime[0:ti+1]/60, np.percentile(zvort_mv[:,c2],25,axis=1), np.percentile(zvort_mv[:,c2],75,axis=1),
                          alpha=0.2, color='orange')
    if np.any(c3):
        ax3.fill_between(ptime[0:ti+1]/60, np.percentile(zvort_mv[:,c3],25,axis=1), np.percentile(zvort_mv[:,c3],75,axis=1),
                          alpha=0.2, color='r')
    if np.any(c4):
        ax3.fill_between(ptime[0:ti+1]/60, np.percentile(zvort_mv[:,c4],25,axis=1), np.percentile(zvort_mv[:,c4],75,axis=1),
                          alpha=0.1, color='k')
    ax3.axhline(y=0, color='k', linestyle='--', linewidth=1)
    ax3.set_ylabel('zvort (1/s)')
    ax3.set_xlim(tl)
    ax3.set_ylim([-0.03,0.04])
    ax3.set_title(f"{box_name} median trajectory \u03B6")
    
    
    ax4.plot(ptime[0:ti+1]/60, np.median(b_mv[:,c0],axis=1), '-b', linewidth=1)
    ax4.plot(ptime[0:ti+1]/60, np.median(b_mv[:,c1],axis=1), '-g', linewidth=1)
    ax4.plot(ptime[0:ti+1]/60, np.median(b_mv[:,c2],axis=1), 'orange', linewidth=1)
    ax4.plot(ptime[0:ti+1]/60, np.median(b_mv[:,c3],axis=1), '-r', linewidth=1)
    ax4.plot(ptime[0:ti+1]/60, np.median(b_mv[:,c4],axis=1), '-k', linewidth=1)
    if np.any(c0):
        ax4.fill_between(ptime[0:ti+1]/60, np.percentile(b_mv[:,c0],25,axis=1), np.percentile(b_mv[:,c0],75,axis=1),
                         alpha=0.2, color='b')
    if np.any(c1):
        ax4.fill_between(ptime[0:ti+1]/60, np.percentile(b_mv[:,c1],25,axis=1), np.percentile(b_mv[:,c1],75,axis=1),
                         alpha=0.2, color='g')
    if np.any(c2):
        ax4.fill_between(ptime[0:ti+1]/60, np.percentile(b_mv[:,c2],25,axis=1), np.percentile(b_mv[:,c2],75,axis=1),
                          alpha=0.2, color='orange')
    if np.any(c3):
        ax4.fill_between(ptime[0:ti+1]/60, np.percentile(b_mv[:,c3],25,axis=1), np.percentile(b_mv[:,c3],75,axis=1),
                          alpha=0.2, color='r')
    if np.any(c4):
        ax4.fill_between(ptime[0:ti+1]/60, np.percentile(b_mv[:,c4],25,axis=1), np.percentile(b_mv[:,c4],75,axis=1),
                          alpha=0.1, color='k')
    ax4.axhline(y=0, color='k', linestyle='--', linewidth=1)
    ax4.set_ylabel('buoyancy')
    ax4.set_xlim(tl)
    ax4.set_xlabel('Time (min)')
    ax4.set_ylim([-0.4,0.1])
    ax4.set_title(f"{box_name} median trajectory buoyancy")
    
    if figsave:
        plt.savefig(ip+f"imgs_trajectories/{img_str}_timeseries-median_{stime/60:.0f}min.png", dpi=300)
    
    plt.show()





#%% Find parcels in the area of interest and interpolate pgfz (or whatever variable) to parcel locations

from scipy.interpolate import interpn,interp2d,RegularGridInterpolator

fp = '/Volumes/Promise_Pegasus_70TB/merger/merger-125m'

ds = nc.Dataset(fp+'/cm1out_pdata2.nc')
ptime = ds.variables['time'][:].data
pid = ds.variables['xh'][:].data
x = ds.variables['x'][:].data/1000
y = ds.variables['y'][:].data/1000
z = ds.variables['z'][:].data/1000
# u = ds.variables['u'][:].data
# v = ds.variables['v'][:].data
w = ds.variables['w'][:].data
# dbz = ds.variables['dbz'][:].data
zvort = ds.variables['zvort'][:].data
# wspd = (ds.variables['u'][:].data**2 + ds.variables['v'][:].data**2)**0.5
ds.close()

t_m = 210 # time of merger
i = np.where(ptime==t_m*60)[0][0]

x1,x2 = -20,-16
y1,y2 = -75,-71

pids = pid[(x[i,:]>=x1) & (x[i,:]<=x2) & (y[i,:]>=y1) & (y[i,:]<=y2)]
xp = x[i, (x[i,:]>=x1) & (x[i,:]<=x2) & (y[i,:]>=y1) & (y[i,:]<=y2)]
yp = y[i, (x[i,:]>=x1) & (x[i,:]<=x2) & (y[i,:]>=y1) & (y[i,:]<=y2)]
zp = z[i, (x[i,:]>=x1) & (x[i,:]<=x2) & (y[i,:]>=y1) & (y[i,:]<=y2)]
# xp = np.array([x[:,int(iprcl)] for iprcl in pids]).transpose()
# yp = np.array([y[:,int(iprcl)] for iprcl in pids]).transpose()
# zp = np.array([z[:,int(iprcl)] for iprcl in pids]).transpose()

if False:
    stimes = np.linspace(10860,14400,60)
    pgfz_dn = np.zeros(shape=(len(stimes),len(pids)), dtype=float)
    for fnum in np.arange(14,74):
        print(f"Interpolating file {fnum}...")
        ds = nc.Dataset(fp+f"/pp/dyn_{fnum:06d}.nc")
        stime = ds.variables['time'][:].data
        i = np.where(ptime==stime)[0][0]
        xh = ds.variables['xh'][:].data
        yh = ds.variables['yh'][:].data
        zf = ds.variables['zf'][:].data
        
        gridpts = (zf, yh, xh)
        plocs = (zp[i,:], yp[i,:], xp[i,:])
        pgfz_dn[fnum-14,:] = interpn(gridpts, ds.variables['pgfz_dn'][:].data, plocs, method='slinear')
        ds.close()
    
    
    f = RegularGridInterpolator((stimes,pids), pgfz_dn, method='slinear', bounds_error=False, fill_value=None)
    vpgdn_p = f((ptime, pids))

#%% 3D plots

xlims = [-40,10]
ylims = [-110,-60]
zlims = [0,10]
vlim = 5

iy1 = np.where(yh >= ylims[0])[0][0]
iy2 = np.where(yh >= ylims[1])[0][0]
ix1 = np.where(xh >= xlims[0])[0][0]
ix2 = np.where(xh >= xlims[1])[0][0]

xx,yy = np.meshgrid(xh[ix1:ix2], yh[iy1:iy2], indexing='xy')


fig = plt.figure(figsize=(9,6))

ax1 = fig.add_subplot(1,1,1, projection='3d')
ax1.contour(xx, yy, dbz[iy1:iy2,ix1:ix2], levels=[30], zdir='z', offset=0, colors='k', linewidths=1)
ax1.scatter(x_mv[:,::qix]/1000, y_mv[:,::qix]/1000, z_mv[:,::qix]/1000, s=1, c=b_mv[:,::qix], marker='.', cmap='pyart_HomeyerRainbow', vmin=-0.3, vmax=0.1)
ax1.set_xlim(xlims)
ax1.set_ylim(ylims)
ax1.set_zlim(zlims)
ax1.set_xlabel('x (km)', fontsize=14)
ax1.set_ylabel('y (km)', fontsize=14)
ax1.set_zlabel('z (km)', fontsize=14)
ax1.set_title(f"B, all sources", fontsize=14)

plt.tight_layout()
plt.show()


#%% Fix repeated times in merger pdata -- qlcs and supercell are fine --> done!

if False:
    from CM1utils import *
    import xarray as xr
    
    ds = nc.Dataset('/Volumes/Promise_Pegasus_70TB/merger/merger-125m/cm1out_pdata.nc')
    xh = ds.variables['xh'][:].data
    yh = ds.variables['yh'][:].data
    zh = ds.variables['zh'][:].data
    time = np.delete(ds.variables['time'][:].data, [201,202])
    
    dvars = [v for v in list(ds.variables.keys()) if (v not in ['xh','yh','zh','time'])]
    data_vars = dict()
    for v in dvars:
        print(f"Fixing {v}")
        tmp = np.delete(ds.variables[v][:].data, [201,202], axis=0)
        data_vars[v] = xr.DataArray(data=tmp, coords=dict(time=time, xh=xh), dims=['time','xh'], attrs=dict(long_name=ds.variables[v].long_name, units=ds.variables[v].units))
    
    xh = xr.DataArray(data=xh, coords=dict(xh=xh), dims=['xh'], attrs=dict(long_name='parcel ID number', units=''))
    yh = xr.DataArray(data=yh, coords=dict(yh=yh), dims=['yh'], attrs=dict(long_name='', units=''))
    zh = xr.DataArray(data=zh, coords=dict(zh=zh), dims=['zh'], attrs=dict(long_name='', units=''))
    time = xr.DataArray(data=time, coords=dict(time=time), dims=['time'], attrs=dict(long_name='time', units='s'))
    
    coords = {'xh':xh, 'yh':yh, 'zh':zh, 'time':time}
    
    dset = xr.Dataset(data_vars=data_vars, coords=coords, attrs=dict(cm1_version='cm1r19.1', Conventions='COARDS'))
    dset.to_netcdf('/Volumes/Promise_Pegasus_70TB/merger/merger-125m/cm1out_pdata3.nc')
    
    dset.close()
    ds.close()
    
    del data_vars,coords,tmp



#%% Interpolate missing 180.25 s in merger/qlcs/supercell pdata

if False:
    from CM1utils import *
    import xarray as xr
    from scipy.interpolate import RegularGridInterpolator
    
    sims = ['qlcs', 'supercell']
    # sims = ['merger']
    
    for s in sims:
        print(f"... {s} ...")
        ds = nc.Dataset(f"/Volumes/Promise_Pegasus_70TB/merger/{s}-125m/cm1out_pdata.nc")
        ds2 = nc.Dataset(f"/Volumes/Promise_Pegasus_70TB/merger/{s}-125m/cm1out_pdata2.nc", 'r+', clobber=True)
        time = ds.variables['time'][:].data
        pid = ds.variables['xh'][:].data
        yh = ds.variables['yh'][:].data
        zh = ds.variables['zh'][:].data
        newtime = np.linspace(10800, 14400, 241)
        
        ds2.createDimension('xh', size=len(pid))
        ds2.createDimension('yh', size=len(yh))
        ds2.createDimension('zh', size=len(zh))
        ds2.createDimension('time', size=len(newtime))
        
        cx = ds2.createVariable('xh', 'f4', ('xh'))
        cx.units = ''; cx.long_name = 'parcel ID number'; cx[:] = pid[:]
        cy = ds2.createVariable('yh', 'f4', ('yh'))
        cy.units = ''; cy.long_name = ''; cy[:] = yh[:]
        cz = ds2.createVariable('zh', 'f4', ('zh'))
        cz.units = ''; cz.long_name = ''; cz[:] = zh[:]
        ct = ds2.createVariable('time', 'f4', ('time'))
        ct.units = 's'; cy.long_name = 'time'; ct[:] = newtime[:]
        
        dvars = [v for v in list(ds.variables.keys()) if (v not in ['xh','yh','zh','time'])]
        for v in dvars:
            vint = np.zeros(shape=(len(newtime),len(pid)), dtype=float)
            vint[0,:] = ds.variables[v][:].data[0,:]
            vint[2::,:] = ds.variables[v][:].data[1::,:]
            print(f"Interpolating {v}")
            f = RegularGridInterpolator((time[0:5], pid), ds.variables[v][:].data[0:5,:], method='slinear')
            newdata = f((10815, pid))
            vint[1,:] = newdata
            
            print(f"Writing {v}")
            tmp = ds2.createVariable(v, 'f4', ('time', 'xh'))
            tmp.units = ds.variables[v].units
            tmp.long_name = ds.variables[v].long_name
            tmp[:,:] = vint[:,:]
            
            del tmp,vint
        
        ds2.close()
        ds.close()
    



