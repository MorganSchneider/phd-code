#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 14:34:11 2023

@author: morgan.schneider

- PAPER FIGURES -
traj_z+b+src
traj_z+b+src_Q
traj_z+b+src_S
traj_timeseries
traj_vort_timeseries

- Other figures -
Zoomed-in 5-min time series of mid-level parcel w and zvort
3D trajectory plots
Jiang and Dawson 3D trajectories with vorticity vectors
3D trajectory animations
"""

####################
### Load modules ###
####################

from CM1utils import *
from matplotlib import patches

#%% Load parcel file

# Parcel data dimensions (time, pid)

# fp = '/Volumes/Promise_Pegasus_70TB/merger/merger-125m/'
ip = '/Users/morgan.schneider/Documents/merger/merger-125m/'

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
b = ds.variables['b'][:].data
zvort = ds.variables['zvort'][:].data
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
# vpg = ds.variables['vpg'][:].data
ds.close()


#%% Load data and filter trajectories

# Read gridded output
fnum = 53
ds = nc.Dataset(fp+f"cm1out_{fnum:06d}.nc")
stime = ds.variables['time'][:].data[0]
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
zh = ds.variables['z'][:].data
dbz = ds.variables['dbz'][:].data[0,0,:,:]
ds.close()

ds = nc.Dataset(fp+"base/cm1out_000013.nc")
dbz0 = ds.variables['dbz2'][:].data[0,0,:,:]
ds.close()

# Filter trajectories



ti = np.where(ptime == stime)[0][0]

box_name = 'MERGER'

dbfile = open(ip+'boxes_s1.pkl', 'rb')
box = pickle.load(dbfile)
x1 = box['x1_pp'][fnum-14]
x2 = box['x2_pp'][fnum-14]
y1 = box['y1_pp'][fnum-14]
y2 = box['y2_pp'][fnum-14]
dbfile.close()


wvort_cond = (zvort[ti,:]>=0.01) & (w[ti,:]>=5) & (z[ti,:]<=1000) & (z[ti,:]>=100) & (x[ti,:]>=x1*1000) & (x[ti,:]<=x2*1000) & (y[ti,:]>=y1*1000) & (y[ti,:]<=y2*1000)

pids = pid[wvort_cond]
x_mv = x[:, wvort_cond] # was 0:ti+1
y_mv = y[:, wvort_cond]
z_mv = z[:, wvort_cond]
w_mv = w[:, wvort_cond]
zvort_mv = zvort[:, wvort_cond]
b_mv = b[:, wvort_cond]
# u_mv = u[:, wvort_cond]
# v_mv = v[:, wvort_cond]



#%% Actually find GMM clusters

from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_score

n = 2
X1 = np.array([b_mv[0,:], x_mv[0,:], y_mv[0,:]]).transpose()
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


fig,(ax1,ax2) = plt.subplots(1,2,figsize=(12,5))
ax1.scatter(x_mv[:ti,:]/1000, y_mv[:ti,:]/1000, s=1, c=np.tile(labels1, (len(x_mv[:ti,0]),1)), marker='.', cmap='HomeyerRainbow')
ax1.contour(xh, yh, dbz, levels=[30], colors='k', linewidths=1)
ax1.set_xlim(xl)
ax1.set_ylim(yl)
ax1.set_title(f"Clusters, 30 dBZ (SIL = {sil1:.2f})")
p = ax2.scatter(x_mv[0,:]/1000, y_mv[0,:]/1000, s=1, c=labels1, marker='.', cmap='HomeyerRainbow')
ax2.contour(xh, yh, dbz0, levels=[30], colors='k', linewidths=1)
ax2.set_xlim(xl)
ax2.set_ylim(yl)
ax2.set_title(f"Cluster starting points, initial 30 dBZ")
plt.colorbar(p, ax=ax2)
plt.suptitle(f"{box_name} - {n} GMM trajectory clusters")
plt.show()



#%% Set source regions (need to combine this with the other code later)

if np.equal(box_name, 'MERGER (MV1)'):
    dbfile = open(ip+f"traj_clusters/traj_clusters_{stime/60:.0f}min_v2.pkl", 'rb')
    ccs = pickle.load(dbfile)
    cc = ccs['mv1']
    dbfile.close()
    
    
    c0 = (cc == 0) # low-level inflow
    c1 = (cc == 1) # mid-level inflow
    c2 = (cc == 2) # supercell cold pool
    c3 = (cc == 3) # QLCS outflow
    c4 = (cc == 4) # parcels starting in mesocyclone
    
    
    print(f"{len(pids[c0])} parcels from near-sfc inflow ({len(pids[c0])/len(pids)*100:.02f}%)")
    print(f"{len(pids[c1])} parcels from mid-level inflow ({len(pids[c1])/len(pids)*100:.02f}%)")
    print(f"{len(pids[c2])} parcels from supercell cold pool ({len(pids[c2])/len(pids)*100:.02f}%)")
    print(f"{len(pids[c3])} parcels from QLCS outflow ({len(pids[c3])/len(pids)*100:.02f}%)")
    print(f"{len(pids[c4])} parcels in updraft ({len(pids[c4])/len(pids)*100:.02f}%)")
    print(f"{len(pids)} total parcels")
    
    img_str = 'MV1'
    xl = [-60,30]
    yl = [-130,-40]


if np.equal(box_name, 'MERGER (MV2)'):
    dbfile = open(ip+f"traj_clusters/traj_clusters_{stime/60:.0f}min.pkl", 'rb')
    ccs = pickle.load(dbfile)
    cc = ccs['mv2']
    dbfile.close()
    
    c0 = (cc == 0) # low-level inflow
    c1 = (cc == 1) # mid-level inflow
    c2 = (cc == 2) # supercell cold pool
    c3 = (cc == 3) # QLCS low-level outflow
    c4 = (cc == 4) # QLCS outflow aloft
    
    print(f"{len(pids[c0])} parcels from near-sfc inflow ({len(pids[c0])/len(pids)*100:.02f}%)")
    print(f"{len(pids[c1])} parcels from mid-level inflow ({len(pids[c1])/len(pids)*100:.02f}%)")
    print(f"{len(pids[c2])} parcels from supercell cold pool ({len(pids[c2])/len(pids)*100:.02f}%)")
    print(f"{len(pids[c3])} parcels from QLCS low-level outflow ({len(pids[c3])/len(pids)*100:.02f}%)")
    print(f"{len(pids[c4])} parcels from QLCS outflow aloft ({len(pids[c4])/len(pids)*100:.02f}%)")
    print(f"{len(pids)} total parcels")
    
    img_str = 'MV2'
    xl = [-50,30]
    yl = [-100,-20]


if np.equal(box_name, 'QLCS'):
    dbfile = open(ip+f"traj_clusters_{stime/60:.0f}min_v2.pkl", 'rb')
    ccs = pickle.load(dbfile)
    cc = ccs['q']
    dbfile.close()
    
    c0 = (cc == 0) # low-level inflow
    c1 = (cc == 3) # QLCS outflow
    
    print(f"{len(pids[c0])} parcels from low-level inflow ({len(pids[c0])/len(pids)*100:.02f}%)")
    print(f"{len(pids[c1])} parcels from QLCS outflow ({len(pids[c1])/len(pids)*100:.02f}%)")
    print(f"{len(pids)} total parcels")
    
    img_str = 'Q'
    xl = [-50,30]
    yl = [-115,-35]

if np.equal(box_name, 'SUPERCELL'):
    dbfile = open(ip+f"traj_clusters_{stime/60:.0f}min.pkl", 'rb')
    ccs = pickle.load(dbfile)
    cc = ccs['s1']
    dbfile.close()
    
    c0 = (cc == 0) # low-level inflow
    c1 = (cc == 1) # mid-level inflow
    c2 = (cc == 2) # supercell outflow
    
    print(f"{len(pids[c0])} parcels from low-level inflow ({len(pids[c0])/len(pids)*100:.02f}%)")
    print(f"{len(pids[c1])} parcels from mid-level inflow ({len(pids[c1])/len(pids)*100:.02f}%)")
    print(f"{len(pids[c2])} parcels from supercell outflow ({len(pids[c2])/len(pids)*100:.02f}%)")
    print(f"{len(pids)} total parcels")
    
    img_str = 'S'
    xl = [-60,30]
    yl = [-130,-40]



#%% Save source regions and MV parcels to pkl

if False:
    dbfile = open(ip+f"traj_clusters/traj_clusters_{stime/60:.0f}min_v2.pkl", 'wb')
    ccs = {'mv1':cc}
    pickle.dump(ccs, dbfile)
    dbfile.close()

if False:
    dbfile = open(ip+f"traj_clusters/traj_clusters_{stime/60:.0f}min.pkl", 'rb')
    ccs = pickle.load(dbfile)
    dbfile.close()
    
    new_vars = {'mv1':cc}
    ccs.update(new_vars)
    dbfile = open(ip+f"traj_clusters/traj_clusters_{stime/60:.0f}min.pkl", 'wb')
    pickle.dump(ccs, dbfile)
    dbfile.close()


if False:
    cc = labels1
    cc_new = np.zeros(shape=cc.shape, dtype=cc.dtype)
    cc_new[(cc == 0)] = 1
    cc_new[(cc == 1)] = 0
    # cc_new[(cc == 2)] = 2
    # cc_new[(cc == 3)] = 1
    # cc_new[(cc == 4)] = 0
    # cc_new[(labels1 == 3)] = 1
    # cc_new[(cc == 4)] = 3
    
    # cc_new[(cc==1) & (x_mv[0,:]>-15000)] = 0
    # cc_new[(z_mv[0,:]>1100)] = 2
    # cc_new[(x_mv[0,:]<-30000)] = 3
    # cc_new[(cc_new == 1) & (b_mv[0,:] > -0.03)] = 0
    # cc_new[(cc_new == 4) & (b_mv[0,:] < -0.01) & (y_mv[0,:] > -105000)] = 1
    # cc_new[(labels1 == 1) & (x_mv[0,:] > -10000)] = 0
    
    cc_old = cc
    cc = cc_new

if False:
    cc_new = np.zeros(shape=cc.shape, dtype=cc.dtype)
    c0 = (cc == 0)
    c1 = (cc == 1)
    c2 = (cc == 2)
    c3 = (cc == 3)
    c4 = (cc == 4)
    
    cc_new[(cc == 0)] = 0
    cc_new[(cc == 1)] = 1
    cc_new[(cc == 2)] = 2
    cc_new[(cc == 3)] = 3
    cc_new[(cc == 4)] = 4
    
    cc_new[(z_mv[0,:] >= 1500)] = 1
    
    cc_old = cc
    cc = cc_new
    
    print(f"{len(pids[(cc == 1)])} parcels from mid-level inflow ({len(pids[(cc == 1)])/len(pids)*100:.02f}%)")
    

if False:
    img_str = 'S'
    if fnum == 28:
        dbfile = open(f"/Users/morgan.schneider/Documents/merger/traj_{img_str}.pkl", 'wb')
        traj_vars = {'pids':pids, 'x':x_mv, 'y':y_mv, 'z':z_mv, 'w':w_mv, 'zvort':zvort_mv,
                     'b':b_mv}
        traj = {f"{stime/60:.0f}min":traj_vars}
        pickle.dump(traj, dbfile)
        dbfile.close()
    else:
        dbfile = open(f"/Users/morgan.schneider/Documents/merger/traj_{img_str}.pkl", 'rb')
        traj = pickle.load(dbfile)
        dbfile.close()
        
        # traj_time = traj[f"{stime/60:.0f}min"]
        # traj_time.update({'u':u_mv, 'v':v_mv})
        # traj.update({f"{stime/60:.0f}min":traj_time})
        
        traj_time = {'pids':pids, 'x':x_mv, 'y':y_mv, 'z':z_mv, 'w':w_mv, 'zvort':zvort_mv,
                     'b':b_mv}
        traj.update({f"{stime/60:.0f}min":traj_time})
        # traj[f"{stime/60:.0f}min"].update({'u':u_mv, 'v':v_mv})
        dbfile = open(f"/Users/morgan.schneider/Documents/merger/traj_{img_str}.pkl", 'wb')
        pickle.dump(traj, dbfile)
        dbfile.close()


if False:
    cc = np.zeros(shape=(len(pids),), dtype=int)
    cc[(z_mv[0,:] > 1000)] = 1
    
    # dbfile = open(ip+f"traj_clusters_{stime/60:.0f}min_v2.pkl", 'rb')
    # ccs = pickle.load(dbfile)
    # dbfile.close()
    
    # new_vars = {'q':cc}
    # ccs.update(new_vars)
    dbfile = open(ip+f"traj_clusters_{stime/60:.0f}min_v2.pkl", 'wb')
    pickle.dump(cc, dbfile)
    dbfile.close()
    
    
    

#%% Make random plots

qix = int(np.round(len(pids)/90))

figsave = False

xl = [-60,30]
yl = [-130,-40]

# Plot trajectories over dbz
if False:
    qix = 20
    
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    plot_cfill(xf, yf, np.ma.masked_array(dbz, dbz<1), 'dbz', ax, datalims=[1,70], xlims=xl, ylims=yl)
    ax.set_xlabel('E-W distance (km)')
    ax.set_ylabel('N-S distance (km)')
    ax.set_title(f"{box_name} - 10 m $Z_H$ at {stime/60:.0f} min")
    ax.plot(x_mv[:,::qix]/1000, y_mv[:,::qix]/1000, '-k', linewidth=1)
    if figsave:
        plt.savefig(ip+f"imgs_trajectories/{img_str}_traj_{stime/60:.0f}min.png", dpi=300)


# Plot trajectories by buoyancy and height, traj starts by initial buoyancy and height
if True:
    qix = 1
    
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    p = ax.scatter(x_mv[:ti,::qix]/1000, y_mv[:ti,::qix]/1000, s=1, c=b_mv[:ti,::qix], marker='.', cmap='HomeyerRainbow', vmin=-0.2, vmax=0.1)
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
    p = ax.scatter(x_mv[:ti,::qix]/1000, y_mv[:ti,::qix]/1000, s=1, c=z_mv[:ti,::qix]/1000, marker='.', cmap='HomeyerRainbow', vmin=0, vmax=2)
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
    p = ax.scatter(x_mv[0,:]/1000, y_mv[0,:]/1000, s=5, c=b_mv[0,:], marker='.', cmap='HomeyerRainbow', vmin=-0.2, vmax=0.1)
    ax.contour(xh, yh, dbz, levels=[30], colors='k', linewidths=1)
    ax.set_xlim(xl)
    ax.set_ylim(yl)
    ax.set_title(f"{box_name} - Starting buoyancy for parcels at {stime/60:.0f} min")
    plt.colorbar(p, ax=ax)
    ax.set_xlabel('E-W distance (km)')
    ax.set_ylabel('N-S distance (km)')
    if figsave:
        plt.savefig(ip+f"imgs_trajectories/{img_str}_traj-buoystart_{stime/60:.0f}min.png", dpi=300)
    
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    p = ax.scatter(x_mv[0,:]/1000, y_mv[0,:]/1000, s=5, c=z_mv[0,:]/1000, marker='.', cmap='HomeyerRainbow', vmin=0, vmax=2)
    ax.contour(xh, yh, dbz, levels=[30], colors='k', linewidths=1)
    ax.set_xlim(xl)
    ax.set_ylim(yl)
    ax.set_title(f"{box_name} - Starting height for parcels at {stime/60:.0f} min")
    plt.colorbar(p, ax=ax)
    ax.set_xlabel('E-W distance (km)')
    ax.set_ylabel('N-S distance (km)')
    if figsave:
        plt.savefig(ip+f"imgs_trajectories/{img_str}_traj-zstart_{stime/60:.0f}min.png", dpi=300)
    
    # fig,ax = plt.subplots(1,1,figsize=(8,6))
    # p = ax.scatter(x_mv[0,:]/1000, y_mv[0,:]/1000, s=5, c=w_mv[0,:], marker='.', cmap='HomeyerRainbow', vmin=-10, vmax=10)
    # ax.contour(xh, yh, dbz0, levels=[30], colors='k', linewidths=1)
    # ax.set_xlim(xl)
    # ax.set_ylim(yl)
    # ax.set_title(f"{box_name} - Starting vertical velocity for parcels at {stime/60:.0f} min")
    # plt.colorbar(p, ax=ax)
    # ax.set_xlabel('E-W distance (km)')
    # ax.set_ylabel('N-S distance (km)')
    # if figsave:
    #     plt.savefig(ip+f"imgs_trajectories/{img_str}_traj-wstart_{stime/60:.0f}min.png", dpi=300)
    
    plt.show()


# Plot trajectories colored by source
if False:
    fig,(ax1,ax2) = plt.subplots(1,2,figsize=(12,5))
    ax1.scatter(x_mv[:ti,:]/1000, y_mv[:ti,:]/1000, s=1, c=np.tile(cc, (len(x_mv[:ti,0]),1)), marker='.', cmap='HomeyerRainbow', vmin=0, vmax=4)
    ax1.contour(xh, yh, dbz, levels=[30], colors='k', linewidths=1)
    ax1.set_xlim(xl)
    ax1.set_ylim(yl)
    ax1.set_title(f"Trajectory clusters, 30 dBZ")
    p = ax2.scatter(x_mv[0,:]/1000, y_mv[0,:]/1000, s=1, c=cc, marker='.', cmap='HomeyerRainbow', vmin=0, vmax=4)
    ax2.contour(xh, yh, dbz0, levels=[30], colors='k', linewidths=1)
    # ax2.contour(xh, yh, winterp, levels=[5], colors='r')
    ax2.set_xlim(xl)
    ax2.set_ylim(yl)
    ax2.set_title(f"Trajectory starting points, initial 30 dBZ")
    plt.colorbar(p, ax=ax2)
    plt.suptitle(f"{box_name} - Trajectory clusters at {stime/60:.0f} min")
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
    ax.set_title(f"{box_name} - Source A - 10 m $Z_H$ at {stime/60:.0f} min")
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
    p = ax.scatter(x_src[:,::qix]/1000, y_src[:,::qix]/1000, s=1, c=z_src[:,::qix]/1000, marker='.', cmap='HomeyerRainbow', vmin=0, vmax=3.5)
    ax.contour(xh, yh, dbz0, levels=[30], colors='k', linewidths=1)
    ax.set_xlabel('E-W distance (km)')
    ax.set_ylabel('N-S distance (km)')
    ax.set_xlim(xl)
    ax.set_ylim(yl)
    ax.set_title(f"{box_name} - Trajectories by $z$ at {stime/60:.0f} min")
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
    ax.set_title(f"{box_name} - Trajectories by $w$ at {stime/60:.0f} min")
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
    ax.set_title(f"{box_name} - Trajectories by zvort at {stime/60:.0f} min")
    plt.colorbar(p, ax=ax)
    if figsave:
        plt.savefig(ip+f"imgs_trajectories/{img_str}_src1_traj-zvort_{stime/60:.0f}min.png", dpi=300)
    
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    p = ax.scatter(x_src[:,::qix]/1000, y_src[:,::qix]/1000, s=1, c=b_src[:,::qix], marker='.', cmap='HomeyerRainbow', vmin=-0.4, vmax=0.1)
    ax.contour(xh, yh, dbz, levels=[30], colors='k', linewidths=1)
    ax.set_xlabel('E-W distance (km)')
    ax.set_ylabel('N-S distance (km)')
    ax.set_xlim(xl)
    ax.set_ylim(yl)
    ax.set_title(f"{box_name} - Trajectories by buoyancy at {stime/60:.0f} min")
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
        ax.set_title(f"{box_name} - 10 m $Z_H$ at {stime2/60:.0f} min")
        ax.plot(x_mv[0:tj+1,::qix]/1000, y_mv[0:tj+1,::qix]/1000, '-k', linewidth=1)
        plt.show()
        if figsave:
            plt.savefig(ip+f"imgs_trajectories/{img_str}_traj_{stime2/60:.0f}min.png", dpi=300)
    
    # Trajectories colored by z, w, zvort, buoyancy
    if False:
        qix = 1
        
        fig,ax = plt.subplots(1,1,figsize=(8,6))
        p = ax.scatter(x_mv[0:tj+1,::qix]/1000, y_mv[0:tj+1,::qix]/1000, s=1, c=z_mv[0:tj+1,::qix]/1000, marker='.', cmap='HomeyerRainbow', vmin=0, vmax=3.5)
        ax.set_xlabel('E-W distance (km)')
        ax.set_ylabel('N-S distance (km)')
        ax.set_xlim(xl)
        ax.set_ylim(yl)
        ax.set_title(f"{box_name} - Trajectories by $z$ at {stime2/60:.0f} min")
        plt.show()
        if figsave:
            plt.savefig(ip+f"imgs_trajectories/{img_str}_traj-z_{stime2/60:.0f}min.png", dpi=300)
        
        fig,ax = plt.subplots(1,1,figsize=(8,6))
        p = ax.scatter(x_mv[0:tj+1,::qix]/1000, y_mv[0:tj+1,::qix]/1000, s=1, c=w_mv[0:tj+1,::qix], marker='.', cmap='RdBu_r', vmin=-10, vmax=10)
        ax.set_xlabel('E-W distance (km)')
        ax.set_ylabel('N-S distance (km)')
        ax.set_xlim(xl)
        ax.set_ylim(yl)
        ax.set_title(f"{box_name} - Trajectories by $w$ at {stime2/60:.0f} min")
        plt.show()
        if figsave:
            plt.savefig(ip+f"imgs_trajectories/{img_str}_traj-w_{stime2/60:.0f}min.png", dpi=300)
        
        fig,ax = plt.subplots(1,1,figsize=(8,6))
        p = ax.scatter(x_mv[0:tj+1,::qix]/1000, y_mv[0:tj+1,::qix]/1000, s=1, c=zvort_mv[0:tj+1,::qix], marker='.', cmap='RdBu_r', vmin=-0.04, vmax=0.04)
        ax.set_xlabel('E-W distance (km)')
        ax.set_ylabel('N-S distance (km)')
        ax.set_xlim(xl)
        ax.set_ylim(yl)
        ax.set_title(f"{box_name} - Trajectories by zvort at {stime2/60:.0f} min")
        plt.show()
        if figsave:
            plt.savefig(ip+f"imgs_trajectories/{img_str}_traj-zvort_{stime2/60:.0f}min.png", dpi=300)
        
        fig,ax = plt.subplots(1,1,figsize=(8,6))
        p = ax.scatter(x_mv[0:tj+1,::qix]/1000, y_mv[0:tj+1,::qix]/1000, s=1, c=b_mv[0:tj+1,::qix], marker='.', cmap='HomeyerRainbow', vmin=-0.4, vmax=0.1)
        ax.set_xlabel('E-W distance (km)')
        ax.set_ylabel('N-S distance (km)')
        ax.set_xlim(xl)
        ax.set_ylim(yl)
        ax.set_title(f"{box_name} - Trajectories by buoyancy at {stime2/60:.0f} min")
        plt.show()
        if figsave:
            plt.savefig(ip+f"imgs_trajectories/{img_str}_traj-buoy_{stime2/60:.0f}min.png", dpi=300)


#%% Time series plotting (single time)

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


#%% Rearrange and change around cluster numbers for convenience

if False:
    for i in [195,200,205,210,215,220,225,230,235,239]:
        dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/traj_clusters/traj_clusters_{i:.0f}min.pkl", 'rb')
        cc = pickle.load(dbfile)
        cc_old = cc['mv1']
        cc_mv1 = np.zeros(shape=cc_old.shape, dtype=cc_old.dtype)
        cc_mv1[(cc_old == 0)] = 0
        cc_mv1[(cc_old == 1)] = 2
        cc_mv1[(cc_old == 2)] = 1
        cc_mv1[(cc_old == 3)] = 3
        cc_mv1[(cc_old == 4)] = 4
        
        cc_old2 = cc['mv2']
        cc_mv2 = np.zeros(shape=cc_old2.shape, dtype=cc_old2.dtype)
        cc_mv2[(cc_old2 == 0)] = 0
        cc_mv2[(cc_old2 == 1)] = 2
        cc_mv2[(cc_old2 == 2)] = 1
        cc_mv2[(cc_old2 == 3)] = 3
        cc_mv2[(cc_old2 == 4)] = 4
        dbfile.close()
        
        cc_new = {'mv1':cc_mv1, 'mv2':cc_mv2}
        dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/traj_clusters/traj_clusters_{i:.0f}min_v2.pkl", 'wb')
        pickle.dump(cc_new, dbfile)
        dbfile.close()
        
        
        
        dbfile = open(f"/Users/morgan.schneider/Documents/merger/qlcs-125m/traj_clusters/traj_clusters_{i:.0f}min.pkl", 'rb')
        cc = pickle.load(dbfile)
        cc_oldq = cc['q']
        cc_q = np.zeros(shape=cc_oldq.shape, dtype=cc_oldq.dtype)
        cc_q[(cc_oldq == 0)] = 0
        cc_q[(cc_oldq == 1)] = 3
        
        cc_new = {'q':cc_q}
        dbfile = open(f"/Users/morgan.schneider/Documents/merger/qlcs-125m/traj_clusters/traj_clusters_{i:.0f}min_v2.pkl", 'wb')
        pickle.dump(cc_new, dbfile)
        dbfile.close()
    



#%% Make big multipanel trajectory figure ***PAPER FIG***

import matplotlib as mpl
cm = mpl.colors.ListedColormap(['gold', 'red', 'deepskyblue', 'mediumblue'])

use_SR = True


dbfile = open('/Users/morgan.schneider/Documents/merger/traj_MV1.pkl', 'rb')
traj = pickle.load(dbfile)
dbfile.close()

dbfile = open('/Users/morgan.schneider/Documents/merger/merger-125m/boxes_s1.pkl', 'rb')
box = pickle.load(dbfile)
xb1 = box['x1_pp']
xb2 = box['x2_pp']
yb1 = box['y1_pp']
yb2 = box['y2_pp']
dbfile.close()


fnums = [33, 43, 53]

ds = nc.Dataset('/Volumes/Promise_Pegasus_70TB/merger/merger-125m/cm1out_pdata.nc')
ptime = ds.variables['time'][:].data
ds.close()


if use_SR:
    from scipy.interpolate import interp1d
    dbfile = open('/Users/morgan.schneider/Documents/merger/merger-125m/storm_motion.pkl', 'rb')
    sm = pickle.load(dbfile)
    u_storm = sm['u_storm']
    v_storm = sm['v_storm']
    dbfile.close
    stimes = np.linspace(10800, 14400, 61)
    
    # Interpolate storm motion from model output timesto parcel output times (60 s to 15 s)
    fu = interp1d(stimes, u_storm)
    u_storm_interp = fu(ptime)
    fv = interp1d(stimes, v_storm)
    v_storm_interp = fv(ptime)
    
    xl = [-52,32]; yl = [-92,-8]
else:
    xl = [-62,22]; yl = [-122,-38]


ds = nc.Dataset('/Volumes/Promise_Pegasus_70TB/merger/merger-125m/cm1out_000033.nc')
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
ix1 = np.where(xh >= xl[0])[0][0]; ix2 = np.where(xh >= xl[1])[0][1]; ix = slice(ix1,ix2)
iy1 = np.where(yh >= yl[0])[0][0]; iy2 = np.where(yh >= yl[1])[0][1]; iy = slice(iy1,iy2)
ds.close()




fig,axs = plt.subplots(3, 3, figsize=(9,8), subplot_kw=dict(box_aspect=1), sharex=True, sharey=True, layout='constrained')


for i in range(3):
    ds = nc.Dataset(f"/Volumes/Promise_Pegasus_70TB/merger/merger-125m/cm1out_{fnums[i]:06d}.nc")
    stime = ds.variables['time'][0]
    t = stime/60
    dbz = ds.variables['dbz'][:].data[0,0,iy,ix]
    ds.close()
    
    dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/traj_clusters/traj_clusters_{t:.0f}min_v2.pkl", 'rb')
    cc = pickle.load(dbfile)
    cc_mv = cc['mv1']
    if i == 2:
        cc_mv[(cc_mv == 4)] = 0
    dbfile.close()
    
    ti = np.where(ptime == stime)[0][0]
    
    pids_mv = traj[f"{t:.0f}min"]['pids'][(cc_mv <= 3)]
    x_mv = traj[f"{t:.0f}min"]['x'][:ti+1, (cc_mv <= 3)]
    y_mv = traj[f"{t:.0f}min"]['y'][:ti+1, (cc_mv <= 3)]
    z_mv = traj[f"{t:.0f}min"]['z'][:ti+1, (cc_mv <= 3)]
    w_mv = traj[f"{t:.0f}min"]['w'][:ti+1, (cc_mv <= 3)]
    b_mv = traj[f"{t:.0f}min"]['b'][:ti+1, (cc_mv <= 3)]
    cc_mv = cc_mv[(cc_mv <= 3)]
    
    
    
    # 3x3 panel trajectories colored by z, B, and source for merger
    
    x0 = x_mv[:,(cc_mv==0)]; y0 = y_mv[:,(cc_mv==0)]; z0 = z_mv[:,(cc_mv==0)]; b0 = b_mv[:,(cc_mv==0)]; c0 = cc_mv[(cc_mv==0)]
    x1 = x_mv[:,(cc_mv==1)]; y1 = y_mv[:,(cc_mv==1)]; z1 = z_mv[:,(cc_mv==1)]; b1 = b_mv[:,(cc_mv==1)]; c1 = cc_mv[(cc_mv==1)]
    x2 = x_mv[:,(cc_mv==2)]; y2 = y_mv[:,(cc_mv==2)]; z2 = z_mv[:,(cc_mv==2)]; b2 = b_mv[:,(cc_mv==2)]; c2 = cc_mv[(cc_mv==2)]
    x3 = x_mv[:,(cc_mv==3)]; y3 = y_mv[:,(cc_mv==3)]; z3 = z_mv[:,(cc_mv==3)]; b3 = b_mv[:,(cc_mv==3)]; c3 = cc_mv[(cc_mv==3)]
    if i == 0: # was ::50, ::8 for 195 min
        xp = np.append(x0[:,::51], x2[:,::5], axis=1)
        yp = np.append(y0[:,::51], y2[:,::5], axis=1)
        zp = np.append(z0[:,::51], z2[:,::5], axis=1)
        bp = np.append(b0[:,::51], b2[:,::5], axis=1)
        cc = np.append(c0[::51], c2[::5])
    if i == 1:
        xp = np.append(np.append(np.append(x0[:,::15], x1[:,::3], axis=1), x2[:,::3], axis=1), x3[:,:], axis=1)
        yp = np.append(np.append(np.append(y0[:,::15], y1[:,::3], axis=1), y2[:,::3], axis=1), y3[:,:], axis=1)
        zp = np.append(np.append(np.append(z0[:,::15], z1[:,::3], axis=1), z2[:,::3], axis=1), z3[:,:], axis=1)
        bp = np.append(np.append(np.append(b0[:,::15], b1[:,::3], axis=1), b2[:,::3], axis=1), b3[:,:], axis=1)
        cc = np.append(np.append(np.append(c0[::15], c1[::3]), c2[::3]), c3[:])
    if i == 2: # was ::3, :, :, ::4 for 225 min
        xp = np.append(np.append(np.append(x0[:,::3], x1[:,::2], axis=1), x2[:,:], axis=1), x3[:,::6], axis=1)
        yp = np.append(np.append(np.append(y0[:,::3], y1[:,::2], axis=1), y2[:,:], axis=1), y3[:,::6], axis=1)
        zp = np.append(np.append(np.append(z0[:,::3], z1[:,::2], axis=1), z2[:,:], axis=1), z3[:,::6], axis=1)
        bp = np.append(np.append(np.append(b0[:,::3], b1[:,::2], axis=1), b2[:,:], axis=1), b3[:,::6], axis=1)
        cc = np.append(np.append(np.append(c0[::3], c1[::2]), c2[:]), c3[::6])
    
    if use_SR:
        xp,yp = get_SR_positions(xp, yp, u_storm_interp[:ti+1], v_storm_interp[:ti+1], 15)
    
    
    
    # 15*i+14 for 195/210/225, 10*i+19 for 200/210/220
    
    # p1 = axs[0,i].scatter(x_mv/1000, y_mv/1000, s=1, c=z_mv/1000, marker='.', cmap='HomeyerRainbow', vmin=0, vmax=3)
    p1 = axs[0,i].scatter(xp/1000, yp/1000, s=0.5, c=zp/1000, marker='.', cmap='HomeyerRainbow', vmin=0, vmax=3)
    axs[0,i].contour(xh[ix], yh[iy], dbz, levels=[30], colors='k', linewidths=1)
    r1 = patches.Rectangle((xb1[10*i+19],yb1[10*i+19]), 10, 10, fc='none', ec='k', lw=2)
    axs[0,i].add_patch(r1)
    axs[0,i].set_xlim(xl)
    axs[0,i].set_ylim(yl)
    axs[0,i].set_title(f"{t:.0f} min", fontsize=16)
    
    # p2 = axs[1,i].scatter(x_mv/1000, y_mv/1000, s=1, c=b_mv, marker='.', cmap='HomeyerRainbow', vmin=-0.3, vmax=0.1)
    p2 = axs[1,i].scatter(xp/1000, yp/1000, s=0.5, c=bp, marker='.', cmap='HomeyerRainbow', vmin=-0.3, vmax=0.1)
    axs[1,i].contour(xh[ix], yh[iy], dbz, levels=[30], colors='k', linewidths=1)
    r2 = patches.Rectangle((xb1[10*i+19],yb1[10*i+19]), 10, 10, fc='none', ec='k', lw=2)
    axs[1,i].add_patch(r2)
    axs[1,i].set_xlim(xl)
    axs[1,i].set_ylim(yl)
    
    p3 = axs[2,i].scatter([0,0,0,0], [0,0,0,0], s=1, c=[0,1,2,3], marker='.', cmap=cm, vmin=-0.5, vmax=3.5)
    # axs[2,i].scatter(x_mv/1000, y_mv/1000, s=1, c=np.tile(cc_mv,(len(x_mv[:,0]),1)), marker='.', cmap=cm, vmin=-0.5, vmax=3.5)
    axs[2,i].scatter(xp/1000, yp/1000, s=0.5, c=np.tile(cc,(len(xp[:,0]),1)), marker='.', cmap=cm, vmin=-0.5, vmax=3.5)
    axs[2,i].contour(xh[ix], yh[iy], dbz, levels=[30], colors='k', linewidths=1)
    r3 = patches.Rectangle((xb1[10*i+19],yb1[10*i+19]), 10, 10, fc='none', ec='k', lw=2)
    axs[2,i].add_patch(r3)
    axs[2,i].set_xlim(xl)
    axs[2,i].set_ylim(yl)
    axs[2,i].set_xlabel('x (km)', fontsize=14)
    


if False:
    l1, = axs[2,0].plot([xl[0]-2,xl[0]-1], [yl[0]-2,yl[0]-1], color='gold')
    l2, = axs[2,0].plot([xl[0]-2,xl[0]-1], [yl[0]-2,yl[0]-1], color='red')
    l3, = axs[2,0].plot([xl[0]-2,xl[0]-1], [yl[0]-2,yl[0]-1], color='deepskyblue')
    l4, = axs[2,0].plot([xl[0]-2,xl[0]-1], [yl[0]-2,yl[0]-1], color='mediumblue')
    axs[2,0].legend(handles=[l1,l2,l3,l4],
                    labels=['Low-level inflow', 'Mid-level inflow',
                            'Supercell outflow', 'QLCS outflow'], loc='upper right', fontsize=8)

cb1 = plt.colorbar(p1, ax=axs[0,2], extend='both', ticks=np.linspace(0,3,7))
cb1.set_label("Height (km)", fontsize=12)
# tl = cb1.ax.get_yticklabels()
# cb1.ax.set_yticklabels(tl, ha='right')
# cb1.ax.yaxis.set_tick_params(pad=30)

cb2 = plt.colorbar(p2, ax=axs[1,2], extend='both', ticks=np.linspace(-0.3,0.1,5))
cb2.set_label("Buoyancy (m s$^{-2}$)", fontsize=12)
# tl = cb2.ax.get_yticklabels()
# cb2.ax.set_yticklabels(tl, ha='right')
# cb2.ax.yaxis.set_tick_params(pad=30)

cb3 = plt.colorbar(p3, ax=axs[2,2], ticks=[0,1,2,3])
cb3.ax.set_yticklabels([' Low-level\n inflow', ' Mid-level\n inflow', 
                        ' Supercell\n outflow', ' QLCS\n outflow'])
# cb3.set_label("Source region", fontsize=12)

axs[0,0].set_yticks(np.linspace(yl[0]+2, yl[1]-2, 5))
axs[1,0].set_yticks(np.linspace(yl[0]+2, yl[1]-2, 5))
axs[2,0].set_yticks(np.linspace(yl[0]+2, yl[1]-2, 5))

axs[0,0].set_ylabel('y (km)', fontsize=14)
axs[1,0].set_ylabel('y (km)', fontsize=14)
axs[2,0].set_ylabel('y (km)', fontsize=14)

axs[2,0].set_xlabel('x (km)', fontsize=14)
axs[2,1].set_xlabel('x (km)', fontsize=14)
axs[2,2].set_xlabel('x (km)', fontsize=14)

figsave = False

if figsave:
    fig.savefig('/Users/morgan.schneider/Documents/merger/traj_z+B+src_SR.png', dpi=300)


#%% Same plot but for qlcsonly ***PAPER FIG***

import matplotlib as mpl
cm = mpl.colors.ListedColormap(['gold', 'mediumblue'])

use_SR = True


dbfile = open(f"/Users/morgan.schneider/Documents/merger/traj_Q.pkl", 'rb')
traj = pickle.load(dbfile)
dbfile.close()

dbfile = open(f"/Users/morgan.schneider/Documents/merger/qlcs-125m/boxes_q.pkl", 'rb')
box = pickle.load(dbfile)
xb1 = box['x1_pp']
xb2 = box['x2_pp']
yb1 = box['y1_pp']
yb2 = box['y2_pp']
dbfile.close()


fnums = [33, 43, 53]

ds = nc.Dataset('/Volumes/Promise_Pegasus_70TB/merger/qlcs-125m/cm1out_pdata.nc')
ptime = ds.variables['time'][:].data
ds.close()


if use_SR:
    from scipy.interpolate import interp1d
    dbfile = open('/Users/morgan.schneider/Documents/merger/qlcs-125m/storm_motion.pkl', 'rb')
    sm = pickle.load(dbfile)
    u_storm = sm['u_storm']
    v_storm = sm['v_storm']
    dbfile.close
    stimes = np.linspace(10800, 14400, 61)
    
    # Interpolate storm motion from model output timesto parcel output times (60 s to 15 s)
    fu = interp1d(stimes, u_storm)
    u_storm_interp = fu(ptime)
    fv = interp1d(stimes, v_storm)
    v_storm_interp = fv(ptime)
    
    xl = [-52,32]; yl = [-102,-18]
else:
    xl = [-62,22]; yl = [-122,-38]


ds = nc.Dataset('/Volumes/Promise_Pegasus_70TB/merger/qlcs-125m/cm1out_000033.nc')
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
ix1 = np.where(xh >= xl[0])[0][0]; ix2 = np.where(xh >= xl[1])[0][1]; ix = slice(ix1,ix2)
iy1 = np.where(yh >= yl[0])[0][0]; iy2 = np.where(yh >= yl[1])[0][1]; iy = slice(iy1,iy2)
ds.close()




fig,axs = plt.subplots(3, 3, figsize=(9,8), subplot_kw=dict(box_aspect=1), sharex=True, sharey=True, layout='constrained')


for i in range(3):
    ds = nc.Dataset(f"/Volumes/Promise_Pegasus_70TB/merger/qlcs-125m/cm1out_{fnums[i]:06d}.nc")
    stime = ds.variables['time'][0]
    t = stime/60
    dbz = ds.variables['dbz'][:].data[0,0,iy,ix]
    ds.close()
    
    dbfile = open(f"/Users/morgan.schneider/Documents/merger/qlcs-125m/traj_clusters_{t:.0f}min_v2.pkl", 'rb')
    cc = pickle.load(dbfile)
    cc_mv = cc['q']
    dbfile.close()
    
    ti = np.where(ptime == stime)[0][0]
    
    pids_mv = traj[f"{t:.0f}min"]['pids']
    x_mv = traj[f"{t:.0f}min"]['x'][:ti+1,:]
    y_mv = traj[f"{t:.0f}min"]['y'][:ti+1,:]
    z_mv = traj[f"{t:.0f}min"]['z'][:ti+1,:]
    w_mv = traj[f"{t:.0f}min"]['w'][:ti+1,:]
    b_mv = traj[f"{t:.0f}min"]['b'][:ti+1,:]

    
    
    # 3x3 panel trajectories colored by z, B, and source for qlcsonly
    x0 = x_mv[:,(cc_mv==0)]; y0 = y_mv[:,(cc_mv==0)]; z0 = z_mv[:,(cc_mv==0)]; b0 = b_mv[:,(cc_mv==0)]; c0 = cc_mv[(cc_mv==0)]
    x3 = x_mv[:,(cc_mv==3)]; y3 = y_mv[:,(cc_mv==3)]; z3 = z_mv[:,(cc_mv==3)]; b3 = b_mv[:,(cc_mv==3)]; c3 = cc_mv[(cc_mv==3)]
    
    
    if i == 0: # was ::50, ::10 for total ~35
        xp = np.append(x0[:,::40], x3[:,::7], axis=1)
        yp = np.append(y0[:,::40], y3[:,::7], axis=1)
        zp = np.append(z0[:,::40], z3[:,::7], axis=1)
        bp = np.append(b0[:,::40], b3[:,::7], axis=1)
        cc = np.append(c0[::40], c3[::7])
    if i == 1: # was ::30, ::8
        xp = np.append(x0[:,::30], x3[:,::7], axis=1)
        yp = np.append(y0[:,::30], y3[:,::7], axis=1)
        zp = np.append(z0[:,::30], z3[:,::7], axis=1)
        bp = np.append(b0[:,::30], b3[:,::7], axis=1)
        cc = np.append(c0[::30], c3[::7])
    if i == 2: # was ::8, ::8
        xp = np.append(x0[:,::7], x3[:,::7], axis=1)
        yp = np.append(y0[:,::7], y3[:,::7], axis=1)
        zp = np.append(z0[:,::7], z3[:,::7], axis=1)
        bp = np.append(b0[:,::7], b3[:,::7], axis=1)
        cc = np.append(c0[::7], c3[::7])
    
    if use_SR:
        xp,yp = get_SR_positions(xp, yp, u_storm_interp[:ti+1], v_storm_interp[:ti+1], 15)
    
    
    
    p1 = axs[0,i].scatter(xp/1000, yp/1000, s=0.5, c=zp/1000, marker='.', cmap='HomeyerRainbow', vmin=0, vmax=3)
    axs[0,i].contour(xh[ix], yh[iy], dbz, levels=[30], colors='k', linewidths=1)
    r1 = patches.Rectangle((xb1[10*i+19],yb1[10*i+19]), 10, 10, fc='none', ec='k', lw=2)
    axs[0,i].add_patch(r1)
    axs[0,i].set_xlim(xl)
    axs[0,i].set_ylim(yl)
    axs[0,i].set_title(f"{t:.0f} min", fontsize=16)
    
    p2 = axs[1,i].scatter(xp/1000, yp/1000, s=0.5, c=bp, marker='.', cmap='HomeyerRainbow', vmin=-0.3, vmax=0.1)
    axs[1,i].contour(xh[ix], yh[iy], dbz, levels=[30], colors='k', linewidths=1)
    r2 = patches.Rectangle((xb1[10*i+19],yb1[10*i+19]), 10, 10, fc='none', ec='k', lw=2)
    axs[1,i].add_patch(r2)
    axs[1,i].set_xlim(xl)
    axs[1,i].set_ylim(yl)
    
    p3 = axs[2,i].scatter([0,0], [0,0], s=1, c=[0,3], marker='.', cmap=cm, vmin=-1.5, vmax=4.5)
    axs[2,i].scatter(xp/1000, yp/1000, s=0.5, c=np.tile(cc,(len(xp[:,0]),1)), marker='.', cmap=cm, vmin=-1.5, vmax=4.5)
    axs[2,i].contour(xh[ix], yh[iy], dbz, levels=[30], colors='k', linewidths=1)
    r3 = patches.Rectangle((xb1[10*i+19],yb1[10*i+19]), 10, 10, fc='none', ec='k', lw=2)
    axs[2,i].add_patch(r3)
    axs[2,i].set_xlim(xl)
    axs[2,i].set_ylim(yl)
    axs[2,i].set_xlabel('x (km)', fontsize=14)

if False:
    l1, = axs[2,0].plot([xl[0]-2,xl[0]-1], [yl[0]-2,yl[0]-1], color='gold')
    l2, = axs[2,0].plot([xl[0]-2,xl[0]-1], [yl[0]-2,yl[0]-1], color='mediumblue')
    axs[2,0].legend(handles=[l1,l2], labels=['Low-level inflow', 'Low-level outflow'], loc='upper right', fontsize=8)

cb1 = plt.colorbar(p1, ax=axs[0,2], extend='both', ticks=np.linspace(0,3,7))
cb1.set_label("Height (km)", fontsize=12)
# tl = cb1.ax.get_yticklabels()
# cb1.ax.set_yticklabels(tl, ha='right')
# cb1.ax.yaxis.set_tick_params(pad=30)

cb2 = plt.colorbar(p2, ax=axs[1,2], extend='both', ticks=np.linspace(-0.3,0.1,5))
cb2.set_label("Buoyancy (m s$^{-2}$)", fontsize=12)
# tl = cb2.ax.get_yticklabels()
# cb2.ax.set_yticklabels(tl, ha='right')
# cb2.ax.yaxis.set_tick_params(pad=30)

cb3 = plt.colorbar(p3, ax=axs[2,2], ticks=[0,3])
cb3.ax.set_yticklabels([' Low-level\n inflow', ' QLCS\n outflow'])
# cb3.set_label("Source region", fontsize=12)

axs[0,0].set_yticks(np.linspace(yl[0]+2, yl[1]-2, 5))
axs[1,0].set_yticks(np.linspace(yl[0]+2, yl[1]-2, 5))
axs[2,0].set_yticks(np.linspace(yl[0]+2, yl[1]-2, 5))

axs[0,0].set_ylabel('y (km)', fontsize=14)
axs[1,0].set_ylabel('y (km)', fontsize=14)
axs[2,0].set_ylabel('y (km)', fontsize=14)

axs[2,0].set_xlabel('x (km)', fontsize=14)
axs[2,1].set_xlabel('x (km)', fontsize=14)
axs[2,2].set_xlabel('x (km)', fontsize=14)

figsave = True

if figsave:
    fig.savefig('/Users/morgan.schneider/Documents/merger/traj_z+B+src_QLCS_SR.png', dpi=300)




#%% Same plot but for supercell ***MAYBE PAPER FIG***

import matplotlib as mpl
cm = mpl.colors.ListedColormap(['gold', 'red', 'deepskyblue'])

use_SR = True


dbfile = open(f"/Users/morgan.schneider/Documents/merger/traj_S.pkl", 'rb')
traj = pickle.load(dbfile)
dbfile.close()

dbfile = open(f"/Users/morgan.schneider/Documents/merger/supercell-125m/boxes_s1.pkl", 'rb')
box = pickle.load(dbfile)
xb1 = box['x1_pp']
xb2 = box['x2_pp']
yb1 = box['y1_pp']
yb2 = box['y2_pp']
dbfile.close()


fnums = [33, 43, 53]

ds = nc.Dataset('/Volumes/Promise_Pegasus_70TB/merger/supercell-125m/cm1out_pdata.nc')
ptime = ds.variables['time'][:].data
ds.close()


if use_SR:
    from scipy.interpolate import interp1d
    dbfile = open('/Users/morgan.schneider/Documents/merger/supercell-125m/storm_motion.pkl', 'rb')
    sm = pickle.load(dbfile)
    u_storm = sm['u_storm']
    v_storm = sm['v_storm']
    dbfile.close
    stimes = np.linspace(10800, 14400, 61)
    
    # Interpolate storm motion from model output timesto parcel output times (60 s to 15 s)
    fu = interp1d(stimes, u_storm)
    u_storm_interp = fu(ptime)
    fv = interp1d(stimes, v_storm)
    v_storm_interp = fv(ptime)
    
    xl = [-52,32]; yl = [-102,-18]
else:
    xl = [-62,22]; yl = [-122,-38]



ds = nc.Dataset('/Volumes/Promise_Pegasus_70TB/merger/supercell-125m/cm1out_000033.nc')
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
ix1 = np.where(xh >= xl[0])[0][0]; ix2 = np.where(xh >= xl[1])[0][1]; ix = slice(ix1,ix2)
iy1 = np.where(yh >= yl[0])[0][0]; iy2 = np.where(yh >= yl[1])[0][1]; iy = slice(iy1,iy2)
ds.close()




fig,axs = plt.subplots(3, 3, figsize=(9,8), subplot_kw=dict(box_aspect=1), sharex=True, sharey=True, layout='constrained')


for i in range(3):
    ds = nc.Dataset(f"/Volumes/Promise_Pegasus_70TB/merger/supercell-125m/cm1out_{fnums[i]:06d}.nc")
    stime = ds.variables['time'][0]
    t = stime/60
    dbz = ds.variables['dbz'][:].data[0,0,iy,ix]
    ds.close()
    
    dbfile = open(f"/Users/morgan.schneider/Documents/merger/supercell-125m/traj_clusters_{t:.0f}min_v2.pkl", 'rb')
    cc = pickle.load(dbfile)
    cc_mv = cc['s1']
    dbfile.close()
    
    ti = np.where(ptime == stime)[0][0]
    
    pids_mv = traj[f"{t:.0f}min"]['pids']
    x_mv = traj[f"{t:.0f}min"]['x'][:ti+1,:]
    y_mv = traj[f"{t:.0f}min"]['y'][:ti+1,:]
    z_mv = traj[f"{t:.0f}min"]['z'][:ti+1,:]
    w_mv = traj[f"{t:.0f}min"]['w'][:ti+1,:]
    b_mv = traj[f"{t:.0f}min"]['b'][:ti+1,:]

    
    
    # 3x3 panel trajectories colored by z, B, and source for qlcsonly
    x0 = x_mv[:,(cc_mv==0)]; y0 = y_mv[:,(cc_mv==0)]; z0 = z_mv[:,(cc_mv==0)]; b0 = b_mv[:,(cc_mv==0)]; c0 = cc_mv[(cc_mv==0)]
    x1 = x_mv[:,(cc_mv==1)]; y1 = y_mv[:,(cc_mv==1)]; z1 = z_mv[:,(cc_mv==1)]; b1 = b_mv[:,(cc_mv==1)]; c1 = cc_mv[(cc_mv==1)]
    x2 = x_mv[:,(cc_mv==2)]; y2 = y_mv[:,(cc_mv==2)]; z2 = z_mv[:,(cc_mv==2)]; b2 = b_mv[:,(cc_mv==2)]; c2 = cc_mv[(cc_mv==2)]
    
    if i == 0:
        xp = np.append(x0[:,::21], np.append(x1, x2, axis=1), axis=1)
        yp = np.append(y0[:,::21], np.append(y1, y2, axis=1), axis=1)
        zp = np.append(z0[:,::21], np.append(z1, z2, axis=1), axis=1)
        bp = np.append(b0[:,::21], np.append(b1, b2, axis=1), axis=1)
        cc = np.append(c0[::21], np.append(c1, c2))
    if i == 1:
        xp = np.append(x0[:,::7], np.append(x1, x2, axis=1), axis=1)
        yp = np.append(y0[:,::7], np.append(y1, y2, axis=1), axis=1)
        zp = np.append(z0[:,::7], np.append(z1, z2, axis=1), axis=1)
        bp = np.append(b0[:,::7], np.append(b1, b2, axis=1), axis=1)
        cc = np.append(c0[::7], np.append(c1, c2))
    if i == 2:
        xp = np.append(x0, np.append(x1, x2, axis=1), axis=1)
        yp = np.append(y0, np.append(y1, y2, axis=1), axis=1)
        zp = np.append(z0, np.append(z1, z2, axis=1), axis=1)
        bp = np.append(b0, np.append(b1, b2, axis=1), axis=1)
        cc = np.append(c0, np.append(c1, c2))
    
    if use_SR:
        xp,yp = get_SR_positions(xp, yp, u_storm_interp[:ti+1], v_storm_interp[:ti+1], 15)
    
    
    
    p1 = axs[0,i].scatter(xp/1000, yp/1000, s=0.5, c=zp/1000, marker='.', cmap='HomeyerRainbow', vmin=0, vmax=3)
    axs[0,i].contour(xh[ix], yh[iy], dbz, levels=[30], colors='k', linewidths=1)
    r1 = patches.Rectangle((xb1[10*i+19],yb1[10*i+19]), 10, 10, fc='none', ec='k', lw=2)
    axs[0,i].add_patch(r1)
    axs[0,i].set_xlim(xl)
    axs[0,i].set_ylim(yl)
    axs[0,i].set_title(f"{t:.0f} min", fontsize=16)
    
    p2 = axs[1,i].scatter(xp/1000, yp/1000, s=0.5, c=bp, marker='.', cmap='HomeyerRainbow', vmin=-0.3, vmax=0.1)
    axs[1,i].contour(xh[ix], yh[iy], dbz, levels=[30], colors='k', linewidths=1)
    r2 = patches.Rectangle((xb1[10*i+19],yb1[10*i+19]), 10, 10, fc='none', ec='k', lw=2)
    axs[1,i].add_patch(r2)
    axs[1,i].set_xlim(xl)
    axs[1,i].set_ylim(yl)
    
    p3 = axs[2,i].scatter([0,0], [0,0], s=1, c=[0,3], marker='.', cmap=cm, vmin=-0.5, vmax=2.5)
    axs[2,i].scatter(xp/1000, yp/1000, s=0.5, c=np.tile(cc,(len(xp[:,0]),1)), marker='.', cmap=cm, vmin=-0.5, vmax=2.5)
    axs[2,i].contour(xh[ix], yh[iy], dbz, levels=[30], colors='k', linewidths=1)
    r3 = patches.Rectangle((xb1[10*i+19],yb1[10*i+19]), 10, 10, fc='none', ec='k', lw=2)
    axs[2,i].add_patch(r3)
    axs[2,i].set_xlim(xl)
    axs[2,i].set_ylim(yl)
    axs[2,i].set_xlabel('x (km)', fontsize=14)

if False:
    l1, = axs[2,0].plot([xl[0]-2,xl[0]-1], [yl[0]-2,yl[0]-1], color='gold')
    l2, = axs[2,0].plot([xl[0]-2,xl[0]-1], [yl[0]-2,yl[0]-1], color='red')
    l3, = axs[2,0].plot([xl[0]-2,xl[0]-1], [yl[0]-2,yl[0]-1], color='deepskyblue')
    axs[2,0].legend(handles=[l1,l2,l3], labels=['Low-level inflow', 'Mid-level inflow', 'Supercell outflow'], loc='upper right', fontsize=8)

cb1 = plt.colorbar(p1, ax=axs[0,2], extend='both', ticks=np.linspace(0,3,7))
cb1.set_label("Height (km)", fontsize=12)
# tl = cb1.ax.get_yticklabels()
# cb1.ax.set_yticklabels(tl, ha='right')
# cb1.ax.yaxis.set_tick_params(pad=30)

cb2 = plt.colorbar(p2, ax=axs[1,2], extend='both', ticks=np.linspace(-0.3,0.1,5))
cb2.set_label("Buoyancy (m s$^{-2}$)", fontsize=12)
# tl = cb2.ax.get_yticklabels()
# cb2.ax.set_yticklabels(tl, ha='right')
# cb2.ax.yaxis.set_tick_params(pad=30)

cb3 = plt.colorbar(p3, ax=axs[2,2], ticks=[0,1,2])
cb3.ax.set_yticklabels([' Low-level\n inflow', ' Mid-level\n inflow', ' Supercell\n outflow'])
# cb3.set_label("Source region", fontsize=12)

axs[0,0].set_yticks(np.linspace(yl[0]+2, yl[1]-2, 5))
axs[1,0].set_yticks(np.linspace(yl[0]+2, yl[1]-2, 5))
axs[2,0].set_yticks(np.linspace(yl[0]+2, yl[1]-2, 5))

axs[0,0].set_ylabel('y (km)', fontsize=14)
axs[1,0].set_ylabel('y (km)', fontsize=14)
axs[2,0].set_ylabel('y (km)', fontsize=14)

axs[2,0].set_xlabel('x (km)', fontsize=14)
axs[2,1].set_xlabel('x (km)', fontsize=14)
axs[2,2].set_xlabel('x (km)', fontsize=14)

figsave = False

if figsave:
    fig.savefig('/Users/morgan.schneider/Documents/merger/traj_z+B+src_SUP_SR.png', dpi=300)



#%% One bigass time series plot of z/B/w/zvort for all sources (don't use)

from matplotlib.ticker import MultipleLocator


dbfile = open(f"/Users/morgan.schneider/Documents/merger/traj_MV1.pkl", 'rb')
traj_mv1 = pickle.load(dbfile)
dbfile.close()


fnums = [33, 43, 53]
times = [200, 210, 220]

ds = nc.Dataset('/Volumes/Promise_Pegasus_70TB/merger/merger-125m/cm1out_pdata.nc')
ptime = ds.variables['time'][:].data/60
ds.close()

n = np.where(ptime == 220)[0][0]



fig,ax = plt.subplots(4, 1, figsize=(15,11), sharex=True, layout='constrained')

# for a in range(len(ax)):
#     ax[a].axhline(0, color='darkgray', linestyle='-', linewidth=2)


for i in range(3):
    t = times[i]
    it1 = np.where(ptime >= t-10)[0][0]
    it2 = np.where(ptime > t)[0][0]
    it = slice(it1,it2)
    
    dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/traj_clusters/traj_clusters_{t:.0f}min_v2.pkl", 'rb')
    cc = pickle.load(dbfile)
    cc_mv1 = cc['mv1']
    if i == 2:
        cc_mv1[(cc_mv1 == 4)] = 0
    dbfile.close()
    
    
    ### MV1 ###
    c0_mv1 = (cc_mv1 == 0)
    c1_mv1 = (cc_mv1 == 1)
    c2_mv1 = (cc_mv1 == 2)
    c3_mv1 = (cc_mv1 == 3)
    
    z_mv1 = traj_mv1[f"{t:.0f}min"]['z']
    b_mv1 = traj_mv1[f"{t:.0f}min"]['b']
    w_mv1 = traj_mv1[f"{t:.0f}min"]['w']
    zvort_mv1 = traj_mv1[f"{t:.0f}min"]['zvort']
    
    # colorblind-friendly colors
    col1 = 'gold'
    col2 = 'red'
    col3 = 'deepskyblue'
    col4 = 'mediumblue'
    
    
    # QLCS cold pool
    if np.any(cc_mv1 == 3):
        c = c3_mv1
        ax[0].fill_between(ptime[it], np.percentile(z_mv1[it,c], 25, axis=1)/1000,
                             np.percentile(z_mv1[it,c], 75, axis=1)/1000, alpha=0.25, color=col4)
        ax[1].fill_between(ptime[it], np.percentile(b_mv1[it,c], 25, axis=1),
                             np.percentile(b_mv1[it,c], 75, axis=1), alpha=0.25, color=col4)
        ax[2].fill_between(ptime[it], np.percentile(w_mv1[it,c], 25, axis=1),
                             np.percentile(w_mv1[it,c], 75, axis=1), alpha=0.25, color=col4)
        ax[3].fill_between(ptime[it], np.percentile(zvort_mv1[it,c], 25, axis=1),
                             np.percentile(zvort_mv1[it,c], 75, axis=1), alpha=0.25, color=col4)
    
    # supercell cold pool
    if np.any(cc_mv1 == 2):
        c = c2_mv1
        ax[0].fill_between(ptime[it], np.percentile(z_mv1[it,c], 25, axis=1)/1000,
                             np.percentile(z_mv1[it,c], 75, axis=1)/1000, alpha=0.3, color=col3)
        ax[1].fill_between(ptime[it], np.percentile(b_mv1[it,c], 25, axis=1),
                             np.percentile(b_mv1[it,c], 75, axis=1), alpha=0.3, color=col3)
        ax[2].fill_between(ptime[it], np.percentile(w_mv1[it,c], 25, axis=1),
                             np.percentile(w_mv1[it,c], 75, axis=1), alpha=0.3, color=col3)
        ax[3].fill_between(ptime[it], np.percentile(zvort_mv1[it,c], 25, axis=1),
                             np.percentile(zvort_mv1[it,c], 75, axis=1), alpha=0.3, color=col3)
    
    # low-level inflow
    if np.any(cc_mv1 == 0):
        c = c0_mv1
        ax[0].fill_between(ptime[it], np.percentile(z_mv1[it,c], 25, axis=1)/1000,
                             np.percentile(z_mv1[it,c], 75, axis=1)/1000, alpha=0.3, color=col1)
        ax[1].fill_between(ptime[it], np.percentile(b_mv1[it,c], 25, axis=1),
                             np.percentile(b_mv1[it,c], 75, axis=1), alpha=0.3, color=col1)
        ax[2].fill_between(ptime[it], np.percentile(w_mv1[it,c], 25, axis=1),
                             np.percentile(w_mv1[it,c], 75, axis=1), alpha=0.3, color=col1)
        ax[3].fill_between(ptime[it], np.percentile(zvort_mv1[it,c], 25, axis=1),
                             np.percentile(zvort_mv1[it,c], 75, axis=1), alpha=0.3, color=col1)
    
    # mid-level inflow
    if np.any(cc_mv1 == 1):
        c = c1_mv1
        ax[0].fill_between(ptime[it], np.percentile(z_mv1[it,c], 25, axis=1)/1000,
                             np.percentile(z_mv1[it,c], 75, axis=1)/1000, alpha=0.2, color=col2)
        ax[1].fill_between(ptime[it], np.percentile(b_mv1[it,c], 25, axis=1),
                             np.percentile(b_mv1[it,c], 75, axis=1), alpha=0.2, color=col2)
        ax[2].fill_between(ptime[it], np.percentile(w_mv1[it,c], 25, axis=1),
                             np.percentile(w_mv1[it,c], 75, axis=1), alpha=0.2, color=col2)
        ax[3].fill_between(ptime[it], np.percentile(zvort_mv1[it,c], 25, axis=1),
                             np.percentile(zvort_mv1[it,c], 75, axis=1), alpha=0.2, color=col2)
    
    # QLCS cold pool
    s4, = ax[0].plot(ptime[it], np.median(z_mv1[it,c3_mv1],axis=1)/1000, col4, linewidth=2)
    ax[1].plot(ptime[it], np.median(b_mv1[it,c3_mv1],axis=1), col4, linewidth=2)
    ax[2].plot(ptime[it], np.median(w_mv1[it,c3_mv1],axis=1), col4, linewidth=2)
    ax[3].plot(ptime[it], np.median(zvort_mv1[it,c3_mv1],axis=1), col4, linewidth=2)
    # supercell cold pool
    s3, = ax[0].plot(ptime[it], np.median(z_mv1[it,c2_mv1],axis=1)/1000, col3, linewidth=2)
    col3 = 'dodgerblue'
    ax[0].plot(ptime[it], np.median(z_mv1[it,c2_mv1],axis=1)/1000, col3, linewidth=2)
    ax[1].plot(ptime[it], np.median(b_mv1[it,c2_mv1],axis=1), col3, linewidth=2)
    ax[2].plot(ptime[it], np.median(w_mv1[it,c2_mv1],axis=1), col3, linewidth=2)
    ax[3].plot(ptime[it], np.median(zvort_mv1[it,c2_mv1],axis=1), col3, linewidth=2)
    # low-level environment
    s1, = ax[0].plot(ptime[it], np.median(z_mv1[it,c0_mv1],axis=1)/1000, col1, linewidth=2)
    ax[1].plot(ptime[it], np.median(b_mv1[it,c0_mv1],axis=1), col1, linewidth=2)
    ax[2].plot(ptime[it], np.median(w_mv1[it,c0_mv1],axis=1), col1, linewidth=2)
    ax[3].plot(ptime[it], np.median(zvort_mv1[it,c0_mv1],axis=1), col1, linewidth=2)
    # mid-level environment
    s2, = ax[0].plot(ptime[it], np.median(z_mv1[it,c1_mv1],axis=1)/1000, col2, linewidth=2)
    # col2 = 'maroon'
    # ax[0].plot(ptime[it], np.median(z_mv1[it,c1_mv1],axis=1)/1000, col2, linewidth=2)
    ax[1].plot(ptime[it], np.median(b_mv1[it,c1_mv1],axis=1), col2, linewidth=2)
    ax[2].plot(ptime[it], np.median(w_mv1[it,c1_mv1],axis=1), col2, linewidth=2)
    ax[3].plot(ptime[it], np.median(zvort_mv1[it,c1_mv1],axis=1), col2, linewidth=2)
    

ax[0].legend(handles=[s1, s2, s3, s4],
               labels=['Low-level inflow', 'Mid-level inflow', 'Supercell outflow', 'QLCS outflow'],
               loc='upper left', bbox_to_anchor=(0.13,0.99), fontsize=14)

ax[3].set_xlabel('Time (min)', fontsize=24)
ax[3].tick_params(axis='x', which='major', labelsize=22)
ax[3].set_xlim([190,220])
ax[0].set_ylim([0,3.25])
ax[1].set_ylim([-0.4,0.1])
ax[2].set_ylim([-15,15])
ax[3].set_ylim([-0.02,0.05])
# ax[2].set_yticks([-15,-10,-5,0,5,10,15])
# ax[3].set_yticks([-0.04,-0.02,0,0.02,0.04])
ax[0].set_ylabel('Height AGL', fontsize=20)
ax[1].set_ylabel("Buoyancy", fontsize=20)
ax[2].set_ylabel("Vertical velocity", fontsize=20)
ax[3].set_ylabel("Vertical vorticity", fontsize=20)

for i in range(len(ax)):
    ax[i].axvline(200, color='k', linestyle='--', linewidth=3)
    ax[i].axvline(210, color='k', linestyle='--', linewidth=3)
    ax[i].axhline(0, color='k', linestyle='--', linewidth=1.5)
    ax[i].tick_params(axis='y', which='major', labelsize=15)
    ax[i].grid(visible=True, which='major', color='darkgray', linestyle='-')
    ax[i].grid(visible=True, which='minor', color='lightgray', linestyle='--')

ax[0].yaxis.set_major_locator(MultipleLocator(1))
ax[0].yaxis.set_minor_locator(MultipleLocator(0.5))
ax[1].yaxis.set_major_locator(MultipleLocator(0.1))
ax[2].yaxis.set_major_locator(MultipleLocator(5))
ax[3].yaxis.set_major_locator(MultipleLocator(0.02))
ax[3].yaxis.set_minor_locator(MultipleLocator(0.01))
ax[3].xaxis.set_major_locator(MultipleLocator(5))
ax[3].xaxis.set_minor_locator(MultipleLocator(1))

# was 3.2
ax[0].text(193.5, 3.5, '200 min', fontsize=24) # 185 for 195 min
ax[0].text(203.5, 3.5, '210 min', fontsize=24)
ax[0].text(213.5, 3.5, '220 min', fontsize=24) # 215 for 225 min

fig.align_ylabels(ax[:])

# plt.show()

figsave = False

if figsave:
    fig.savefig('/Users/morgan.schneider/Documents/merger/traj_timeseries_old.png', dpi=300)


#%% Bigass time series v2 ***ACTUAL PAPER FIG***

from matplotlib.ticker import MultipleLocator


dbfile = open(f"/Users/morgan.schneider/Documents/merger/traj_MV1.pkl", 'rb')
traj_mv1 = pickle.load(dbfile)
dbfile.close()


fnums = [33, 43, 53]
times = [200, 210, 220]

ds = nc.Dataset('/Volumes/Promise_Pegasus_70TB/merger/merger-125m/cm1out_pdata.nc')
ptime = ds.variables['time'][:].data/60
ds.close()

n = np.where(ptime == 220)[0][0]

# colorblind-friendly colors
col1 = 'gold'
col2 = 'red'
col3 = 'deepskyblue'
col4 = 'mediumblue'



fig,ax = plt.subplots(4, 3, figsize=(16.5,12))

for i in range(3):
    t = times[i]
    it1 = np.where(ptime >= t-15)[0][0]
    it2 = np.where(ptime > t)[0][0]
    it = slice(it1,it2)
    
    dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/traj_clusters/traj_clusters_{t:.0f}min_v2.pkl", 'rb')
    cc = pickle.load(dbfile)
    cc_mv1 = cc['mv1']
    if i == 2:
        cc_mv1[(cc_mv1 == 4)] = 0
    dbfile.close()
    
    
    ### MV1 ###
    c0_mv1 = (cc_mv1 == 0)
    c1_mv1 = (cc_mv1 == 1)
    c2_mv1 = (cc_mv1 == 2)
    c3_mv1 = (cc_mv1 == 3)
    
    z_mv1 = traj_mv1[f"{t:.0f}min"]['z']
    b_mv1 = traj_mv1[f"{t:.0f}min"]['b']
    w_mv1 = traj_mv1[f"{t:.0f}min"]['w']
    zvort_mv1 = traj_mv1[f"{t:.0f}min"]['zvort']
    
    
    # QLCS cold pool
    if np.any(cc_mv1 == 3):
        c = c3_mv1
        ax[0,i].fill_between(ptime[it], np.percentile(z_mv1[it,c], 25, axis=1)/1000,
                             np.percentile(z_mv1[it,c], 75, axis=1)/1000, alpha=0.25, color=col4)
        ax[1,i].fill_between(ptime[it], np.percentile(b_mv1[it,c], 25, axis=1),
                             np.percentile(b_mv1[it,c], 75, axis=1), alpha=0.25, color=col4)
        ax[2,i].fill_between(ptime[it], np.percentile(w_mv1[it,c], 25, axis=1),
                             np.percentile(w_mv1[it,c], 75, axis=1), alpha=0.25, color=col4)
        ax[3,i].fill_between(ptime[it], np.percentile(zvort_mv1[it,c], 25, axis=1),
                             np.percentile(zvort_mv1[it,c], 75, axis=1), alpha=0.25, color=col4)
    
    # supercell cold pool
    if np.any(cc_mv1 == 2):
        c = c2_mv1
        ax[0,i].fill_between(ptime[it], np.percentile(z_mv1[it,c], 25, axis=1)/1000,
                             np.percentile(z_mv1[it,c], 75, axis=1)/1000, alpha=0.3, color=col3)
        ax[1,i].fill_between(ptime[it], np.percentile(b_mv1[it,c], 25, axis=1),
                             np.percentile(b_mv1[it,c], 75, axis=1), alpha=0.3, color=col3)
        ax[2,i].fill_between(ptime[it], np.percentile(w_mv1[it,c], 25, axis=1),
                             np.percentile(w_mv1[it,c], 75, axis=1), alpha=0.3, color=col3)
        ax[3,i].fill_between(ptime[it], np.percentile(zvort_mv1[it,c], 25, axis=1),
                             np.percentile(zvort_mv1[it,c], 75, axis=1), alpha=0.3, color=col3)
    
    # low-level inflow
    if np.any(cc_mv1 == 0):
        c = c0_mv1
        ax[0,i].fill_between(ptime[it], np.percentile(z_mv1[it,c], 25, axis=1)/1000,
                             np.percentile(z_mv1[it,c], 75, axis=1)/1000, alpha=0.3, color=col1)
        ax[1,i].fill_between(ptime[it], np.percentile(b_mv1[it,c], 25, axis=1),
                             np.percentile(b_mv1[it,c], 75, axis=1), alpha=0.3, color=col1)
        ax[2,i].fill_between(ptime[it], np.percentile(w_mv1[it,c], 25, axis=1),
                             np.percentile(w_mv1[it,c], 75, axis=1), alpha=0.3, color=col1)
        ax[3,i].fill_between(ptime[it], np.percentile(zvort_mv1[it,c], 25, axis=1),
                             np.percentile(zvort_mv1[it,c], 75, axis=1), alpha=0.3, color=col1)
    
    # mid-level inflow
    if np.any(cc_mv1 == 1):
        c = c1_mv1
        ax[0,i].fill_between(ptime[it], np.percentile(z_mv1[it,c], 25, axis=1)/1000,
                             np.percentile(z_mv1[it,c], 75, axis=1)/1000, alpha=0.2, color=col2)
        ax[1,i].fill_between(ptime[it], np.percentile(b_mv1[it,c], 25, axis=1),
                             np.percentile(b_mv1[it,c], 75, axis=1), alpha=0.2, color=col2)
        ax[2,i].fill_between(ptime[it], np.percentile(w_mv1[it,c], 25, axis=1),
                             np.percentile(w_mv1[it,c], 75, axis=1), alpha=0.2, color=col2)
        ax[3,i].fill_between(ptime[it], np.percentile(zvort_mv1[it,c], 25, axis=1),
                             np.percentile(zvort_mv1[it,c], 75, axis=1), alpha=0.2, color=col2)
    
    # QLCS cold pool
    s4, = ax[0,i].plot(ptime[it], np.median(z_mv1[it,c3_mv1],axis=1)/1000, col4, linewidth=2)
    ax[1,i].plot(ptime[it], np.median(b_mv1[it,c3_mv1],axis=1), col4, linewidth=2)
    ax[2,i].plot(ptime[it], np.median(w_mv1[it,c3_mv1],axis=1), col4, linewidth=2)
    ax[3,i].plot(ptime[it], np.median(zvort_mv1[it,c3_mv1],axis=1), col4, linewidth=2)
    # supercell cold pool
    s3, = ax[0,i].plot(ptime[it], np.median(z_mv1[it,c2_mv1],axis=1)/1000, col3, linewidth=2)
    col3 = 'dodgerblue'
    ax[0,i].plot(ptime[it], np.median(z_mv1[it,c2_mv1],axis=1)/1000, col3, linewidth=2)
    ax[1,i].plot(ptime[it], np.median(b_mv1[it,c2_mv1],axis=1), col3, linewidth=2)
    ax[2,i].plot(ptime[it], np.median(w_mv1[it,c2_mv1],axis=1), col3, linewidth=2)
    ax[3,i].plot(ptime[it], np.median(zvort_mv1[it,c2_mv1],axis=1), col3, linewidth=2)
    # low-level environment
    s1, = ax[0,i].plot(ptime[it], np.median(z_mv1[it,c0_mv1],axis=1)/1000, col1, linewidth=2)
    ax[1,i].plot(ptime[it], np.median(b_mv1[it,c0_mv1],axis=1), col1, linewidth=2)
    ax[2,i].plot(ptime[it], np.median(w_mv1[it,c0_mv1],axis=1), col1, linewidth=2)
    ax[3,i].plot(ptime[it], np.median(zvort_mv1[it,c0_mv1],axis=1), col1, linewidth=2)
    # mid-level environment
    s2, = ax[0,i].plot(ptime[it], np.median(z_mv1[it,c1_mv1],axis=1)/1000, col2, linewidth=2)
    # col2 = 'maroon'
    # ax[0,i].plot(ptime[it], np.median(z_mv1[it,c1_mv1],axis=1)/1000, col2, linewidth=2)
    ax[1,i].plot(ptime[it], np.median(b_mv1[it,c1_mv1],axis=1), col2, linewidth=2)
    ax[2,i].plot(ptime[it], np.median(w_mv1[it,c1_mv1],axis=1), col2, linewidth=2)
    ax[3,i].plot(ptime[it], np.median(zvort_mv1[it,c1_mv1],axis=1), col2, linewidth=2)
    

ax[0,0].legend(handles=[s1, s2, s3, s4],
               labels=['Low-level inflow', 'Mid-level inflow', 'Supercell outflow', 'QLCS outflow'],
               loc='upper right', fontsize=14)

ax[3,1].set_xlabel('Time (min)', fontsize=20)

ax[0,0].set_ylabel('Height AGL', fontsize=18)
ax[1,0].set_ylabel("Buoyancy", fontsize=18)
ax[2,0].set_ylabel("Vertical velocity", fontsize=18)
ax[3,0].set_ylabel("Vertical vorticity", fontsize=18)

for i in range(3):
    ax[0,i].set_ylim([0,3.25])
    ax[1,i].set_ylim([-0.4,0.1])
    ax[2,i].set_ylim([-15,15])
    ax[3,i].set_ylim([-0.02,0.05])
    ax[3,i].tick_params(axis='x', which='major', labelsize=16)
    
    ax[0,i].yaxis.set_major_locator(MultipleLocator(1))
    ax[0,i].yaxis.set_minor_locator(MultipleLocator(0.5))
    ax[1,i].yaxis.set_major_locator(MultipleLocator(0.1))
    ax[2,i].yaxis.set_major_locator(MultipleLocator(5))
    ax[3,i].yaxis.set_major_locator(MultipleLocator(0.02))
    ax[3,i].yaxis.set_minor_locator(MultipleLocator(0.01))
    for j in range(4):
        ax[j,i].xaxis.set_major_locator(MultipleLocator(5))
        ax[j,i].xaxis.set_minor_locator(MultipleLocator(1))
        ax[j,i].axhline(0, color='k', linestyle='--', linewidth=1.5)
        ax[j,i].tick_params(axis='y', which='major', labelsize=15)
        ax[j,i].grid(visible=True, which='major', color='darkgray', linestyle='-')
        ax[j,i].grid(visible=True, which='minor', color='lightgray', linestyle='--')
        if i > 0:
            ax[j,i].set_yticklabels([])
            # ax[j,i].axvline(ptime[it][0], color='k', linestyle='--', linewidth=3)

for j in range(4):
    ax[j,0].set_xlim([185,200])
    ax[j,1].set_xlim([195,210])
    ax[j,2].set_xlim([205,220])
    if j < 3:
        for i in range(3):
            ax[j,i].set_xticklabels([])

ax[3,0].set_xticks([185, 190, 195, 200])
ax[3,0].set_xticklabels(["    185", "190", "195", "200   "])

ax[3,1].set_xticks([195, 200, 205, 210])
ax[3,1].set_xticklabels(["    195", "200", "205", "210   "])

ax[3,2].set_xticks([205, 210, 215, 220])
ax[3,2].set_xticklabels(["    205", "210", "215", "220"])

# for j in range(3):
#     ticks = ax[3,j].xaxis.get_majorticklabels()
#     for i in range(len(ticks)):
#         if i == 0:
#             ticks[i].set_horizontalalignment('left')
        # elif i == len(ticks)-1:
        #     ticks[i].set_horizontalalignment('right')
        # else:
        #     ticks[i].set_horizontalalignment('center')

# ax[2].set_yticks([-15,-10,-5,0,5,10,15])
# ax[3].set_yticks([-0.04,-0.02,0,0.02,0.04])

ax[0,0].text(190, 3.5, '200 min', fontsize=22) # 185 for 195 min
ax[0,1].text(200, 3.5, '210 min', fontsize=22)
ax[0,2].text(210, 3.5, '220 min', fontsize=22) # 215 for 225 min

fig.align_ylabels(ax[:])

plt.subplots_adjust(wspace=0.06, hspace=0.15)
# plt.show()

figsave = False

if figsave:
    fig.savefig('/Users/morgan.schneider/Documents/merger/traj_timeseries.png', dpi=300)


#%% zoomed-in time series of w and zeta for mid-level

ds = nc.Dataset('/Volumes/Promise_Pegasus_70TB/merger/merger-125m/cm1out_pdata.nc')
ptime = ds.variables['time'][:].data/60
ds.close()

t1 = 210
span = 15
itmp1 = np.where(ptime == t1-span)[0][0]
itmp2 = np.where(ptime > t1)[0][0]
it1 = slice(itmp1,itmp2)

t2 = 220
itmp1 = np.where(ptime == t2-span)[0][0]
itmp2 = np.where(ptime > t2)[0][0]
it2 = slice(itmp1,itmp2)

dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/traj_clusters/traj_clusters_{t1:.0f}min_v2.pkl", 'rb')
tmp = pickle.load(dbfile)
cc1 = tmp['mv1']
dbfile.close()

dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/traj_clusters/traj_clusters_{t2:.0f}min_v2.pkl", 'rb')
tmp = pickle.load(dbfile)
cc2 = tmp['mv1']
dbfile.close()


w_ml1 = traj_mv1[f"{t1:.0f}min"]['w'][it1,(cc1==1)]
zvort_ml1 = traj_mv1[f"{t1:.0f}min"]['zvort'][it1,(cc1==1)]
w_ml2 = traj_mv1[f"{t2:.0f}min"]['w'][it2,(cc2==1)]
zvort_ml2 = traj_mv1[f"{t2:.0f}min"]['zvort'][it2,(cc2==1)]


# colorblind-friendly colors
col1 = 'red'
col2 = 'gold'
col3 = 'deepskyblue'
col4 = 'mediumblue'

fig,axs = plt.subplots(2, 2, figsize=(10,6))

axs[0,0].fill_between(ptime[it1], np.percentile(w_ml1, 25, axis=1),
                     np.percentile(w_ml1, 75, axis=1), alpha=0.3, color=col3)
axs[1,0].fill_between(ptime[it1], np.percentile(zvort_ml1, 25, axis=1),
                     np.percentile(zvort_ml1, 75, axis=1), alpha=0.3, color=col1)
axs[0,1].fill_between(ptime[it2], np.percentile(w_ml2, 25, axis=1),
                     np.percentile(w_ml2, 75, axis=1), alpha=0.3, color=col3)
axs[1,1].fill_between(ptime[it2], np.percentile(zvort_ml2, 25, axis=1),
                     np.percentile(zvort_ml2, 75, axis=1), alpha=0.3, color=col1)
axs[0,0].plot(ptime[it1], np.median(w_ml1,axis=1), col3, linewidth=2)
axs[1,0].plot(ptime[it1], np.median(zvort_ml1,axis=1), col1, linewidth=2)
axs[0,1].plot(ptime[it2], np.median(w_ml2,axis=1), col3, linewidth=2)
axs[1,1].plot(ptime[it2], np.median(zvort_ml2,axis=1), col1, linewidth=2)



for i in range(2): # column loop
    for j in range(2): # row loop
        if span == 5:
            axs[j,i].xaxis.set_major_locator(MultipleLocator(1))
        elif span == 10:
            axs[j,i].xaxis.set_major_locator(MultipleLocator(2))
            axs[j,i].xaxis.set_minor_locator(MultipleLocator(1))
        elif span > 10:
            axs[j,i].xaxis.set_major_locator(MultipleLocator(5))
            axs[j,i].xaxis.set_minor_locator(MultipleLocator(1))
        
        if i == 0: # if left column (210 min)
            axs[j,i].set_xlim([t1-span,t1])
        elif i == 1: # if right column (220 min)
            axs[j,i].set_xlim([t2-span,t2])
            axs[j,i].set_yticklabels([])
        
        if j == 0: # if top row (w)
            axs[j,i].set_xticklabels([])
            axs[j,i].set_ylim([-16,12])
            axs[j,i].yaxis.set_major_locator(MultipleLocator(4))
        elif j == 1: # if bottom row (zeta)
            axs[j,i].set_ylim([-0.01,0.05])
            axs[j,i].tick_params(axis='x', which='major', labelsize=12)
            axs[j,i].yaxis.set_major_locator(MultipleLocator(0.01))
        
        axs[j,i].axhline(0, color='k', linestyle='--', linewidth=1.5)
        axs[j,i].tick_params(axis='y', which='major', labelsize=12)
        axs[j,i].grid(visible=True, which='major', color='darkgray', linestyle='-')
        axs[j,i].grid(visible=True, which='minor', color='lightgray', linestyle='--')

axs[0,0].set_ylabel('Vertical velocity', fontsize=15)
axs[1,0].set_ylabel('Vertical vorticity', fontsize=15)
axs[1,0].set_xlabel('                                                        Time (min)', fontsize=16)

if span == 5:
    axs[1,0].set_xticks([205,206,207,208,209,210])
    axs[1,1].set_xticks([215,216,217,218,219,220])
    axs[1,0].set_xticklabels(["   205", "206", "207", "208", "209", "210    "])
    axs[1,1].set_xticklabels(["    215", "216", "217", "218", "219", "220"])
elif span == 10:
    axs[1,0].set_xticks([200,202,204,206,208,210])
    axs[1,1].set_xticks([210,212,214,216,218,220])
    axs[1,0].set_xticklabels(["   200", "202", "204", "206", "208", "210    "])
    axs[1,1].set_xticklabels(["    210", "212", "214", "216", "218", "220"])
elif span == 15:
    axs[1,0].set_xticks([195,200,205,210])
    axs[1,1].set_xticks([205,210,215,220])
    axs[1,0].set_xticklabels(["   195", "200", "205", "210    "])
    axs[1,1].set_xticklabels(["    205", "210", "215", "220"])

axs[0,0].text(201, 13, '210 min', fontsize=16) # 185 for 195 min
axs[0,1].text(211, 13, '220 min', fontsize=16)

fig.align_ylabels(axs[:,0])

plt.subplots_adjust(wspace=0.05, hspace=0.1)
# plt.show()

figsave = False

if figsave:
    fig.savefig('/Users/morgan.schneider/Documents/merger/traj_mv1_timeseries_zoomed.png', dpi=300)
    


#%% Load stuff for 3D plots of MERGER MV1 trajectories

fnum = 43

xl = [-40,20]
yl = [-110,-30]

ds = nc.Dataset(f"/Volumes/Promise_Pegasus_70TB/merger/merger-125m/cm1out_{fnum:06d}.nc")
stime = ds.variables['time'][:].data[0]
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
z = ds.variables['z'][:].data
ix1 = np.where(xh >= xl[0])[0][0]
ix2 = np.where(xh >= xl[1])[0][0]
iy1 = np.where(yh >= yl[0])[0][0]
iy2 = np.where(yh >= yl[1])[0][0]
ix = slice(ix1,ix2+1)
iy = slice(iy1,iy2+1)
xx,yy = np.meshgrid(xh[ix], yh[iy], indexing='xy')
# x3,y3,z3 = np.meshgrid(xh[ix], yh[iy], z, indexing='xy')
dbz = ds.variables['dbz'][:].data[0,0,iy,ix]
ds.close()

dbfile = open(f"/Users/morgan.schneider/Documents/merger/traj_MV1.pkl", 'rb')
traj_mv1 = pickle.load(dbfile)
dbfile.close()

dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/traj_clusters/traj_clusters_{stime/60:.0f}min_v2.pkl", 'rb')
cc = pickle.load(dbfile)
cc_mv1 = cc['mv1']
dbfile.close()

ds = nc.Dataset('/Volumes/Promise_Pegasus_70TB/merger/merger-125m/cm1out_pdata.nc')
ptime = ds.variables['time'][:].data
ds.close()

ti = np.where(ptime == stime)[0][0]
ti0 = np.where(ptime == stime-1800)[0][0]

x_mv1 = traj_mv1[f"{stime/60:.0f}min"]['x'][:,(cc_mv1 <= 3)]
y_mv1 = traj_mv1[f"{stime/60:.0f}min"]['y'][:,(cc_mv1 <= 3)]
z_mv1 = traj_mv1[f"{stime/60:.0f}min"]['z'][:,(cc_mv1 <= 3)]
# b_mv1 = traj_mv1[f"{stime/60:.0f}min"]['b'][:,(cc_mv1 <= 3)]
w_mv1 = traj_mv1[f"{stime/60:.0f}min"]['w'][:,(cc_mv1 <= 3)]
zvort_mv1 = traj_mv1[f"{stime/60:.0f}min"]['zvort'][:,(cc_mv1 <= 3)]
# vpg_mv1 = traj_mv1[f"{stime/60:.0f}min"]['vpg'][:,(cc_mv1 <= 3)]
cc_mv1 = cc_mv1[(cc_mv1 <= 3)]

c0_mv1 = (cc_mv1 == 0)
c1_mv1 = (cc_mv1 == 1)
c2_mv1 = (cc_mv1 == 2)
c3_mv1 = (cc_mv1 == 3)


x_ml = x_mv1[:,c1_mv1]
y_ml = y_mv1[:,c1_mv1]
z_ml = z_mv1[:,c1_mv1]
w_ml = w_mv1[:,c1_mv1]
zvort_ml = zvort_mv1[:,c1_mv1]

dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/hvort_traj_{stime/60:.0f}min.pkl", 'rb')
vort_traj = pickle.load(dbfile)
xvort_ml = vort_traj['xvort_ml']
yvort_ml = vort_traj['yvort_ml']
hvort_ml = vort_traj['hvort_ml']
vort_ml = vort_traj['vort_ml']
vort_sw = vort_traj['vort_sw_ml']
vort_cw = vort_traj['vort_cw_ml_signed']
dbfile.close()


import matplotlib as mpl
cm = mpl.colors.ListedColormap(['red', 'gold', 'deepskyblue', 'mediumblue'])

#%% Make 3D trajectory plots for MERGER MV1 colored by whatever variables

xlims = [-35,15]
ylims = [-82,-32]
zlims = [0,3]

ixx = slice(np.where(xh[ix]>=xlims[0])[0][0], np.where(xh[ix]>=xlims[1])[0][1])
iyy = slice(np.where(yh[iy]>=ylims[0])[0][0], np.where(yh[iy]>=ylims[1])[0][1])

view_angle = 45 # azim viewing angles: -5 for E, -30 for SE, -120 for SW
if view_angle == -5:
    view_dir = 'E'
elif view_angle == -30:
    view_dir = 'SE'
elif view_angle == -120:
    view_dir = 'SW'

qix = 4

ip = '/Users/morgan.schneider/Documents/merger/merger-125m/imgs_trajectories/'

figsave = False

if False:
    fig = plt.figure(figsize=(9,6))
    ax1 = fig.add_subplot(1,1,1, projection='3d')
    ax1.contour(xx[iyy,ixx], yy[iyy,ixx], dbz[iyy,ixx], levels=[30], zdir='z', offset=0, colors='k', linewidths=1)
    ax1.scatter(x_mv1[:,::qix]/1000, y_mv1[:,::qix]/1000, z_mv1[:,::qix]/1000, s=30, c=np.tile(cc_mv1[::qix],(len(x_mv1[:,0]),1)), cmap=cm, vmin=0, vmax=3)
    ax1.set_xlim(xlims)
    ax1.set_ylim(ylims)
    ax1.set_zlim(zlims)
    ax1.set_xlabel('x (km)', fontsize=14)
    ax1.set_ylabel('y (km)', fontsize=14)
    # ax1.set_zlabel('z (km)', fontsize=14)
    ax1.set_title(f"Sources, {stime/60:.0f} min \n view from {view_dir}", fontsize=14)
    ax1.view_init(elev=15, azim=view_angle) # azim viewing angles: -5 for E, -30 for SE, -120 for SW
    plt.tight_layout()
    if figsave:
        plt.savefig(ip+f"MV1_sources3d_view{view_angle:.0f}_{stime/60:.0f}min.png", dpi=300)


x_median = np.median(x_ml[:ti+1,:], axis=1)
y_median = np.median(y_ml[:ti+1,:], axis=1)
z_median = np.median(z_ml[:ti+1,:], axis=1)
xvort_median = np.median(xvort_ml[:ti+1,:], axis=1)
yvort_median = np.median(yvort_ml[:ti+1,:], axis=1)
zvort_median = np.median(zvort_ml[:ti+1,:], axis=1)


if True:
    from mpl_toolkits.mplot3d import Axes3D
    qit = 20
    
    fig = plt.figure(figsize=(9,6))
    ax1 = fig.add_subplot(1,1,1, projection='3d')
    ax1.contour(xx[iyy,ixx], yy[iyy,ixx], dbz[iyy,ixx], levels=[30], zdir='z', offset=0, colors='k', linewidths=1)
    # ax1.scatter(x_ml[:ti+1,:]/1000, y_ml[:ti+1,:]/1000, z_ml[:ti+1,:]/1000, s=10, marker='.', c=hvort_ml[:ti+1,:], cmap='HomeyerRainbow', vmin=0, vmax=0.08)
    ax1.scatter(x_median[:ti+1]/1000, y_median[:ti+1]/1000, z_median[:ti+1]/1000, s=10, marker='.', c=np.median(vort_ml[:ti+1,:],axis=1), cmap='HomeyerRainbow', vmin=0, vmax=0.1)
    ax1.quiver(x_median[::qit]/1000, y_median[::qit]/1000, z_median[::qit]/1000, xvort_median[::qit], yvort_median[::qit], 0*zvort_median[::qit], color='b', length=100*np.mean(hvort_ml[:ti+1,:]/vort_ml[:ti+1,:]), arrow_length_ratio=0.1, normalize=False)
    ax1.quiver(x_median[::qit]/1000, y_median[::qit]/1000, z_median[::qit]/1000, 0*xvort_median[::qit], 0*yvort_median[::qit], zvort_median[::qit], color='r', length=100*np.mean(np.abs(zvort_ml[:ti+1,:])/vort_ml[:ti+1,:]), arrow_length_ratio=0.1, normalize=False)
    ax1.set_xlim(xlims)
    ax1.set_ylim(ylims)
    ax1.set_zlim(zlims)
    ax1.set_xlabel('x (km)', fontsize=14)
    ax1.set_ylabel('y (km)', fontsize=14)
    # ax1.set_zlabel('z (km)', fontsize=14)
    ax1.set_title(f"Sources, {stime/60:.0f} min \n view from {view_dir}", fontsize=14)
    ax1.view_init(elev=25, azim=view_angle) # azim viewing angles: -5 for E, -30 for SE, -120 for SW
    plt.tight_layout()
    if figsave:
        plt.savefig(ip+f"traj_hvort3d_view{view_angle:.0f}_{stime/60:.0f}min.png", dpi=300)
    
    qit = 8
    
    fig = plt.figure(figsize=(9,6))
    ax1 = fig.add_subplot(1,1,1, projection='3d')
    ax1.contour(xx[iyy,ixx], yy[iyy,ixx], dbz[iyy,ixx], levels=[30], zdir='z', offset=0, colors='k', linewidths=1)
    # ax1.scatter(x_ml[:ti+1,:]/1000, y_ml[:ti+1,:]/1000, z_ml[:ti+1,:]/1000, s=10, marker='.', c=hvort_ml[:ti+1,:], cmap='HomeyerRainbow', vmin=0, vmax=0.08)
    ax1.scatter(x_median/1000, y_median/1000, z_median/1000, s=10, marker='.', c=np.median(vort_ml[:ti+1,:],axis=1), cmap='HomeyerRainbow', vmin=0, vmax=0.1)
    ax1.quiver(x_median[::qit]/1000, y_median[::qit]/1000, z_median[::qit]/1000, xvort_median[::qit], yvort_median[::qit], zvort_median[::qit], color='k', linewidth=2, length=1, arrow_length_ratio=0.2, normalize=True)
    ax1.set_xlim(xlims)
    ax1.set_ylim(ylims)
    ax1.set_zlim(zlims)
    ax1.set_xlabel('x (km)', fontsize=14)
    ax1.set_ylabel('y (km)', fontsize=14)
    # ax1.set_zlabel('z (km)', fontsize=14)
    ax1.set_title(f"Sources, {stime/60:.0f} min \n view from {view_dir}", fontsize=14)
    ax1.view_init(elev=25, azim=view_angle) # azim viewing angles: -5 for E, -30 for SE, -120 for SW
    plt.tight_layout()
    if figsave:
        plt.savefig(ip+f"traj_hvort3d_view{view_angle:.0f}_{stime/60:.0f}min.png", dpi=300)

#%% Define shit for 3D vorticity vectors (ripped from Jiang & Dawson 2023)

# https://github.itap.purdue.edu/STorMLab/Supercell_surface_boundaries/tree/master/QinScript/JupyterNotebookForSSBPaper2022
# _Fig13_TornadogenesisParcels_DistributionVort

from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d.axes3d import Axes3D
from mpl_toolkits.mplot3d.proj3d import proj_transform

class Arrow3D(FancyArrowPatch):
    def __init__(self, x, y, z, dx, dy, dz, *args, **kwargs):
        super().__init__((0, 0), (0, 0), *args, **kwargs)
        self._xyz = (x, y, z)
        self._dxdydz = (dx, dy, dz)
    
    def draw(self, renderer):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)

        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        super().draw(renderer)
    
    def do_3d_projection(self, renderer=None):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)

        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        return np.min(zs)
        
def _arrow3D(ax, x, y, z, dx, dy, dz ,lengthscale=10, *args, **kwargs):
    '''Add an 3d arrow to an `Axes3D` instance.'''
    for i in range(x.shape[0]):
        dx1 = dx[i]*lengthscale
        dy1 = dy[i]*lengthscale
        dz1 = dz[i]*lengthscale
        arrow = Arrow3D(x[i], y[i], z[i], dx1, dy1, dz1, *args, **kwargs)
        ax.add_artist(arrow)

setattr(Axes3D, 'arrow3D', _arrow3D)

#%% Plot the 3D trajectories and vorticity vectors (Jiang & Dawson 2023)

if stime/60 == 210:
    xlims = [-25,15]
    ylims = [-95,-55]
    zlims = [0,3]
elif stime/60 == 225:
    xlims = [-20,20]
    ylims = [-80,-40]
    zlims = [0,3]

ixx = slice(np.where(xh[ix]>=xlims[0])[0][0], np.where(xh[ix]>=xlims[1])[0][0])
iyy = slice(np.where(yh[iy]>=ylims[0])[0][0], np.where(yh[iy]>=ylims[1])[0][0])

x_median = np.median(x_ml[ti0:ti+1,:], axis=1)
y_median = np.median(y_ml[ti0:ti+1,:], axis=1)
z_median = np.median(z_ml[ti0:ti+1,:], axis=1)
xvort_median = np.median(xvort_ml[ti0:ti+1,:], axis=1)
yvort_median = np.median(yvort_ml[ti0:ti+1,:], axis=1)
zvort_median = np.median(zvort_ml[ti0:ti+1,:], axis=1)

qit = 2

view_angle = 25 # 150, 165, 190

fig = plt.figure(figsize=(9,6))
ax = fig.add_subplot(projection='3d')
ax.contour(xx[iyy,ixx], yy[iyy,ixx], dbz[iyy,ixx], levels=[30], zdir='z', offset=0, colors='gray', linewidths=1)
# p = ax.scatter3D(x_median/1000, y_median/1000, z_median/1000, s=5, marker='.',
#                  c=np.median(vort_ml[ti0:ti+1,:],axis=1), cmap='HomeyerRainbow', vmin=0, vmax=0.1)
p = ax.scatter3D(x_median/1000, y_median/1000, z_median/1000, s=5, marker='.', color='dimgray')
# a = ax.arrow3D(x_median[::qit].reshape(-1)/1000, y_median[::qit].reshape(-1)/1000, z_median[::qit].reshape(-1)/1000,
#            xvort_median[::qit].reshape(-1)*13, yvort_median[::qit].reshape(-1)*13, zvort_median[::qit].reshape(-1),
#             lengthscale=20, mutation_scale=7, edgecolor='k', facecolor='r', linewidth=0.5)
a1 = ax.arrow3D(x_median[::qit].reshape(-1)/1000, y_median[::qit].reshape(-1)/1000, z_median[::qit].reshape(-1)/1000,
           xvort_median[::qit].reshape(-1)*13, yvort_median[::qit].reshape(-1)*13, 0*zvort_median[::qit].reshape(-1),
            lengthscale=25, mutation_scale=7, edgecolor='k', facecolor='b', linewidth=0.5)
a2 = ax.arrow3D(x_median[::qit].reshape(-1)/1000, y_median[::qit].reshape(-1)/1000, z_median[::qit].reshape(-1)/1000,
           xvort_median[::qit].reshape(-1)*0, yvort_median[::qit].reshape(-1)*0, zvort_median[::qit].reshape(-1),
            lengthscale=25, mutation_scale=7, edgecolor='k', facecolor='r', linewidth=0.5)
ax.view_init(elev=5, azim=view_angle)
ax.xaxis.pane.fill = False
ax.yaxis.pane.fill = False
ax.zaxis.pane.fill = False
ax.xaxis._axinfo["grid"].update({"linewidth":0.5, "color" : "gray","linestyle": '--'})
ax.yaxis._axinfo["grid"].update({"linewidth":0.5, "color" : "gray","linestyle": '--'})
ax.zaxis._axinfo["grid"].update({"linewidth":0.5, "color" : "gray","linestyle": '--'})
ax.set_xlabel('x (km)',fontsize=12)
ax.set_ylabel('y (km)',fontsize=12)
ax.zaxis.set_rotate_label(False) 
ax.set_zlabel('Height (km)',fontsize=12,labelpad=5,rotation=90)
ax.set_xlim(xlims)
ax.set_ylim(ylims)
ax.set_zlim(zlims)
ax.set_title(f"Vorticity vectors, {stime/60:.0f} min (view from ENE)", fontsize=12)

figsave = False
if figsave:
    plt.savefig(f"/Users/morgan.schneider/Documents/merger/traj_vort3d_view{view_angle:.0f}_{stime/60:.0f}min.png", dpi=300)
    # plt.savefig(f"/Users/morgan.schneider/Documents/merger/traj_vort3d_sideview_{stime/60:.0f}min.png", dpi=300)


# fig = plt.figure(figsize=(9,6))
# ax1 = fig.add_subplot(1,1,1, projection='3d')
# ax1.contour(xx[iyy,ixx], yy[iyy,ixx], dbz[iyy,ixx], levels=[30], zdir='z', offset=0, colors='k', linewidths=1)
# ax1.scatter(x_mv1[:,c1_mv1]/1000, y_mv1[:,c1_mv1]/1000, z_mv1[:,c1_mv1]/1000, s=1, c=z_mv1[:,c1_mv1]/1000, marker='.', cmap='HomeyerRainbow', vmin=0, vmax=3)
# ax1.set_xlim(xlims)
# ax1.set_ylim(ylims)
# ax1.set_zlim(zlims)
# ax1.set_xlabel('x (km)', fontsize=14)
# ax1.set_ylabel('y (km)', fontsize=14)
# # ax1.set_zlabel('z (km)', fontsize=14)
# ax1.set_title(f"Mid-level environmental inflow, height AGL", fontsize=14)
# ax1.view_init(elev=15, azim=-120) # azim viewing angles: -5 for E, -30 for SE, -120 for SW
# plt.tight_layout()
# plt.show()


# fig = plt.figure(figsize=(9,6))
# ax1 = fig.add_subplot(1,1,1, projection='3d')
# ax1.contour(xx[iyy,ixx], yy[iyy,ixx], dbz[iyy,ixx], levels=[30], zdir='z', offset=0, colors='k', linewidths=1)
# ax1.scatter(x_mv1[:,c1_mv1]/1000, y_mv1[:,c1_mv1]/1000, z_mv1[:,c1_mv1]/1000, s=1, c=w_mv1[:,c1_mv1], marker='.', cmap='RdBu_r', vmin=-15, vmax=15)
# ax1.set_xlim(xlims)
# ax1.set_ylim(ylims)
# ax1.set_zlim(zlims)
# ax1.set_xlabel('x (km)', fontsize=14)
# ax1.set_ylabel('y (km)', fontsize=14)
# ax1.set_zlabel('z (km)', fontsize=14)
# ax1.set_title(f"Mid-level environmental inflow, w", fontsize=14)
# ax1.view_init(elev=15, azim=-120) # azim viewing angles: -5 for E, -30 for SE, -120 for SW
# plt.tight_layout()
# plt.show()



#%% Calculate storm motion and save to pickle


ip = '/Users/morgan.schneider/Documents/merger/supercell-125m/'


# Calculate storm motion for streamwise/crosswise vorticity
if 'u_storm' not in locals():
    dbfile = open(ip+'boxes_s1.pkl', 'rb') # interchange with boxes_q
    box = pickle.load(dbfile)
    x1 = 1000 * np.array(box['x1_pp'])
    y1 = 1000 * np.array(box['y1_pp'])
    dbfile.close()
    
    for i in range(len(x1)-1):
        if x1[i+1] == x1[i]:
            if (i != 0):
                x1[i] = (x1[i+1] + x1[i-1]) / 2
            else:
                x1[i+1] = (x1[i+2] + x1[i]) / 2
        
        if y1[i+1] == y1[i]:
            if (i != 0):
                y1[i] = (y1[i+1] + y1[i-1]) / 2
            else:
                y1[i+1] = (y1[i+2] + y1[i]) / 2
    
    u_s = np.zeros(shape=(61,), dtype=float)
    v_s = np.zeros(shape=(61,), dtype=float)
    u_s[1:] = np.gradient(x1, np.linspace(10860, 14400, 60))
    v_s[1:] = np.gradient(y1, np.linspace(10860, 14400, 60))
    
    if u_s[1] == u_s[2]:
        u_s[0] = u_s[1]
    else:
        u_s[0] = u_s[1] - np.diff(u_s)[1]
    
    if v_s[1] == v_s[2]:
        v_s[0] = v_s[1]
    else:
        v_s[0] = v_s[1] - np.diff(v_s)[1]
    
    u_storm = u_s;  v_storm = v_s
    u_storm[1:-1] = movmean(u_s,3)[1:-1]
    v_storm[1:-1] = movmean(v_s,3)[1:-1]


if False:
    data = {'u_storm':u_storm, 'v_storm':v_storm}
    save_to_pickle(data, ip+'storm_motion.pkl')
    
    # f = open('/Users/morgan.schneider/Documents/cm1_pp_v1.1/scripts/storm_motion.input', 'w')
    # f.close()
    # f = open('/Users/morgan.schneider/Documents/cm1_pp_v1.1/scripts/storm_motion.input', 'a')
    # for i in range(len(u_storm)):
    #     f.write(f"{u_storm[i]:.3f}\t{v_storm[i]:.3f}\n")
    # f.close()







#%% Interpolate vorticity to mid-level parcels

from scipy.interpolate import RegularGridInterpolator

mvtime = 220

dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/traj_clusters/traj_clusters_{mvtime}min_v2.pkl", 'rb')
cc = pickle.load(dbfile)
cc_mv1 = cc['mv1']
dbfile.close()

dbfile = open(f"/Users/morgan.schneider/Documents/merger/traj_MV1.pkl", 'rb')
traj_mv1 = pickle.load(dbfile)
pids_ml = traj_mv1[f"{mvtime}min"]['pids'][(cc_mv1 == 1)]
x_ml = traj_mv1[f"{mvtime}min"]['x'][:,(cc_mv1 == 1)]/1000
y_ml = traj_mv1[f"{mvtime}min"]['y'][:,(cc_mv1 == 1)]/1000
z_ml = traj_mv1[f"{mvtime}min"]['z'][:,(cc_mv1 == 1)]/1000
u_ml = traj_mv1[f"{mvtime}min"]['u'][:,(cc_mv1 == 1)]
v_ml = traj_mv1[f"{mvtime}min"]['v'][:,(cc_mv1 == 1)]
zvort_ml = traj_mv1[f"{mvtime}min"]['zvort'][:,(cc_mv1 == 1)]
dbfile.close()

ds = nc.Dataset('/Volumes/Promise_Pegasus_70TB/merger/merger-125m/cm1out_pdata.nc')
ptime = ds.variables['time'][:].data
ds.close()

ds = nc.Dataset(fp+"cm1out_000014.nc")
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
zh = ds.variables['z'][:].data
ds.close()

x_median = np.median(x_ml, axis=1)
y_median = np.median(y_ml, axis=1)
z_median = np.median(z_ml, axis=1)

xlims = [-30,5]
ylims = [-130,-25]
zlims = [0,4.1]

iz = slice(0, np.where(zh >= zlims[1])[0][1])
ix = slice(np.where(xh >= xlims[0])[0][0], np.where(xh >= xlims[1])[0][1])
iy = slice(np.where(yh >= ylims[0])[0][0], np.where(yh >= ylims[1])[0][1])


dbfile = open('/Users/morgan.schneider/Documents/merger/merger-125m/storm_motion.pkl', 'rb')
sm = pickle.load(dbfile)
u_storm = sm['u_storm']
v_storm = sm['v_storm']
dbfile.close()



# Interpolate model grid x/y/h/3d vorticity to trajectories
if False:
    stimes = np.zeros(shape=(61,), dtype=float)
    xvort_ml = np.zeros(shape=(61,len(pids_ml)), dtype=float)
    yvort_ml = np.zeros(shape=(61,len(pids_ml)), dtype=float)
    hvort_ml = np.zeros(shape=(61,len(pids_ml)), dtype=float)
    vort_ml = np.zeros(shape=(61,len(pids_ml)), dtype=float)
    # xvort_median = np.zeros(shape=(61,), dtype=float)
    # yvort_median = np.zeros(shape=(61,), dtype=float)
    # hvort_median = np.zeros(shape=(61,), dtype=float)
    # vort_median = np.zeros(shape=(61,), dtype=float)
    
    for k in np.arange(13,74):
        print(f"cm1out_{k:06d}")
        fp = '/Volumes/Promise_Pegasus_70TB/merger/merger-125m/'
        
        ds = nc.Dataset(fp+f"cm1out_{k:06d}.nc")
        stime = ds.variables['time'][:].data[0]
        xvort = ds.variables['xvort'][:].data[0,iz,iy,ix]
        yvort = ds.variables['yvort'][:].data[0,iz,iy,ix]
        ds.close()
        
        xvort_interp = RegularGridInterpolator((zh[iz], yh[iy], xh[ix]), xvort)
        yvort_interp = RegularGridInterpolator((zh[iz], yh[iy], xh[ix]), yvort)
        
        it = np.where(ptime == stime)[0][0]
        
        for p in range(len(pids_ml)):
            xvort_ml[k-13,p] = xvort_interp((z_ml[it,p], y_ml[it,p], x_ml[it,p]))
            yvort_ml[k-13,p] = yvort_interp((z_ml[it,p], y_ml[it,p], x_ml[it,p]))
            hvort_ml[k-13,p] = np.sqrt(xvort_ml[k-13,p]**2 + yvort_ml[k-13,p]**2)
            vort_ml[k-13,p] = np.sqrt(hvort_ml[k-13,p]**2 + zvort_ml[it,p]**2)
        # xvort_median[k-13] = xvort_interp((z_median[it], y_median[it], x_median[it]))
        # yvort_median[k-13] = yvort_interp((z_median[it], y_median[it], x_median[it]))
        # hvort_median[k-13] = np.median(hvort_ml[k-13,:])
        # vort_median[k-13] = np.median(vort_ml[k-13,:])
        stimes[k-13] = stime/60
    
    ptimes = ptime[0:it+1]/60
    
    xvort_ml_interp = RegularGridInterpolator((stimes, pids_ml), xvort_ml)
    yvort_ml_interp = RegularGridInterpolator((stimes, pids_ml), yvort_ml)
    hvort_ml_interp = RegularGridInterpolator((stimes, pids_ml), hvort_ml)
    vort_ml_interp = RegularGridInterpolator((stimes, pids_ml), vort_ml)
    # xvort_median_interp = RegularGridInterpolator((stimes,), xvort_median)
    # yvort_median_interp = RegularGridInterpolator((stimes,), yvort_median)
    # hvort_median_interp = RegularGridInterpolator((stimes,), hvort_median)
    # vort_median_interp = RegularGridInterpolator((stimes,), vort_median)
    
    pp,tp = np.meshgrid(pids_ml, ptimes, indexing='xy')
    
    xvort_ml = xvort_ml_interp((tp, pp))
    yvort_ml = yvort_ml_interp((tp, pp))
    hvort_ml = hvort_ml_interp((tp, pp))
    vort_ml = vort_ml_interp((tp, pp))
    # xvort_median = xvort_median_interp((ptimes,))
    # yvort_median = yvort_median_interp((ptimes,))
    # hvort_median = hvort_median_interp((ptimes,))
    # vort_median = vort_median_interp((ptimes,))




# Get streamwise/crosswise vorticity and interpolate to trajectories
if False:
    stimes = np.zeros(shape=(61,), dtype=float)
    vort_sw_ml = np.zeros(shape=(61,len(pids_ml)), dtype=float)
    vort_cw_ml = np.zeros(shape=(61,len(pids_ml)), dtype=float)
    # vort_sw_median = np.zeros(shape=(61,), dtype=float)
    # vort_cw_median = np.zeros(shape=(61,), dtype=float)
    
    for k in np.arange(13,74):
        print(f"cm1out_{k:06d}")
        fp = '/Volumes/Promise_Pegasus_70TB/merger/merger-125m/'
        
        ds = nc.Dataset(fp+f"cm1out_{k:06d}.nc")
        stime = ds.variables['time'][:].data[0]
        xvort = ds.variables['xvort'][:].data[0,iz,iy,ix]
        yvort = ds.variables['yvort'][:].data[0,iz,iy,ix]
        ds.close()
        
        u_sr = u_ml - u_storm[k-13]
        v_sr = v_ml - v_storm[k-13]
        
        xvort_interp = RegularGridInterpolator((zh[iz], yh[iy], xh[ix]), xvort)
        yvort_interp = RegularGridInterpolator((zh[iz], yh[iy], xh[ix]), yvort)
        del xvort,yvort
        
        it = np.where(ptime == stime)[0][0]
        
        for p in range(len(pids_ml)):
            usr_ml = u_sr[it,p]
            vsr_ml = v_sr[it,p]
            xvort_ml = xvort_interp((z_ml[it,p], y_ml[it,p], x_ml[it,p]))
            yvort_ml = yvort_interp((z_ml[it,p], y_ml[it,p], x_ml[it,p]))
            vort_sw_ml[k-13,p] = (usr_ml*xvort_ml + vsr_ml*yvort_ml) / np.sqrt(usr_ml**2 + vsr_ml**2)
            vort_cw_ml[k-13,p] = (-vsr_ml*xvort_ml + usr_ml*yvort_ml) / np.sqrt(usr_ml**2 + vsr_ml**2)
        # vort_sw_median[k-13] = np.median(vort_sw_ml[k-13,:])
        # vort_cw_median[k-13] = np.median(vort_cw_ml[k-13,:])
        stimes[k-13] = stime/60
    
    ptimes = ptime[0:it+1]/60
    
    vort_sw_ml_interp = RegularGridInterpolator((stimes, pids_ml), vort_sw_ml)
    vort_cw_ml_interp = RegularGridInterpolator((stimes, pids_ml), vort_cw_ml)
    # vort_sw_median_interp = RegularGridInterpolator((stimes,), vort_sw_median)
    # vort_cw_median_interp = RegularGridInterpolator((stimes,), vort_cw_median)
    
    pp,tp = np.meshgrid(pids_ml, ptimes, indexing='xy')
    
    vort_sw_ml = vort_sw_ml_interp((tp, pp))
    vort_cw_ml = vort_cw_ml_interp((tp, pp))
    # vort_sw_median = vort_sw_median_interp((ptimes,))
    # vort_cw_median = vort_cw_median_interp((ptimes,))



# Interpolate pgf from model grid to trajectories
if False:
    # dbzs = dict()
    
    if 'prs0' not in locals():
        ds = nc.Dataset(fp+"base/cm1out_000001.nc")
        prs0 = ds.variables['prs0'][:].data[0,iz,iy,ix]
        ds.close()
    
    stimes = np.zeros(shape=(61,), dtype=float)
    xpgf_ml = np.zeros(shape=(61,len(pids_ml)), dtype=float)
    ypgf_ml = np.zeros(shape=(61,len(pids_ml)), dtype=float)
    hpgf_ml = np.zeros(shape=(61,len(pids_ml)), dtype=float)
    pgf_ml = np.zeros(shape=(61,len(pids_ml)), dtype=float)
    xpgf_median = np.zeros(shape=(61,), dtype=float)
    ypgf_median = np.zeros(shape=(61,), dtype=float)
    hpgf_median = np.zeros(shape=(61,), dtype=float)
    pgf_median = np.zeros(shape=(61,), dtype=float)
    
    for k in np.arange(13,74):
        print(f"cm1out_{k:06d}")
        fp = '/Volumes/Promise_Pegasus_70TB/merger/merger-125m/'
        
        ds = nc.Dataset(fp+f"cm1out_{k:06d}.nc")
        stime = ds.variables['time'][:].data[0]
        # if 'dbz2' in ds.variables:
        #     dbz = ds.variables['dbz2'][:].data[0,0,iy,ix]
        # else:
        #     dbz = ds.variables['dbz'][:].data[0,0,iy,ix]
        prspert = ds.variables['prs'][:].data[0,iz,iy,ix] - prs0
        xpgf = -1/1.1 * np.gradient(prspert, xh[ix]*1000, axis=2)
        ypgf = -1/1.1 * np.gradient(prspert, yh[iy]*1000, axis=1)
        zpgf = -1/1.1 * np.gradient(prspert, zh[iz]*1000, axis=0)
        pgf = (xpgf**2 + ypgf**2 + zpgf**2)**0.5
        hpgf = (xpgf**2 + ypgf**2)**0.5
        del zpgf
        ds.close()
        
        xpgf_interp = RegularGridInterpolator((zh[iz], yh[iy], xh[ix]), xpgf)
        ypgf_interp = RegularGridInterpolator((zh[iz], yh[iy], xh[ix]), ypgf)
        hpgf_interp = RegularGridInterpolator((zh[iz], yh[iy], xh[ix]), hpgf)
        pgf_interp = RegularGridInterpolator((zh[iz], yh[iy], xh[ix]), pgf)
        
        it = np.where(ptime == stime)[0][0]
        
        for p in range(len(pids_ml)):
            xpgf_ml[k-13,p] = xpgf_interp((z_ml[it,p], y_ml[it,p], x_ml[it,p]))
            ypgf_ml[k-13,p] = ypgf_interp((z_ml[it,p], y_ml[it,p], x_ml[it,p]))
            hpgf_ml[k-13,p] = hpgf_interp((z_ml[it,p], y_ml[it,p], x_ml[it,p]))
            pgf_ml[k-13,p] = pgf_interp((z_ml[it,p], y_ml[it,p], x_ml[it,p]))
        xpgf_median[k-13] = xpgf_interp((z_median[it], y_median[it], x_median[it]))
        ypgf_median[k-13] = ypgf_interp((z_median[it], y_median[it], x_median[it]))
        hpgf_median[k-13] = hpgf_interp((z_median[it], y_median[it], x_median[it]))
        pgf_median[k-13] = pgf_interp((z_median[it], y_median[it], x_median[it]))
        stimes[k-13] = stime/60
    
    ptimes = ptime[0:it+1]/60
    
    xpgf_ml_interp = RegularGridInterpolator((stimes, pids_ml), xpgf_ml)
    ypgf_ml_interp = RegularGridInterpolator((stimes, pids_ml), ypgf_ml)
    hpgf_ml_interp = RegularGridInterpolator((stimes, pids_ml), hpgf_ml)
    pgf_ml_interp = RegularGridInterpolator((stimes, pids_ml), pgf_ml)
    xpgf_median_interp = RegularGridInterpolator((stimes,), xpgf_median)
    ypgf_median_interp = RegularGridInterpolator((stimes,), ypgf_median)
    hpgf_median_interp = RegularGridInterpolator((stimes,), hpgf_median)
    pgf_median_interp = RegularGridInterpolator((stimes,), pgf_median)
    
    pp,tp = np.meshgrid(pids_ml, ptimes, indexing='xy')
    
    xpgf_ml = xpgf_ml_interp((tp, pp))
    ypgf_ml = ypgf_ml_interp((tp, pp))
    hpgf_ml = hpgf_ml_interp((tp, pp))
    pgf_ml = pgf_ml_interp((tp, pp))
    xpgf_median = xpgf_median_interp((ptimes,))
    ypgf_median = ypgf_median_interp((ptimes,))
    hpgf_median = hpgf_median_interp((ptimes,))
    pgf_median = pgf_median_interp((ptimes,))

    
    
    
#%% Save interpolations to pickle

if False:
    # dbfile = open(ip+'hppga_traj_210min.pkl', 'wb')
    # save_vars = {'xpgf_ml':xpgf_ml, 'ypgf_ml':ypgf_ml, 'hpgf_ml':hpgf_ml, 'pgf_ml':pgf_ml,
    #               'xpgf_median':xpgf_median, 'ypgf_median':ypgf_median,
    #               'hpgf_median':hpgf_median, 'pgf_median':pgf_median}
    # pickle.dump(save_vars, dbfile)
    # dbfile.close()
    
    dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/hvort_traj_{mvtime}min.pkl", 'wb')
    save_vars = {'xvort_ml':xvort_ml, 'yvort_ml':yvort_ml,
                 'hvort_ml':hvort_ml, 'vort_ml':vort_ml}
    pickle.dump(save_vars, dbfile)
    dbfile.close()

if False:
    # dbfile = open(ip+f"hppga_traj_{mvtime}min.pkl", 'rb')
    # tmp = pickle.load(dbfile)
    # dbfile.close()
    
    # save_vars = {'xpgf_ml':xpgf_ml, 'ypgf_ml':ypgf_ml, 'hpgf_ml':hpgf_ml, 'pgf_ml':pgf_ml,
    #               'xpgf_median':xpgf_median, 'ypgf_median':ypgf_median,
    #               'hpgf_median':hpgf_median, 'pgf_median':pgf_median}
    # tmp.update(save_vars)
    # dbfile = open(ip+f"hppga_traj_{mvtime}min.pkl", 'wb')
    # pickle.dump(tmp, dbfile)
    # dbfile.close()
    
    dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/hvort_traj_{mvtime}min.pkl", 'rb')
    tmp = pickle.load(dbfile)
    dbfile.close()
    
    save_vars = {'vort_sw_ml':vort_sw_ml, 'vort_cw_ml':vort_cw_ml}
    tmp.update(save_vars)
    dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/hvort_traj_{mvtime}min.pkl", 'wb')
    pickle.dump(tmp, dbfile)
    dbfile.close()


# correcting vort_cw maybe?
if False:
    for i in [210,220,225]:
        dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/hvort_traj_{i}min.pkl", 'rb')
        tmp = pickle.load(dbfile)
        dbfile.close()
        
        vort_cw_wrong = tmp['vort_cw_ml_signed']
        vort_cw_right = -1 * vort_cw_wrong
        
        dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/hvort_traj_{i}min.pkl", 'wb')
        tmp.update({'vort_cw_ml':vort_cw_right})
        pickle.dump(tmp, dbfile)
        dbfile.close()




#%% Plot interpolated vorticity time series for mid-level parcels (don't use)

from matplotlib.ticker import MultipleLocator

dbfile = open(f"/Users/morgan.schneider/Documents/merger/traj_MV1.pkl", 'rb')
traj = pickle.load(dbfile)
dbfile.close()

times = [210, 220]

ds = nc.Dataset('/Volumes/Promise_Pegasus_70TB/merger/merger-125m/cm1out_pdata.nc')
ptime = ds.variables['time'][:].data/60
ds.close()

# n = np.where(ptime == 220)[0][0]


fig,axs = plt.subplots(4, 1, figsize=(9,8), sharex=True, layout='constrained')

for a in range(len(axs)):
    axs[a].axhline(0, color='darkgray', linestyle='-', linewidth=2.25)


for i in range(len(times)):
    t = times[i]
    it1 = np.where(ptime >= t-10)[0][0]
    it2 = np.where(ptime > t)[0][0]
    it = slice(it1,it2)
    
    dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/traj_clusters/traj_clusters_{t:.0f}min_v2.pkl", 'rb')
    c = pickle.load(dbfile)
    cc = c['mv1']
    dbfile.close()
    
    w_ml = traj[f"{t:.0f}min"]['w'][:,(cc == 1)]
    zvort_ml = traj[f"{t:.0f}min"]['zvort'][:,(cc == 1)]
    
    dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/hvort_traj_{t:.0f}min.pkl", 'rb')
    vort_traj = pickle.load(dbfile)
    hvort_ml = vort_traj['hvort_ml']
    vort_ml = vort_traj['vort_ml']
    vort_sw = vort_traj['vort_sw_ml']
    vort_cw = vort_traj['vort_cw_ml']
    dbfile.close()
    
    
    
    # colorblind-friendly colors
    # col = ['red', 'gold', 'deepskyblue', 'mediumblue']
    col = ['black', 'midnightblue', 'mediumblue', 'red', 'purple', 'deepskyblue']
    al = [0.2, 0.2, 0.3, 0.3, 0.3, 0.3]
    
    
    axs[0].fill_between(ptime[it], np.percentile(w_ml[it,:], 25, axis=1),
                        np.percentile(w_ml[it,:], 75, axis=1), alpha=al[0], color=col[0])
    
    axs[1].fill_between(ptime[it], np.percentile(vort_ml[it,:], 25, axis=1),
                        np.percentile(vort_ml[it,:], 75, axis=1), alpha=al[1], color=col[1])
    
    axs[2].fill_between(ptime[it], np.percentile(hvort_ml[it,:], 25, axis=1),
                        np.percentile(hvort_ml[it,:], 75, axis=1), alpha=al[2], color=col[2])
    
    axs[2].fill_between(ptime[it], np.percentile(zvort_ml[it,:], 25, axis=1),
                        np.percentile(zvort_ml[it,:], 75, axis=1), alpha=al[3], color=col[3])
    
    axs[3].fill_between(ptime[it], np.percentile(vort_sw[it,:], 25, axis=1),
                        np.percentile(vort_sw[it,:], 75, axis=1), alpha=al[4], color=col[4])
    
    axs[3].fill_between(ptime[it], np.percentile(vort_cw[it,:], 25, axis=1),
                        np.percentile(vort_cw[it,:], 75, axis=1), alpha=al[5], color=col[5])
        
    s1, = axs[0].plot(ptime[it], np.median(w_ml[it,:], axis=1), col[0], linewidth=2)
    s2, = axs[1].plot(ptime[it], np.median(vort_ml[it,:], axis=1), col[1], linewidth=2)
    s3, = axs[2].plot(ptime[it], np.median(hvort_ml[it,:], axis=1), col[2], linewidth=2)
    s4, = axs[2].plot(ptime[it], np.median(zvort_ml[it,:], axis=1), col[3], linewidth=2)
    s5, = axs[3].plot(ptime[it], np.median(vort_sw[it,:], axis=1), col[4], linewidth=2)
    s6, = axs[3].plot(ptime[it], np.median(vort_cw[it,:], axis=1), col[5], linewidth=2)
    

axs[3].set_xlabel('Time (min)', fontsize=16)
axs[3].tick_params(axis='x', which='major', labelsize=14)
axs[3].set_xlim([200,220])

axs[0].legend(handles=[s1], labels=["w"], loc=3, fontsize=12)
axs[1].legend(handles=[s2], labels=["Total |\u03c9|"], loc=3, fontsize=12)
axs[2].legend(handles=[s3,s4], labels=["|\u03c9$_H$|", "\u03B6"], loc=3, fontsize=12, ncol=2)
axs[3].legend(handles=[s5,s6], labels=["Streamwise", "Crosswise"], loc=3, fontsize=12, ncol=2)
axs[0].set_ylim([-15, 15]) # w
axs[1].set_ylim([-0.02, 0.08]) # total -- [-0.02, 0.1] for 225
axs[2].set_ylim([-0.04, 0.06]) # horiz and vert -- [-0.05, 0.09] for 225
axs[3].set_ylim([-0.06, 0.06]) # streamwise and crosswise -- [-0.07, 0.07] for 225
# axs[0].set_yticks([-15, -10, -5, 0, 5, 10, 15])
# # axs[1].set_yticks([-0.04, -0.02, 0, 0.02, 0.04, 0.06, 0.08])
# # axs[1].set_yticks([-0.05, -0.025, 0, 0.025, 0.05, 0.075, 0.1])
# axs[1].set_yticks([-0.06, -0.03, 0, 0.03, 0.06, 0.09])
# axs[2].set_yticks([0, 0.025, 0.05, 0.075, 0.1])
# axs[3].set_yticks([-0.06, -0.04, -0.02, 0, 0.02, 0.04, 0.06])
    

# axs[0].set_title("Mid-level parcel vorticity for merger MV at 210/225 min", fontsize=14)
axs[0].text(203.5, 16, '210 min', fontsize=18)
axs[0].text(213.5, 16, '220 min', fontsize=18)

for i in range(len(axs)):
    axs[i].axvline(210, color='k', linestyle='--', linewidth=2)
    axs[i].tick_params(axis='y', which='major', labelsize=12)
    axs[i].grid(visible=True, which='major', color='darkgray', linestyle='-')
    axs[i].grid(visible=True, which='minor', color='lightgray', linestyle='--')

axs[0].yaxis.set_major_locator(MultipleLocator(5))
axs[1].yaxis.set_major_locator(MultipleLocator(0.02))
axs[2].yaxis.set_major_locator(MultipleLocator(0.02))
axs[3].yaxis.set_major_locator(MultipleLocator(0.02))
# axs[3].yaxis.set_minor_locator(MultipleLocator(0.01))
axs[3].xaxis.set_major_locator(MultipleLocator(5))
axs[3].xaxis.set_minor_locator(MultipleLocator(1))

# plt.show()

figsave = False

if figsave:
    plt.savefig(f"/Users/morgan.schneider/Documents/merger/traj_vort_timeseries_210+220_old.png", dpi=300)


#%% Vorticity time series v2 ***PAPER FIG***

from matplotlib.ticker import MultipleLocator

dbfile = open(f"/Users/morgan.schneider/Documents/merger/traj_MV1.pkl", 'rb')
traj = pickle.load(dbfile)
dbfile.close()

times = [210, 220]

ds = nc.Dataset('/Volumes/Promise_Pegasus_70TB/merger/merger-125m/cm1out_pdata.nc')
ptime = ds.variables['time'][:].data/60
ds.close()

# n = np.where(ptime == 220)[0][0]


fig,axs = plt.subplots(4, 2, figsize=(9,8))

for i in range(2):
    for j in range(4):
        axs[j,i].axhline(0, color='darkgray', linestyle='-', linewidth=2.25)

for i in range(len(times)):
    t = times[i]
    it1 = np.where(ptime >= t-15)[0][0]
    it2 = np.where(ptime > t)[0][0]
    it = slice(it1,it2)
    
    dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/traj_clusters/traj_clusters_{t:.0f}min_v2.pkl", 'rb')
    c = pickle.load(dbfile)
    cc = c['mv1']
    dbfile.close()
    
    w_ml = traj[f"{t:.0f}min"]['w'][:,(cc == 1)]
    zvort_ml = traj[f"{t:.0f}min"]['zvort'][:,(cc == 1)]
    
    dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/hvort_traj_{t:.0f}min.pkl", 'rb')
    vort_traj = pickle.load(dbfile)
    hvort_ml = vort_traj['hvort_ml']
    vort_ml = vort_traj['vort_ml']
    vort_sw = vort_traj['vort_sw_ml']
    vort_cw = vort_traj['vort_cw_ml']
    dbfile.close()
    
    
    # colorblind-friendly colors
    col = ['black', 'midnightblue', 'mediumblue', 'red', 'deepskyblue', 'purple']
    al = [0.2, 0.2, 0.3, 0.3, 0.3, 0.3]
    
    
    axs[0,i].fill_between(ptime[it], np.percentile(w_ml[it,:], 25, axis=1),
                        np.percentile(w_ml[it,:], 75, axis=1), alpha=al[0], color=col[0])
    
    axs[1,i].fill_between(ptime[it], np.percentile(vort_ml[it,:], 25, axis=1),
                        np.percentile(vort_ml[it,:], 75, axis=1), alpha=al[1], color=col[1])
    
    axs[2,i].fill_between(ptime[it], np.percentile(hvort_ml[it,:], 25, axis=1),
                        np.percentile(hvort_ml[it,:], 75, axis=1), alpha=al[2], color=col[2])
    
    axs[2,i].fill_between(ptime[it], np.percentile(zvort_ml[it,:], 25, axis=1),
                        np.percentile(zvort_ml[it,:], 75, axis=1), alpha=al[3], color=col[3])
    
    axs[3,i].fill_between(ptime[it], np.percentile(vort_sw[it,:], 25, axis=1),
                        np.percentile(vort_sw[it,:], 75, axis=1), alpha=al[4], color=col[4])
    
    axs[3,i].fill_between(ptime[it], np.percentile(vort_cw[it,:], 25, axis=1),
                        np.percentile(vort_cw[it,:], 75, axis=1), alpha=al[5], color=col[5])
        
    s1, = axs[0,i].plot(ptime[it], np.median(w_ml[it,:], axis=1), col[0], linewidth=2)
    s2, = axs[1,i].plot(ptime[it], np.median(vort_ml[it,:], axis=1), col[1], linewidth=2)
    s3, = axs[2,i].plot(ptime[it], np.median(hvort_ml[it,:], axis=1), col[2], linewidth=2)
    s4, = axs[2,i].plot(ptime[it], np.median(zvort_ml[it,:], axis=1), col[3], linewidth=2)
    s5, = axs[3,i].plot(ptime[it], np.median(vort_sw[it,:], axis=1), col[4], linewidth=2)
    s6, = axs[3,i].plot(ptime[it], np.median(vort_cw[it,:], axis=1), col[5], linewidth=2)


axs[0,0].legend(handles=[s1], labels=["w"], loc=3, fontsize=11)
axs[1,0].legend(handles=[s2], labels=["Total |\u03c9|"], loc=3, fontsize=11)
axs[2,0].legend(handles=[s3,s4], labels=["|\u03c9$_H$|", "\u03B6"], loc=3, fontsize=11, ncol=2)
axs[3,0].legend(handles=[s5,s6], labels=["\u03c9$_{SW}$", "\u03c9$_{CW}$"], loc=3, fontsize=11, ncol=2)
axs[3,0].set_xlabel('                                                        Time (min)', fontsize=14)

axs[0,0].set_ylabel("Vertical\n velocity ", fontsize=12)
axs[1,0].set_ylabel("Total\n vorticity ", fontsize=12)
axs[2,0].set_ylabel("Horiz. and\n vert. vorticity ", fontsize=12)
axs[3,0].set_ylabel("Horiz. vorticity\n components ", fontsize=12)

# axs[0,0].set_ylabel("Vertical\n velocity ", fontsize=12)
# axs[1,0].set_ylabel("Vorticity\n Total |\u03c9| ", fontsize=12)
# axs[2,0].set_ylabel("Vorticity\n |\u03c9$_H$| and \u03B6 ", fontsize=12)
# axs[3,0].set_ylabel("Vorticity\n \u03c9$_H$ components ", fontsize=12)


for i in range(2):
    axs[0,i].set_ylim([-15, 15]) # w
    axs[1,i].set_ylim([-0.02, 0.08]) # total -- [-0.02, 0.1] for 225
    axs[2,i].set_ylim([-0.04, 0.06]) # horiz and vert -- [-0.05, 0.09] for 225
    axs[3,i].set_ylim([-0.06, 0.06]) # streamwise and crosswise -- [-0.07, 0.07] for 225
    axs[3,i].tick_params(axis='x', which='major', labelsize=12)
    
    axs[0,i].yaxis.set_major_locator(MultipleLocator(5))
    axs[1,i].yaxis.set_major_locator(MultipleLocator(0.02))
    axs[2,i].yaxis.set_major_locator(MultipleLocator(0.02))
    axs[3,i].yaxis.set_major_locator(MultipleLocator(0.02))
    # axs[3,i].yaxis.set_minor_locator(MultipleLocator(0.01))
    
    for j in range(4):
        axs[j,i].xaxis.set_major_locator(MultipleLocator(5))
        axs[j,i].xaxis.set_minor_locator(MultipleLocator(1))
        axs[j,i].tick_params(axis='y', which='major', labelsize=10)
        axs[j,i].grid(visible=True, which='major', color='darkgray', linestyle='-')
        axs[j,i].grid(visible=True, which='minor', color='lightgray', linestyle='--')
        if i > 0:
            axs[j,i].set_yticklabels([])

for j in range(4):
    axs[j,0].set_xlim([195,210])
    axs[j,1].set_xlim([205,220])
    if j < 3:
        for i in range(2):
            axs[j,i].set_xticklabels([])
    
# axs[0].set_yticks([-15, -10, -5, 0, 5, 10, 15])
# # axs[1].set_yticks([-0.04, -0.02, 0, 0.02, 0.04, 0.06, 0.08])
# # axs[1].set_yticks([-0.05, -0.025, 0, 0.025, 0.05, 0.075, 0.1])
# axs[1].set_yticks([-0.06, -0.03, 0, 0.03, 0.06, 0.09])
# axs[2].set_yticks([0, 0.025, 0.05, 0.075, 0.1])
# axs[3].set_yticks([-0.06, -0.04, -0.02, 0, 0.02, 0.04, 0.06])

axs[3,0].set_xticks([195, 200, 205, 210])
axs[3,0].set_xticklabels(["    195", "200", "205", "210   "])

axs[3,1].set_xticks([205, 210, 215, 220])
axs[3,1].set_xticklabels(["    205", "210", "215", "220   "])

# axs[0].set_title("Mid-level parcel vorticity for merger MV at 210/225 min", fontsize=14)
axs[0,0].text(200.5, 16, '210 min', fontsize=16)
axs[0,1].text(210.5, 16, '220 min', fontsize=16)

fig.align_ylabels(axs[:,0])

plt.subplots_adjust(wspace=0.06, hspace=0.15)

# plt.show()

figsave = False

if figsave:
    plt.savefig(f"/Users/morgan.schneider/Documents/merger/traj_vort_timeseries.png", dpi=300)


#%% 3D plot animations colored by whatever variables - mid level MERGER source

# data = ypgf_ml; datalims = [-0.02,0.02]; cm = 'seismic'; title_str = 'y-direction PPGA'; file_str = 'YPPGA'
# data = ppga_ml; datalims = [0,0.1]; cm = 'HomeyerRainbow'; title_str = 'Total PPGA'; file_str = 'PPGA'
# data = hppga_ml; datalims = [0,0.1]; cm = 'HomeyerRainbow'; title_str = 'Horizontal PPGA'; file_str = 'HPPGA'
# data = z_ml; datalims = [0,3.5]; cm = 'HomeyerRainbow'; title_str = 'Height AGL'; file_str = 'z'
# data = traj_mv1[f"210min"]['w'][:,(cc_mv1 == 1)]; datalims = [-10,10]; cm = 'coolwarm'; title_str = 'Vertical velocity'; file_str = 'w'
# data = traj_mv1[f"210min"]['b'][:,(cc_mv1 == 1)]; datalims = [-0.2,0.05]; cm = cmocean.cm.curl; title_str = 'Buoyancy'; file_str = 'buoy'
data = traj_mv1[f"210min"]['zvort'][:,(cc_mv1 == 1)]; datalims = [-0.03,0.03]; cm = 'coolwarm'; title_str = 'Vertical vorticity'; file_str = 'zvort'



from matplotlib.animation import FuncAnimation
from matplotlib.colors import TwoSlopeNorm
if datalims[0] == 0:
    norm = TwoSlopeNorm(vcenter=0.5*datalims[1], vmin=datalims[0], vmax=datalims[1])
    datalims = [None,None]
else:
    norm = TwoSlopeNorm(vcenter=0, vmin=datalims[0], vmax=datalims[1])
    datalims = [None,None]

view_angle = -150

fig = plt.figure(figsize=(6,6))
axa = fig.add_subplot(1,1,1, projection='3d')
c = axa.contour(xx, yy, dbzs['0'], levels=[30], zdir='z', offset=0, colors='k', linewidths=1)
axa.scatter(x_ml[0,:], y_ml[0,:], z_ml[0,:], c=data[0,:], s=5, marker='.', norm=norm, cmap=cm, vmin=datalims[0], vmax=datalims[1])
axa.set_xlim(xlims)
axa.set_ylim(ylims)
axa.set_zlim(zlims)
axa.set_xlabel('x (km)', fontsize=14)
axa.set_ylabel('y (km)', fontsize=14)
# axa.set_zlabel('z (km)', fontsize=14)
axa.set_title(f"Mid-level inflow source, MV parcels at 210 min \n {title_str}", fontsize=14)
axa.view_init(elev=15, azim=view_angle)
# plt.colorbar(s, ax=axa, extend='both')
plt.tight_layout()


def animate_traj(i):
    global axa
    axa.clear()
    
    c = axa.contour(xx, yy, dbzs[f"{i}"], levels=[30], zdir='z', offset=0, colors='k', linewidths=1)
    axa.scatter(x_ml[0:i+1,:], y_ml[0:i+1,:], z_ml[0:i+1,:], c=data[0:i+1,:], s=5, marker='.', norm=norm, cmap=cm, vmin=datalims[0], vmax=datalims[1])
    axa.set_xlim(xlims)
    axa.set_ylim(ylims)
    axa.set_zlim(zlims)
    axa.set_xlabel('x (km)', fontsize=14)
    axa.set_ylabel('y (km)', fontsize=14)
    # axa.set_zlabel('z (km)', fontsize=14)
    axa.set_title(f"Mid-level inflow source, MV parcels at 210 min \n {title_str}", fontsize=14)
    axa.view_init(elev=15, azim=view_angle)
    # plt.colorbar(s, ax=axa, extend='both')
    plt.tight_layout()

figsave = True

anim = FuncAnimation(fig, animate_traj, frames=len(dbzs), interval=50, repeat=False, blit=False)
if figsave:
    anim.save(ip+f"MV1_traj_{file_str}.gif", dpi=300)
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
    



