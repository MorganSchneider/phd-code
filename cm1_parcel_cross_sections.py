#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 15:53:42 2024

@author: morgan.schneider
"""

####################
### Load modules ###
####################

from CM1utils import *

#%% Load parcel data

mv_time = 210
fnum = 43

fp = '/Volumes/Promise_Pegasus_70TB/merger/merger-125m/'
ip = '/Users/morgan.schneider/Documents/merger/merger-125m/'

# Read parcel data
ds = nc.Dataset(fp+'cm1out_pdata.nc')
ptime = ds.variables['time'][:].data
xp = ds.variables['x'][:].data
yp = ds.variables['y'][:].data
zp = ds.variables['z'][:].data
ds.close()

dbfile = open('/Users/morgan.schneider/Documents/merger/traj_MV1.pkl', 'rb')
traj = pickle.load(dbfile)
pids_mv = traj[f"{mv_time}min"]['pids']
x_mv = xp[:,int(pids_mv)]/1000
y_mv = yp[:,int(pids_mv)]/1000
z_mv = zp[:,int(pids_mv)]/1000
# x_mv = traj[f"{mv_time}min"]['x']/1000
# y_mv = traj[f"{mv_time}min"]['y']/1000
# z_mv = traj[f"{mv_time}min"]['z']/1000
# w_mv = traj[f"{mv_time}min"]['w']
# b_mv = traj[f"{mv_time}min"]['b']
# zvort_mv = traj[f"{mv_time}min"]['zvort']
# vpg_mv = traj[f"{mv_time}min"]['vpg']
dbfile.close()

dbfile = open(ip+f"traj_clusters_{mv_time}min_v2.pkl", 'rb')
c = pickle.load(dbfile)
cc = c['mv1']
dbfile.close()

# Mid-level source
x_ll = x_mv[:,(cc==0)]
y_ll = y_mv[:,(cc==0)]
z_ll = z_mv[:,(cc==0)]

x_ml = x_mv[:,(cc==1)]
y_ml = y_mv[:,(cc==1)]
z_ml = z_mv[:,(cc==1)]

x_sc = x_mv[:,(cc==2)]
y_sc = y_mv[:,(cc==2)]
z_sc = z_mv[:,(cc==2)]

x_qc = x_mv[:,(cc==3)]
y_qc = y_mv[:,(cc==3)]
z_qc = z_mv[:,(cc==3)]

#%% Load model grid data

fn = 43
# 181 min -> 14 | 195 min -> 28 | 210 min -> 43 | 225 min -> 58 | 240 min -> 73

if fn == 13:
    fname = f"base/cm1out_{fn:06d}.nc"
else:
    fname = f"cm1out_{fn:06d}.nc"

# Read output file
ds = nc.Dataset(fp+fname)
stime = ds.variables['time'][:].data[0]
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
z = ds.variables['z'][:].data

iz5 = np.where(z >= 5)[0][0]
iz = slice(0,iz5+1)

ti = np.where(ptime == stime)[0][0]
y_median = np.median(y_ml[ti,:])
idy = (np.abs(yh - y_median)).argmin()
# x_median = np.median(x_ml[ti,:])
# idx = (np.abs(xh - x_median)).argmin()

xlims = [-30,10] #[-55,25]
# ylims = [-60,20] #[-100,-20]
# ylims = [yh[idy]-2, yh[idy]+2]

ix = slice(np.where(xh >= xlims[0])[0][0], np.where(xh >= xlims[1])[0][1])
# iy = slice(np.where(yh >= ylims[0])[0][0], np.where(yh >= ylims[1])[0][1])
iy = slice(idy-4,idy+5)

prs = ds.variables['prs'][:].data[0,iz,idy,ix]
winterp = ds.variables['winterp'][:].data[0,iz,idy,ix]
# uinterp = ds.variables['uinterp'][:].data[0,iz,iy,ix]
# vinterp = ds.variables['vinterp'][:].data[0,iz,iy,ix]
# zvort = ds.variables['zvort'][:].data[0,iz,iy,ix]
# S_N = np.gradient(uinterp, xh[ix]*1000, axis=2) - np.gradient(vinterp, yh[iy]*1000, axis=1)
# S_S = np.gradient(vinterp, xh[ix]*1000, axis=2) + np.gradient(uinterp, yh[iy]*1000, axis=1)
# OW = np.min((S_N**2 + S_S**2 - zvort**2), axis=1)
# del S_N,S_S,zvort,uinterp,vinterp

if fn == 13:
    thr = ds.variables['th'][:].data[0,iz,idy,ix] * (1 + 0.61*ds.variables['qv'][:].data[0,iz,idy,ix] - 
                (ds.variables['qc'][:].data[0,iz,idy,ix] + ds.variables['qr'][:].data[0,iz,idy,ix] + 
                  ds.variables['qi'][:].data[0,iz,idy,ix] + ds.variables['qs'][:].data[0,iz,idy,ix] + 
                  ds.variables['qg'][:].data[0,iz,idy,ix] + ds.variables['qhl'][:].data[0,iz,idy,ix]))
    thr0 = ds.variables['th0'][:].data[0,iz,idy,ix] * (1 + 0.61*ds.variables['qv0'][:].data[0,iz,idy,ix])
    thrpert = thr - thr0
    # B = 9.8 * (thrpert/thr0)
    del thr,thr0
ds.close()

ds = nc.Dataset(fp+f"base/cm1out_000001.nc")
prspert = (prs - ds.variables['prs0'][:].data[0,iz,idy,ix])/100
del prs
ds.close()

if fn > 13:
    ds = nc.Dataset(fp+f"pp/dyn_{fn:06d}.nc")
    thrpert = ds.variables['thrpert'][:].data[iz,idy,ix]
    # thr0 = np.moveaxis(np.tile(ds.variables['thr0'][:].data, (thrpert1.shape[1], thrpert1.shape[2], 1)), -1, 0)
    # B = 9.8 * (thrpert/thr0)
    # del thr0
    ds.close()


#%% Plot cross sections

# could color the parcels by buoyancy or something

w_lims = [-15,15]
thr_lims = [-10,10]
pp_lims = [-4,4]
# OW_lims = [-0.01,0]
zlims = [0,3.5]

qix = int(np.round(len(xh[(xh>=xlims[0]) & (xh<=xlims[1])])/60)*2) # scale wind vector spacing to number of grid pts
qiz = 3


figsave = False

fig,ax = plt.subplots(1,1,figsize=(9,4))
plot_cfill(xh[ix], z[iz], winterp, 'w', ax, datalims=w_lims, xlims=xlims, ylims=zlims)
ax.contour(xh[ix], z[iz], thrpert, levels=[-5,-1], colors='purple', linewidths=[2,1], linestyles='-')
# ax.contour(xh[ix], z[iz], OW, levels=[-0.01,-0.0025,-0.001], colors='hotpink', linewidths=1, linestyles='-')
ax.scatter(x_ml[ti,:], z_ml[ti,:], s=10, color='k', marker='.')
ax.set_xlabel('x (km)')
ax.set_ylabel('z (km)')
ax.set_title(f"Mid-level parcels in the MV at {mv_time} min \n t={stime/60:.0f} min")
if figsave:
    plt.savefig()

thr_cm = cmocean.tools.lighten(cmocean.cm.curl, 0.8)

fig,ax = plt.subplots(1,1,figsize=(9,4))
plot_cfill(xh[ix], z[iz], thrpert, 'thrpert', ax, datalims=thr_lims, cmap=thr_cm, xlims=xlims, ylims=zlims)
ax.contour(xh[ix], z[iz], prspert, levels=[1,2,3,4,5], colors='r', linewidths=1)
ax.contour(xh[ix], z[iz], prspert, levels=[-5,-4,-3,-2,-1], colors='b', linewidths=1, linestyles='-')
# ax.contour(xh[ix], z[iz], OW, levels=[-0.01,-0.0025,-0.001], colors='hotpink', linewidths=1, linestyles='-')
ax.quiver(xh[ix][::qix], z[iz][::qiz], 0*winterp[::qiz,::qix], winterp[::qiz,::qix], scale=250, width=0.002, pivot='tail')
ax.scatter(x_ml[ti,:], z_ml[ti,:], s=10, color='k', marker='.')
ax.set_xlabel('x (km)')
ax.set_ylabel('z (km)')
ax.set_title(f"Mid-level parcels in the MV at {mv_time} min \n t={stime/60:.0f} min")


plt.show()



#%% Load model grid data in loop

from matplotlib.animation import FuncAnimation

new_pkl = False
update_pkl = False

for fn in np.arange(13,fnum+1):
    if fn == 13:
        fname = f"base/cm1out_{fn:06d}.nc"
    else:
        fname = f"cm1out_{fn:06d}.nc"
    
    print(f"cm1out_{fn:06d}")
    # Read output file
    ds = nc.Dataset(fp+fname)
    stime = ds.variables['time'][:].data[0]
    xh = ds.variables['xh'][:].data
    yh = ds.variables['yh'][:].data
    z = ds.variables['z'][:].data
    
    iz5 = np.where(z >= 5)[0][0]
    iz = slice(0,iz5+1)
    
    ti = np.where(ptime == stime)[0][0]
    y_median = np.median(y_ml[ti,:])
    idy = (np.abs(yh - y_median)).argmin()
    # x_median = np.median(x_ml[ti,:])
    # idx = (np.abs(xh - x_median)).argmin()
    
    xlims = [-30,10] #[-55,25]
    
    ix = slice(np.where(xh >= xlims[0])[0][0], np.where(xh >= xlims[1])[0][1])
    iy = slice(idy-4,idy+5)
    
    # prs = ds.variables['prs'][:].data[0,iz,idy,ix]
    # winterp = ds.variables['winterp'][:].data[0,iz,idy,ix]
    
    uinterp = ds.variables['uinterp'][:].data[0,iz,iy,ix]
    vinterp = ds.variables['vinterp'][:].data[0,iz,iy,ix]
    zvort = ds.variables['zvort'][:].data[0,iz,idy,ix]
    S_N = np.gradient(uinterp, xh[ix]*1000, axis=2) - np.gradient(vinterp, yh[iy]*1000, axis=1)
    S_S = np.gradient(vinterp, xh[ix]*1000, axis=2) + np.gradient(uinterp, yh[iy]*1000, axis=1)
    # OW = np.min((S_N**2 + S_S**2 - zvort**2), axis=1)
    OW = S_N[:,4,:]**2 + S_S[:,4,:]**2 - zvort**2
    del S_N,S_S,uinterp,vinterp
    
    # if fn == 13:
    #     thr = ds.variables['th'][:].data[0,iz,idy,ix] * (1 + 0.61*ds.variables['qv'][:].data[0,iz,idy,ix] - 
    #                 (ds.variables['qc'][:].data[0,iz,idy,ix] + ds.variables['qr'][:].data[0,iz,idy,ix] + 
    #                   ds.variables['qi'][:].data[0,iz,idy,ix] + ds.variables['qs'][:].data[0,iz,idy,ix] + 
    #                   ds.variables['qg'][:].data[0,iz,idy,ix] + ds.variables['qhl'][:].data[0,iz,idy,ix]))
    #     thr0 = ds.variables['th0'][:].data[0,iz,idy,ix] * (1 + 0.61*ds.variables['qv0'][:].data[0,iz,idy,ix])
    #     thrpert = thr - thr0
    #     B = 9.8 * (thrpert/thr0)
    #     del thr,thr0
    
    ds.close()
    
    # if fn > 13:
    #     ds = nc.Dataset(fp+f"pp/dyn_{fn:06d}.nc")
    #     thrpert = ds.variables['thrpert'][:].data[iz,idy,ix]
    #     thr0 = np.tile(ds.variables['thr0'][:].data[iz], (thrpert.shape[1], 1)).transpose()
    #     B = 9.8 * (thrpert/thr0)
    #     del thr0
    #     ds.close()
    
    # ds = nc.Dataset(fp+f"base/cm1out_000001.nc")
    # prspert = (prs - ds.variables['prs0'][:].data[0,iz,idy,ix])/100
    # del prs
    # ds.close()
    
    
    
    if new_pkl:
        dbfile = open(ip+f"cross_sections/MV1_{mv_time}min/cross{stime/60:.0f}.pkl", 'wb')
        dat = {'time':stime, 'w':winterp, 'prspert':prspert, 'thrpert':thrpert, 'B':B}
        pickle.dump(dat, dbfile)
        dbfile.close()
        del winterp,thrpert,prspert,B
        
    if update_pkl:
        dbfile = open(ip+f"cross_sections/MV1_{mv_time}min/cross{stime/60:.0f}.pkl", 'rb')
        dat = pickle.load(dbfile)
        dbfile.close()
        
        dat.update({'OW':OW, 'zvort':zvort})
        dbfile = open(ip+f"cross_sections/MV1_{mv_time}min/cross{stime/60:.0f}.pkl", 'wb')
        pickle.dump(dat, dbfile)
        dbfile.close()
        del OW,zvort
    
    
#%% Animate w, thrpert, prspert cross sections

w_lims = [-20,20]
thr_lims = [-15,15]
pp_lims = [-4,4]
# OW_lims = [-0.01,0]
zlims = [0,3.5]

ds = nc.Dataset(fp+fname)
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
z = ds.variables['z'][:].data
ds.close()

iz5 = np.where(z >= 5)[0][0]
iz = slice(0,iz5+1)

w_levs = np.linspace(w_lims[0], w_lims[1], (w_lims[1]-w_lims[0])+1)
thr_levs = np.linspace(thr_lims[0], thr_lims[1], (thr_lims[1]-thr_lims[0])+1)
pp_levs = np.linspace(pp_lims[0], pp_lims[1], 2*(pp_lims[1]-pp_lims[0])+1)

ti = np.where(ptime == stime)[0][0]
y_median = np.median(y_ml[ti,:])
idy = (np.abs(yh - y_median)).argmin()
# x_median = np.median(x_ml[ti,:])
# idx = (np.abs(xh - x_median)).argmin()

xlims = [-30,10] #[-55,25]

ix = slice(np.where(xh >= xlims[0])[0][0], np.where(xh >= xlims[1])[0][1])
iy = slice(idy-4,idy+5)

qix = int(np.round(len(xh[(xh>=xlims[0]) & (xh<=xlims[1])])/60)*2) # scale wind vector spacing to number of grid pts
qiz = 3


if True:
    dbfile = open(ip+f"cross_sections/MV1_{mv_time}min/cross180.pkl", 'rb')
    cs = pickle.load(dbfile)
    stime = cs['time']
    w_cross = cs['w']
    prs_cross = cs['prspert']
    thr_cross = cs['thrpert']
    OW_cross = cs['OW']
    zvort_cross = cs['zvort']
    B_cross = cs['B']
    dbfile.close()
    
    fig,ax = plt.subplots(1,1,figsize=(9,4))
    plot_contourf(xh[ix], z[iz], w_cross, 'w', ax, levels=w_levs, datalims=w_lims, xlims=xlims, ylims=zlims)
    ax.contour(xh[ix], z[iz], thr_cross, levels=[-5,-1], colors='purple', linewidths=[2,1], linestyles='-')
    # ax.contour(xh[ix], z[iz], OW, levels=[-0.01,-0.0025,-0.001], colors='hotpink', linewidths=1, linestyles='-')
    ax.scatter(x_ml[0,:], z_ml[0,:], s=10, color='k', marker='.')
    ax.set_xlabel('x (km)')
    ax.set_ylabel('z (km)')
    ax.set_title(f"Mid-level parcels in the MV at {mv_time} min \n t={stime/60:.0f} min")
    plt.tight_layout()
    
    def animate_w_cross(i):
        global ax
        ax.clear()
        
        dbfile = open(ip+f"cross_sections/MV1_{mv_time}min/cross{i+180}.pkl", 'rb')
        cs = pickle.load(dbfile)
        stime = cs['time']
        w_cross = cs['w']
        prs_cross = cs['prspert']
        thr_cross = cs['thrpert']
        OW_cross = cs['OW']
        zvort_cross = cs['zvort']
        B_cross = cs['B']
        dbfile.close()
        ti = np.where(ptime == stime)[0][0]
        
        plot_contourf(xh[ix], z[iz], w_cross, 'w', ax, levels=w_levs, datalims=w_lims, xlims=xlims, ylims=zlims, cbar=False)
        ax.contour(xh[ix], z[iz], thr_cross, levels=[-5,-1], colors='purple', linewidths=[2,1], linestyles='-')
        # ax.contour(xh[ix], z[iz], OW, levels=[-0.01,-0.0025,-0.001], colors='hotpink', linewidths=1, linestyles='-')
        ax.scatter(x_ml[ti,:], z_ml[ti,:], s=10, color='k', marker='.')
        ax.set_xlabel('x (km)')
        ax.set_ylabel('z (km)')
        ax.set_title(f"Mid-level parcels in the MV at {mv_time} min \n t={stime/60:.0f} min")
        plt.tight_layout()
    
    figsave = True
    
    anim = FuncAnimation(fig, animate_w_cross, frames=fnum-12, interval=100, repeat=False, blit=False)
    if figsave:
        anim.save(ip+f"MV1_w_cross_{mv_time}min.gif", dpi=300)
    plt.show()



if True:
    dbfile = open(ip+f"cross_sections/MV1_{mv_time}min/cross180.pkl", 'rb')
    cs = pickle.load(dbfile)
    stime = cs['time']
    w_cross = cs['w']
    prs_cross = cs['prspert']
    thr_cross = cs['thrpert']
    OW_cross = cs['OW']
    zvort_cross = cs['zvort']
    B_cross = cs['B']
    dbfile.close()
    
    fig,ax = plt.subplots(1,1,figsize=(9,4))
    plot_contourf(xh[ix], z[iz], thr_cross, 'thrpert', ax, levels=thr_levs, datalims=thr_lims, xlims=xlims, ylims=zlims)
    ax.contour(xh[ix], z[iz], prs_cross, levels=[1,2,3,4,5], colors='r', linewidths=1)
    ax.contour(xh[ix], z[iz], prs_cross, levels=[-5,-4,-3,-2,-1], colors='b', linewidths=1, linestyles='-')
    # ax.contour(xh[ix], z[iz], OW, levels=[-0.01,-0.0025,-0.001], colors='hotpink', linewidths=1, linestyles='-')
    # ax.quiver(xh[ix][::qix], z[iz][::qiz], 0*w_cross[::qiz,::qix], w_cross[::qiz,::qix], scale=250, width=0.002, pivot='tail')
    ax.scatter(x_ml[0,:], z_ml[0,:], s=10, color='k', marker='.')
    ax.set_xlabel('x (km)')
    ax.set_ylabel('z (km)')
    ax.set_title(f"Mid-level parcels in the MV at {mv_time} min \n t={stime/60:.0f} min")
    plt.tight_layout()
    
    def animate_thr_cross(i):
        global ax
        ax.clear()
        
        dbfile = open(ip+f"cross_sections/MV1_{mv_time}min/cross{i+180}.pkl", 'rb')
        cs = pickle.load(dbfile)
        stime = cs['time']
        w_cross = cs['w']
        prs_cross = cs['prspert']
        thr_cross = cs['thrpert']
        OW_cross = cs['OW']
        zvort_cross = cs['zvort']
        B_cross = cs['B']
        dbfile.close()
        ti = np.where(ptime == stime)[0][0]
        
        plot_contourf(xh[ix], z[iz], thr_cross, 'thrpert', ax, levels=thr_levs, datalims=thr_lims, xlims=xlims, ylims=zlims, cbar=False)
        ax.contour(xh[ix], z[iz], prs_cross, levels=[1,2,3,4,5], colors='r', linewidths=1)
        ax.contour(xh[ix], z[iz], prs_cross, levels=[-5,-4,-3,-2,-1], colors='b', linewidths=1, linestyles='-')
        # ax.contour(xh[ix], z[iz], OW, levels=[-0.01,-0.0025,-0.001], colors='hotpink', linewidths=1, linestyles='-')
        # ax.quiver(xh[ix][::qix], z[iz][::qiz], 0*w_cross[::qiz,::qix], w_cross[::qiz,::qix], scale=250, width=0.002, pivot='tail')
        ax.scatter(x_ml[ti,:], z_ml[ti,:], s=10, color='k', marker='.')
        ax.set_xlabel('x (km)')
        ax.set_ylabel('z (km)')
        ax.set_title(f"Mid-level parcels in the MV at {mv_time} min \n t={stime/60:.0f} min")
        plt.tight_layout()
    
    figsave = True
    
    anim = FuncAnimation(fig, animate_thr_cross, frames=fnum-12, interval=100, repeat=False, blit=False)
    if figsave:
        anim.save(ip+f"MV1_thr_cross_{mv_time}min.gif", dpi=300)
    plt.show()


if True:
    dbfile = open(ip+f"cross_sections/MV1_{mv_time}min/cross180.pkl", 'rb')
    cs = pickle.load(dbfile)
    stime = cs['time']
    w_cross = cs['w']
    prs_cross = cs['prspert']
    thr_cross = cs['thrpert']
    OW_cross = cs['OW']
    zvort_cross = cs['zvort']
    B_cross = cs['B']
    dbfile.close()
    
    
    fig,ax = plt.subplots(1,1,figsize=(9,4))
    plot_contourf(xh[ix], z[iz], prs_cross, 'prspert', ax, levels=pp_levs, datalims=pp_lims, xlims=xlims, ylims=zlims)
    ax.contour(xh[ix], z[iz], B_cross, levels=[-0.2,-0.1,0.1,0.2], colors='purple', linewidths=1)
    ax.contour(xh[ix], z[iz], OW_cross, levels=[-0.01,-0.0025,-0.001], colors='k', linewidths=1, linestyles='-')
    ax.scatter(x_ml[0,:], z_ml[0,:], s=10, color='k', marker='.')
    ax.set_xlabel('x (km)')
    ax.set_ylabel('z (km)')
    ax.set_title(f"Mid-level parcels in the MV at {mv_time} min \n t={stime/60:.0f} min")
    plt.tight_layout()
    
    def animate_prs_cross(i):
        global ax
        ax.clear()
        
        dbfile = open(ip+f"cross_sections/MV1_{mv_time}min/cross{i+180}.pkl", 'rb')
        cs = pickle.load(dbfile)
        stime = cs['time']
        w_cross = cs['w']
        prs_cross = cs['prspert']
        thr_cross = cs['thrpert']
        OW_cross = cs['OW']
        zvort_cross = cs['zvort']
        B_cross = cs['B']
        dbfile.close()
        ti = np.where(ptime == stime)[0][0]
        
        plot_contourf(xh[ix], z[iz], prs_cross, 'prspert', ax, levels=pp_levs, datalims=pp_lims, xlims=xlims, ylims=zlims, cbar=False)
        ax.contour(xh[ix], z[iz], B_cross, levels=[-0.2,-0.1,0.1,0.2], colors='purple', linewidths=1)
        ax.contour(xh[ix], z[iz], OW_cross, levels=[-0.01,-0.0025,-0.001], colors='k', linewidths=1, linestyles='-')
        ax.scatter(x_ml[ti,:], z_ml[ti,:], s=10, color='k', marker='.')
        ax.set_xlabel('x (km)')
        ax.set_ylabel('z (km)')
        ax.set_title(f"Mid-level parcels in the MV at {mv_time} min \n t={stime/60:.0f} min")
        plt.tight_layout()
    
    figsave = True
    
    anim = FuncAnimation(fig, animate_prs_cross, frames=fnum-12, interval=150, repeat=False, blit=False)
    if figsave:
        anim.save(ip+f"MV1_prs_cross_{mv_time}min.gif", dpi=300)
    plt.show()


#%% Plan views

from matplotlib.animation import FuncAnimation

new_pkl = False
update_pkl = False

xl = [-50,10]
yl = [-140,-80]

for fn in np.arange(13,fnum+1):
    if fn == 13:
        fname = f"base/cm1out_{fn:06d}.nc"
    else:
        fname = f"cm1out_{fn:06d}.nc"
    
    print(f"cm1out_{fn:06d}")
    # Read output file
    ds = nc.Dataset(fp+fname)
    stime = ds.variables['time'][:].data[0]
    xh = ds.variables['xh'][:].data
    yh = ds.variables['yh'][:].data
    z = ds.variables['z'][:].data
    
    iz500 = np.where(z >= 0.5)[0][0]
    iz = slice(0,iz500+1)
    
    ti = np.where(ptime == stime)[0][0]
    y_median = np.median(y_ml[ti,:])
    idy = (np.abs(yh - y_median)).argmin()
    x_median = np.median(x_ml[ti,:])
    idx = (np.abs(xh - x_median)).argmin()
    
    xlims = xl
    ylims = yl + 1.2*(fn-13)
    
    ix = slice(np.where(xh >= xlims[0])[0][0], np.where(xh >= xlims[1])[0][1])
    iy = slice(np.where(yh >= ylims[0])[0][0], np.where(yh >= ylims[1])[0][1])
    
    # prs = ds.variables['prs'][:].data[0,iz,iy,ix]
    # winterp = ds.variables['winterp'][:].data[0,iz500,iy,ix]
    
    uinterp = ds.variables['uinterp'][:].data[0,iz,iy,ix]
    vinterp = ds.variables['vinterp'][:].data[0,iz,iy,ix]
    zvort = ds.variables['zvort'][:].data[0,iz,iy,ix]
    S_N = np.gradient(uinterp, xh[ix]*1000, axis=2) - np.gradient(vinterp, yh[iy]*1000, axis=1)
    S_S = np.gradient(vinterp, xh[ix]*1000, axis=2) + np.gradient(uinterp, yh[iy]*1000, axis=1)
    OW = S_N**2 + S_S**2 - zvort**2
    # OW = OW[iz500,:,:]
    # zvort = zvort[iz500,:,:]
    OW = np.min(OW, axis=0)
    zvort = np.max(zvort, axis=0)
    del S_N,S_S,uinterp,vinterp
    
    # if fn == 13:
    #     thr = ds.variables['th'][:].data[0,iz,iy,ix] * (1 + 0.61*ds.variables['qv'][:].data[0,iz,iy,ix] - 
    #                 (ds.variables['qc'][:].data[0,iz,iy,ix] + ds.variables['qr'][:].data[0,iz,iy,ix] + 
    #                   ds.variables['qi'][:].data[0,iz,iy,ix] + ds.variables['qs'][:].data[0,iz,iy,ix] + 
    #                   ds.variables['qg'][:].data[0,iz,iy,ix] + ds.variables['qhl'][:].data[0,iz,iy,ix]))
    #     thr0 = ds.variables['th0'][:].data[0,iz,iy,ix] * (1 + 0.61*ds.variables['qv0'][:].data[0,iz,iy,ix])
    #     thrpert = thr - thr0
    #     B = np.mean(9.8 * (thrpert/thr0), axis=0)
    #     del thr,thr0
    #     thrpert = thrpert[0,:,:]
    
    ds.close()
    
    # if fn > 13:
    #     ds = nc.Dataset(fp+f"pp/dyn_{fn:06d}.nc")
    #     thrpert = ds.variables['thrpert'][:].data[iz,iy,ix]
    #     thr0 = np.moveaxis(np.tile(ds.variables['thr0'][:].data[iz], (thrpert.shape[1], thrpert.shape[2], 1)), -1, 0)
    #     B = np.mean(9.8 * (thrpert/thr0), axis=0)
    #     del thr0
    #     thrpert = thrpert[0,:,:]
    #     ds.close()
    
    # ds = nc.Dataset(fp+f"base/cm1out_000001.nc")
    # prspert = np.mean((prs - ds.variables['prs0'][:].data[0,iz,iy,ix])/100, axis=0)
    # del prs
    # ds.close()
    
    
    
    if new_pkl:
        dbfile = open(ip+f"cross_sections/MV1_{mv_time}min/plan{stime/60:.0f}.pkl", 'wb')
        dat = {'time':stime, 'xlims':xlims, 'ylims':ylims, 'w':winterp, 'prspert':prspert, 'thrpert':thrpert, 'B':B}
        pickle.dump(dat, dbfile)
        dbfile.close()
        del winterp,thrpert,prspert,B
        
    if update_pkl:
        dbfile = open(ip+f"cross_sections/MV1_{mv_time}min/plan{stime/60:.0f}.pkl", 'rb')
        dat = pickle.load(dbfile)
        dbfile.close()
        
        dat.update({'OW':OW, 'zvort':zvort})
        dbfile = open(ip+f"cross_sections/MV1_{mv_time}min/plan{stime/60:.0f}.pkl", 'wb')
        pickle.dump(dat, dbfile)
        dbfile.close()
        del OW,zvort

#%% Animate plan views

if True:
    dbfile = open(ip+f"cross_sections/MV1_{mv_time}min/plan180.pkl", 'rb')
    cs = pickle.load(dbfile)
    stime = cs['time']
    xlims = cs['xlims']
    ylims = cs['ylims']
    w_plan = cs['w'] # 500m
    prs_plan = cs['prspert'] # mean 0-500m
    thr_plan = cs['thrpert'] # sfc
    OW_plan = cs['OW'] # min 0-500m
    zvort_plan = cs['zvort'] # max 0-500m
    B_plan = cs['B'] # mean 0-500m
    dbfile.close()
    
    ix = slice(np.where(xh >= xlims[0])[0][0], np.where(xh >= xlims[1])[0][1])
    iy = slice(np.where(yh >= ylims[0])[0][0], np.where(yh >= ylims[1])[0][1])
    
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    plot_contourf(xh[ix], yh[iy], w_plan, 'w', ax, levels=w_levs, datalims=w_lims, xlims=xlims, ylims=ylims)
    ax.contour(xh[ix], yh[iy], thr_plan, levels=[-5,-1], colors='purple', linewidths=1, linestyles=['-','--'])
    # ax.contour(xh[ix], yh[iy], OW_plan, levels=[-0.01,-0.0025,-0.001], colors='hotpink', linewidths=1, linestyles='-')
    # ax.contour(xh[ix], yh[iy], prs_plan, levels=[-4,-3,-2,-1,1,2,3,4], colors='purple', linewidths=1)
    ax.plot([-30,10], [np.median(y_ml[0,:]),np.median(y_ml[0,:])], 'k', linewidth=0.75)
    ax.scatter(x_ll[0,:], y_ll[0,:], s=20, color='gray', marker='.')
    ax.scatter(x_sc[0,:], y_sc[0,:], s=20, color='gray', marker='.')
    ax.scatter(x_qc[0,:], y_qc[0,:], s=20, color='gray', marker='.')
    ax.scatter(x_ml[0,:], y_ml[0,:], s=50, color='k', marker='.')
    ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')
    ax.set_title(f"Mid-level parcels in the MV at {mv_time} min \n t={stime/60:.0f} min")
    plt.tight_layout()
    
    def animate_w_plan(i):
        global ax
        ax.clear()
        
        dbfile = open(ip+f"cross_sections/MV1_{mv_time}min/plan{i+180}.pkl", 'rb')
        cs = pickle.load(dbfile)
        stime = cs['time']
        xlims = cs['xlims']
        ylims = cs['ylims']
        w_plan = cs['w']
        prs_plan = cs['prspert']
        thr_plan = cs['thrpert']
        OW_plan = cs['OW']
        zvort_plan = cs['zvort']
        B_plan = cs['B']
        dbfile.close()
        
        ti = np.where(ptime == stime)[0][0]
        ix = slice(np.where(xh >= xlims[0])[0][0], np.where(xh >= xlims[1])[0][1])
        iy = slice(np.where(yh >= ylims[0])[0][0], np.where(yh >= ylims[1])[0][1])
        
        plot_contourf(xh[ix], yh[iy], w_plan, 'w', ax, levels=w_levs, datalims=w_lims, xlims=xlims, ylims=ylims, cbar=False)
        ax.contour(xh[ix], yh[iy], thr_plan, levels=[-5,-1], colors='purple', linewidths=1, linestyles=['-','--'])
        # ax.contour(xh[ix], yh[iy], OW_plan, levels=[-0.01,-0.0025,-0.001], colors='hotpink', linewidths=1, linestyles='-')
        # ax.contour(xh[ix], yh[iy], prs_plan, levels=[-4,-3,-2,-1,1,2,3,4], colors='purple', linewidths=1)
        ax.plot([-30,10], [np.median(y_ml[ti,:]),np.median(y_ml[ti,:])], 'k', linewidth=0.75)
        ax.scatter(x_ll[ti,:], y_ll[ti,:], s=20, color='gray', marker='.')
        ax.scatter(x_sc[ti,:], y_sc[ti,:], s=20, color='gray', marker='.')
        ax.scatter(x_qc[ti,:], y_qc[ti,:], s=20, color='gray', marker='.')
        ax.scatter(x_ml[ti,:], y_ml[ti,:], s=50, color='k', marker='.')
        ax.set_xlabel('x (km)')
        ax.set_ylabel('y (km)')
        ax.set_title(f"Mid-level parcels in the MV at {mv_time} min \n t={stime/60:.0f} min")
        plt.tight_layout()
    
    figsave = True
    
    anim = FuncAnimation(fig, animate_w_plan, frames=fnum-12, interval=100, repeat=False, blit=False)
    if figsave:
        anim.save(ip+f"MV1_w_plan_{mv_time}min.gif", dpi=300)
    plt.show()
    
    
if True:
    dbfile = open(ip+f"cross_sections/MV1_{mv_time}min/plan180.pkl", 'rb')
    cs = pickle.load(dbfile)
    stime = cs['time']
    xlims = cs['xlims']
    ylims = cs['ylims']
    w_plan = cs['w'] # 500m
    prs_plan = cs['prspert'] # mean 0-500m
    thr_plan = cs['thrpert'] # sfc
    OW_plan = cs['OW'] # min 0-500m
    zvort_plan = cs['zvort'] # max 0-500m
    B_plan = cs['B'] # mean 0-500m
    dbfile.close()
    
    ix = slice(np.where(xh >= xlims[0])[0][0], np.where(xh >= xlims[1])[0][1])
    iy = slice(np.where(yh >= ylims[0])[0][0], np.where(yh >= ylims[1])[0][1])
    
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    plot_contourf(xh[ix], yh[iy], thr_plan, 'thrpert', ax, levels=thr_levs, datalims=thr_lims, xlims=xlims, ylims=ylims)
    # ax.contour(xh[ix], yh[iy], B_plan, levels=[1,2,3,4,5], colors='r', linewidths=1)
    # ax.contour(xh[ix], yh[iy], OW_plan, levels=[-0.01,-0.0025,-0.001], colors='hotpink', linewidths=1, linestyles='-')
    ax.contour(xh[ix], yh[iy], w_plan, levels=[-20,-15,-10,-5,5,10,15,20], colors='purple', linewidths=1)
    ax.plot([-30,10], [np.median(y_ml[0,:]),np.median(y_ml[0,:])], 'k', linewidth=0.75)
    ax.scatter(x_ll[0,:], y_ll[0,:], s=20, color='gray', marker='.')
    ax.scatter(x_sc[0,:], y_sc[0,:], s=20, color='gray', marker='.')
    ax.scatter(x_qc[0,:], y_qc[0,:], s=20, color='gray', marker='.')
    ax.scatter(x_ml[0,:], y_ml[0,:], s=50, color='k', marker='.')
    ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')
    ax.set_title(f"Mid-level parcels in the MV at {mv_time} min \n t={stime/60:.0f} min")
    plt.tight_layout()
    
    def animate_thr_plan(i):
        global ax
        ax.clear()
        
        dbfile = open(ip+f"cross_sections/MV1_{mv_time}min/plan{i+180}.pkl", 'rb')
        cs = pickle.load(dbfile)
        stime = cs['time']
        xlims = cs['xlims']
        ylims = cs['ylims']
        w_plan = cs['w']
        prs_plan = cs['prspert']
        thr_plan = cs['thrpert']
        OW_plan = cs['OW']
        zvort_plan = cs['zvort']
        B_plan = cs['B']
        dbfile.close()
        
        ti = np.where(ptime == stime)[0][0]
        ix = slice(np.where(xh >= xlims[0])[0][0], np.where(xh >= xlims[1])[0][1])
        iy = slice(np.where(yh >= ylims[0])[0][0], np.where(yh >= ylims[1])[0][1])
        
        plot_contourf(xh[ix], yh[iy], thr_plan, 'thrpert', ax, levels=thr_levs, datalims=thr_lims, xlims=xlims, ylims=ylims, cbar=False)
        # ax.contour(xh[ix], yh[iy], B_plan, levels=[1,2,3,4,5], colors='r', linewidths=1)
        # ax.contour(xh[ix], yh[iy], OW_plan, levels=[-0.01,-0.0025,-0.001], colors='hotpink', linewidths=1, linestyles='-')
        ax.contour(xh[ix], yh[iy], w_plan, levels=[-20,-15,-10,-5,5,10,15,20], colors='purple', linewidths=1)
        ax.plot([-30,10], [np.median(y_ml[ti,:]),np.median(y_ml[ti,:])], 'k', linewidth=0.75)
        ax.scatter(x_ll[ti,:], y_ll[ti,:], s=20, color='gray', marker='.')
        ax.scatter(x_sc[ti,:], y_sc[ti,:], s=20, color='gray', marker='.')
        ax.scatter(x_qc[ti,:], y_qc[ti,:], s=20, color='gray', marker='.')
        ax.scatter(x_ml[ti,:], y_ml[ti,:], s=50, color='k', marker='.')
        ax.set_xlabel('x (km)')
        ax.set_ylabel('y (km)')
        ax.set_title(f"Mid-level parcels in the MV at {mv_time} min \n t={stime/60:.0f} min")
        plt.tight_layout()
    
    figsave = True
    
    anim = FuncAnimation(fig, animate_thr_plan, frames=fnum-12, interval=100, repeat=False, blit=False)
    if figsave:
        anim.save(ip+f"MV1_thr_plan_{mv_time}min.gif", dpi=300)
    plt.show()
    
    
    
    
if True:
    dbfile = open(ip+f"cross_sections/MV1_{mv_time}min/plan180.pkl", 'rb')
    cs = pickle.load(dbfile)
    stime = cs['time']
    xlims = cs['xlims']
    ylims = cs['ylims']
    w_plan = cs['w'] # 500m
    prs_plan = cs['prspert'] # mean 0-500m
    thr_plan = cs['thrpert'] # sfc
    OW_plan = cs['OW'] # min 0-500m
    zvort_plan = cs['zvort'] # max 0-500m
    B_plan = cs['B'] # mean 0-500m
    dbfile.close()
    
    ix = slice(np.where(xh >= xlims[0])[0][0], np.where(xh >= xlims[1])[0][1])
    iy = slice(np.where(yh >= ylims[0])[0][0], np.where(yh >= ylims[1])[0][1])
    
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    plot_contourf(xh[ix], yh[iy], prs_plan, 'prspert', ax, levels=pp_levs, datalims=pp_lims, xlims=xlims, ylims=ylims)
    # ax.contour(xh[ix], yh[iy], B_plan, levels=[1,2,3,4,5], colors='r', linewidths=1)
    # ax.contour(xh[ix], yh[iy], OW_plan, levels=[-0.01,-0.0025,-0.001], colors='hotpink', linewidths=1, linestyles='-')
    ax.contour(xh[ix], yh[iy], w_plan, levels=[-20,-15,-10,-5,5,10,15,20], colors='purple', linewidths=1)
    ax.plot([-30,10], [np.median(y_ml[0,:]),np.median(y_ml[0,:])], 'k', linewidth=0.75)
    ax.scatter(x_ll[0,:], y_ll[0,:], s=20, color='gray', marker='.')
    ax.scatter(x_sc[0,:], y_sc[0,:], s=20, color='gray', marker='.')
    ax.scatter(x_qc[0,:], y_qc[0,:], s=20, color='gray', marker='.')
    ax.scatter(x_ml[0,:], y_ml[0,:], s=50, color='k', marker='.')
    ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')
    ax.set_title(f"Mid-level parcels in the MV at {mv_time} min \n t={stime/60:.0f} min")
    plt.tight_layout()
    
    def animate_prs_plan(i):
        global ax
        ax.clear()
        
        dbfile = open(ip+f"cross_sections/MV1_{mv_time}min/plan{i+180}.pkl", 'rb')
        cs = pickle.load(dbfile)
        stime = cs['time']
        xlims = cs['xlims']
        ylims = cs['ylims']
        w_plan = cs['w']
        prs_plan = cs['prspert']
        thr_plan = cs['thrpert']
        OW_plan = cs['OW']
        zvort_plan = cs['zvort']
        B_plan = cs['B']
        dbfile.close()
        
        ti = np.where(ptime == stime)[0][0]
        ix = slice(np.where(xh >= xlims[0])[0][0], np.where(xh >= xlims[1])[0][1])
        iy = slice(np.where(yh >= ylims[0])[0][0], np.where(yh >= ylims[1])[0][1])
        
        plot_contourf(xh[ix], yh[iy], prs_plan, 'prspert', ax, levels=pp_levs, datalims=pp_lims, xlims=xlims, ylims=ylims, cbar=False)
        # ax.contour(xh[ix], yh[iy], B_plan, levels=[1,2,3,4,5], colors='r', linewidths=1)
        # ax.contour(xh[ix], yh[iy], OW_plan, levels=[-0.01,-0.0025,-0.001], colors='hotpink', linewidths=1, linestyles='-')
        ax.contour(xh[ix], yh[iy], w_plan, levels=[-20,-15,-10,-5,5,10,15,20], colors='purple', linewidths=1)
        ax.plot([-30,10], [np.median(y_ml[ti,:]),np.median(y_ml[ti,:])], 'k', linewidth=0.75)
        ax.scatter(x_ll[ti,:], y_ll[ti,:], s=20, color='gray', marker='.')
        ax.scatter(x_sc[ti,:], y_sc[ti,:], s=20, color='gray', marker='.')
        ax.scatter(x_qc[ti,:], y_qc[ti,:], s=20, color='gray', marker='.')
        ax.scatter(x_ml[ti,:], y_ml[ti,:], s=50, color='k', marker='.')
        ax.set_xlabel('x (km)')
        ax.set_ylabel('y (km)')
        ax.set_title(f"Mid-level parcels in the MV at {mv_time} min \n t={stime/60:.0f} min")
        plt.tight_layout()
    
    figsave = True
    
    anim = FuncAnimation(fig, animate_prs_plan, frames=fnum-12, interval=100, repeat=False, blit=False)
    if figsave:
        anim.save(ip+f"MV1_prs_plan_{mv_time}min.gif", dpi=300)
    plt.show()

