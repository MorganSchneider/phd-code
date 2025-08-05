#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 15:53:42 2024

@author: morgan.schneider

- PAPER FIGURES -
traj_thr_cross (single time)
traj_thr_plan (single time)
traj_vort2D_prcl-height (single time)

- Other figures -
Cross section animations (good for talks)
Plan view animations (good for talks)
"""

####################
### Load modules ###
####################

from CM1utils import *

#%% Load parcel data

mv_time = 220
fnum = 53

fp = '/Volumes/Promise_Pegasus_70TB/merger/merger-125m/'
ip = '/Users/morgan.schneider/Documents/merger/merger-125m/'

# Read parcel data
ds = nc.Dataset(fp+'cm1out_pdata.nc')
ptime = ds.variables['time'][:].data
ds.close()

dbfile = open('/Users/morgan.schneider/Documents/merger/traj_MV1.pkl', 'rb')
traj = pickle.load(dbfile)
pids_mv = traj[f"{mv_time}min"]['pids']
x_mv = traj[f"{mv_time}min"]['x']/1000
y_mv = traj[f"{mv_time}min"]['y']/1000
z_mv = traj[f"{mv_time}min"]['z']/1000
u_mv = traj[f"{mv_time}min"]['u']
v_mv = traj[f"{mv_time}min"]['v']
w_mv = traj[f"{mv_time}min"]['w']
zvort_mv = traj[f"{mv_time}min"]['zvort']
dbfile.close()

dbfile = open(ip+f"traj_clusters/traj_clusters_{mv_time}min_v2.pkl", 'rb')
c = pickle.load(dbfile)
cc = c['mv1']
dbfile.close()


# Mid-level source
x_ml = x_mv[:,(cc==1)]
y_ml = y_mv[:,(cc==1)]
z_ml = z_mv[:,(cc==1)]
u_ml = u_mv[:,(cc==1)]
v_ml = v_mv[:,(cc==1)]
w_ml = w_mv[:,(cc==1)]
zvort_ml = zvort_mv[:,(cc==1)]


#%% Load model grid data for one time (useless)

fn = 51 # use 41 (208 min) for 210 min, 56 (223 min) for 225 min, 51 (218 min) for 220 min

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

ti1 = np.where(ptime == stime)[0][0]
y_median = np.median(y_ml[ti1,:])
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



#%% Cross sections save to pkl

ip = '/Users/morgan.schneider/Documents/merger/merger-125m/'

new_pkl = True
update_pkl = False

ds = nc.Dataset(fp+f"base/cm1out_000001.nc")
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
z = ds.variables['z'][:].data

iz5 = np.where(z >= 5)[0][0]
iz = slice(0,iz5+1)
xlims = [-30,10] #[-55,25]
ix = slice(np.where(xh >= xlims[0])[0][0], np.where(xh >= xlims[1])[0][1])

prs0 = ds.variables['prs0'][:].data[0,iz,:,ix]
ds.close()

# ds = nc.Dataset(fp+f"pp/dyn_000014.nc")
# thr0 = ds.variables['thr0'][:].data[iz]
# ds.close()

for fn in np.arange(54,74): # fnum+1:74 for after mv time
    if fn == 13:
        fname = f"base/cm1out_{fn:06d}.nc"
    else:
        fname = f"cm1out_{fn:06d}.nc"
    
    print(f"cm1out_{fn:06d}")
    # Read output file
    ds = nc.Dataset(fp+fname)
    stime = ds.variables['time'][:].data[0]
    
    ti = np.where(ptime == stime)[0][0]
    y_median = np.median(y_ml[ti,:])
    idy = (np.abs(yh - y_median)).argmin()
    # x_median = np.median(x_ml[ti,:])
    # idx = (np.abs(xh - x_median)).argmin()
    
    iy = slice(idy-4,idy+5)
    
    prspert = (ds.variables['prs'][:].data[0,iz,idy,ix] - prs0[:,idy,:])/100
    # winterp = ds.variables['winterp'][:].data[0,iz,idy,ix]
    
    # uinterp = ds.variables['uinterp'][:].data[0,iz,iy,ix]
    # vinterp = ds.variables['vinterp'][:].data[0,iz,iy,ix]
    # zvort = ds.variables['zvort'][:].data[0,iz,idy,ix]
    # S_N = np.gradient(uinterp, xh[ix]*1000, axis=2) - np.gradient(vinterp, yh[iy]*1000, axis=1)
    # S_S = np.gradient(vinterp, xh[ix]*1000, axis=2) + np.gradient(uinterp, yh[iy]*1000, axis=1)
    # # OW = np.min((S_N**2 + S_S**2 - zvort**2), axis=1)
    # OW = S_N[:,4,:]**2 + S_S[:,4,:]**2 - zvort**2
    # del S_N,S_S,uinterp,vinterp
    
    if fn == 13:
        thr = ds.variables['th'][:].data[0,iz,idy,ix] * (1 + 0.61*ds.variables['qv'][:].data[0,iz,idy,ix] - 
                    (ds.variables['qc'][:].data[0,iz,idy,ix] + ds.variables['qr'][:].data[0,iz,idy,ix] + 
                      ds.variables['qi'][:].data[0,iz,idy,ix] + ds.variables['qs'][:].data[0,iz,idy,ix] + 
                      ds.variables['qg'][:].data[0,iz,idy,ix] + ds.variables['qhl'][:].data[0,iz,idy,ix]))
        thr0 = ds.variables['th0'][:].data[0,iz,idy,ix] * (1 + 0.61*ds.variables['qv0'][:].data[0,iz,idy,ix])
        thrpert = thr - thr0
        # B = 9.8 * (thrpert/thr0)
        del thr
    
    ds.close()
    
    if fn > 13:
        ds = nc.Dataset(fp+f"pp/dyn_{fn:06d}.nc")
        thrpert = ds.variables['thrpert'][:].data[iz,idy,ix]
        # thr0 = np.tile(ds.variables['thr0'][:].data[iz], (thrpert.shape[1], 1)).transpose()
        # B = 9.8 * (thrpert/thr0)
        # del thr0
        ds.close()
    
    dat = {'time':stime, 'prspert':prspert, 'thrpert':thrpert}
    fname = ip+f"cross_sections/MV1_{mv_time}min/cross{stime/60:.0f}.pkl"
    save_to_pickle(dat, fname)
    del thrpert,prspert
    
    
    if update_pkl:
        dbfile = open(ip+f"cross_sections/MV1_{mv_time}min/cross{stime/60:.0f}.pkl", 'rb')
        dat = pickle.load(dbfile)
        dbfile.close()
        
        dat.update({'w':winterp, 'B':B, 'OW':OW, 'zvort':zvort})
        dbfile = open(ip+f"cross_sections/MV1_{mv_time}min/cross{stime/60:.0f}.pkl", 'wb')
        pickle.dump(dat, dbfile)
        dbfile.close()
        del winterp,B,OW,zvort


#%% Plan views save to pkl

ip = '/Users/morgan.schneider/Documents/merger/merger-125m/'

new_pkl = False
update_pkl = False

xl = [-40,20]
yl = [-130,-70]

ds = nc.Dataset(fp+f"base/cm1out_000001.nc")
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
z = ds.variables['z'][:].data

iz500 = np.where(z >= 0.5)[0][0]
iz = slice(0,iz500+1)
ix = slice(np.where(xh >= xl[0])[0][0], np.where(xh >= xl[1])[0][1])

# prs0 = ds.variables['prs0'][:].data[0,iz,:,ix]
ds.close()

for fn in [13]:
    if fn == 13:
        fname = f"base/cm1out_{fn:06d}.nc"
    else:
        fname = f"cm1out_{fn:06d}.nc"
    
    print(f"cm1out_{fn:06d}")
    # Read output file
    ds = nc.Dataset(fp+fname)
    stime = ds.variables['time'][:].data[0]
    
    ti = np.where(ptime == stime)[0][0]
    y_median = np.median(y_ml[ti,:])
    idy = (np.abs(yh - y_median)).argmin()
    x_median = np.median(x_ml[ti,:])
    idx = (np.abs(xh - x_median)).argmin()
    
    xlims = xl
    # ylims = yl + 1.2*(fn-14)
    ylims = [yl[0] + 1.2*(fn-14), yl[1] + 1.2*(fn-14)]
    
    iy = slice(np.where(yh >= ylims[0])[0][0], np.where(yh >= ylims[1])[0][1])
    
    # prspert = ds.variables['prs'][:].data[0,iz,iy,ix] - prs0[:,iy,:]
    
    # winterp = ds.variables['winterp'][:].data[0,iz500,iy,ix]
    
    # uinterp = ds.variables['uinterp'][:].data[0,iz,iy,ix]
    # vinterp = ds.variables['vinterp'][:].data[0,iz,iy,ix]
    # zvort = ds.variables['zvort'][:].data[0,iz,iy,ix]
    # S_N = np.gradient(uinterp, xh[ix]*1000, axis=2) - np.gradient(vinterp, yh[iy]*1000, axis=1)
    # S_S = np.gradient(vinterp, xh[ix]*1000, axis=2) + np.gradient(uinterp, yh[iy]*1000, axis=1)
    # OW = S_N**2 + S_S**2 - zvort**2
    # # OW = OW[iz500,:,:]
    # # zvort = zvort[iz500,:,:]
    # OW = np.min(OW, axis=0)
    # zvort = np.max(zvort, axis=0)
    # del S_N,S_S,uinterp,vinterp
    
    if fn == 13:
        thr = ds.variables['th'][:].data[0,0,iy,ix] * (1 + 0.61*ds.variables['qv'][:].data[0,0,iy,ix] - 
                    (ds.variables['qc'][:].data[0,0,iy,ix] + ds.variables['qr'][:].data[0,0,iy,ix] + 
                      ds.variables['qi'][:].data[0,0,iy,ix] + ds.variables['qs'][:].data[0,0,iy,ix] + 
                      ds.variables['qg'][:].data[0,0,iy,ix] + ds.variables['qhl'][:].data[0,0,iy,ix]))
        thr0 = ds.variables['th0'][:].data[0,0,iy,ix] * (1 + 0.61*ds.variables['qv0'][:].data[0,0,iy,ix])
        thrpert = thr - thr0
        # B = np.mean(9.8 * (thrpert/thr0), axis=0)
        # del thr,thr0
    
    ds.close()
    
    if fn > 13:
        ds = nc.Dataset(fp+f"pp/dyn_{fn:06d}.nc")
        thrpert = ds.variables['thrpert'][:].data[0,iy,ix]
        # thr0 = np.moveaxis(np.tile(ds.variables['thr0'][:].data[iz], (thrpert.shape[1], thrpert.shape[2], 1)), -1, 0)
        # B = np.mean(9.8 * (thrpert/thr0), axis=0)
        # del thr0
        ds.close()
    
    
    dat = {'time':stime, 'xlims':xlims, 'ylims':ylims, 'thrpert':thrpert}
    fname = ip+f"cross_sections/MV1_{mv_time}min/plan{stime/60:.0f}.pkl"
    save_to_pickle(dat, fname)
    del thrpert
    
    
    if update_pkl:
        dbfile = open(ip+f"cross_sections/MV1_{mv_time}min/plan{stime/60:.0f}.pkl", 'rb')
        dat = pickle.load(dbfile)
        dbfile.close()
        
        dat.update({'w':winterp, 'B':B, 'OW':OW, 'zvort':zvort})
        dbfile = open(ip+f"cross_sections/MV1_{mv_time}min/plan{stime/60:.0f}.pkl", 'wb')
        pickle.dump(dat, dbfile)
        dbfile.close()
        del winterp,B,OW,zvort

#%% Plot cross sections ***OLD PAPER FIGS HERE***

mv_time = 220

if mv_time == 210:
    fn = 41
elif mv_time == 220:
    fn = 51
elif mv_time == 225:
    fn = 56


fp = '/Volumes/Promise_Pegasus_70TB/merger/merger-125m/'
ip = '/Users/morgan.schneider/Documents/merger/merger-125m/'

# Read parcel data
ds = nc.Dataset(fp+'cm1out_pdata.nc')
ptime = ds.variables['time'][:].data
ds.close()

dbfile = open('/Users/morgan.schneider/Documents/merger/traj_MV1.pkl', 'rb')
traj = pickle.load(dbfile)
pids_mv = traj[f"{mv_time}min"]['pids']
x_mv = traj[f"{mv_time}min"]['x']/1000
y_mv = traj[f"{mv_time}min"]['y']/1000
z_mv = traj[f"{mv_time}min"]['z']/1000
u_mv = traj[f"{mv_time}min"]['u']
v_mv = traj[f"{mv_time}min"]['v']
w_mv = traj[f"{mv_time}min"]['w']
zvort_mv = traj[f"{mv_time}min"]['zvort']
dbfile.close()

dbfile = open(ip+f"traj_clusters/traj_clusters_{mv_time}min_v2.pkl", 'rb')
c = pickle.load(dbfile)
cc = c['mv1']
dbfile.close()


# Mid-level source
x_ml = x_mv[:,(cc==1)]
y_ml = y_mv[:,(cc==1)]
z_ml = z_mv[:,(cc==1)]
u_ml = u_mv[:,(cc==1)]
v_ml = v_mv[:,(cc==1)]
w_ml = w_mv[:,(cc==1)]
zvort_ml = zvort_mv[:,(cc==1)]


ds = nc.Dataset(fp+f"cm1out_{fn:06d}.nc")
stime = ds.variables['time'][:].data[0]
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
z = ds.variables['z'][:].data
iz5 = np.where(z >= 5)[0][0]
iz = slice(0,iz5+1)
xlims = [-30,10] #[-55,25]
ix = slice(np.where(xh >= xlims[0])[0][0], np.where(xh >= xlims[1])[0][1])
ds.close()


# Load saved cross sections
dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/cross_sections/MV1_{mv_time}min/cross{stime/60:.0f}.pkl", 'rb')
cs = pickle.load(dbfile)
# winterp = cs['w']
thrpert = cs['thrpert']
prspert = cs['prspert']
dbfile.close()



w_lims = [-20,20]
thr_lims = [-16,16]
pp_lims = [-4,4]
# OW_lims = [-0.01,0]
zlims = [0,3.5]

w_levs = np.linspace(w_lims[0], w_lims[1], (w_lims[1]-w_lims[0])+1)
thr_levs = np.linspace(thr_lims[0], thr_lims[1], (thr_lims[1]-thr_lims[0])+1)
pp_levs = np.linspace(pp_lims[0], pp_lims[1], 2*(pp_lims[1]-pp_lims[0])+1)

qix = int(np.round(len(xh[(xh>=xlims[0]) & (xh<=xlims[1])])/60)*2) # scale wind vector spacing to number of grid pts
qiz = 3

ti = np.where(ptime == mv_time*60)[0][0]
# ti = np.where(ptime == 240*60)[0][0]
ti0 = np.where(ptime == 180*60)[0][0]

x_median = np.median(x_ml[ti0:ti+1,:], axis=1)
y_median = np.median(y_ml[ti0:ti+1,:], axis=1)
z_median = np.median(z_ml[ti0:ti+1,:], axis=1)
u_median = np.median(u_ml[ti0:ti+1,:], axis=1)
v_median = np.median(v_ml[ti0:ti+1,:], axis=1)
w_median = np.median(w_ml[ti0:ti+1,:], axis=1)
zvort_median = np.median(zvort_ml[ti0:ti+1,:], axis=1)

xlims = [-25, 0] # was [-30, 5]



figsave = False

# CROSS SECTION, 1 PANEL - shade w; contour thrpert; gray parcels
if False:
    fig,ax = plt.subplots(1,1,figsize=(9,4))
    plot_contourf(xh[ix], z[iz], winterp, 'w', ax, levels=w_levs, datalims=w_lims, xlims=xlims, ylims=zlims)
    ax.contour(xh[ix], z[iz], thrpert, levels=[-5,-1], colors='purple', linewidths=[2,1], linestyles='-')
    # ax.contour(xh[ix], z[iz], OW, levels=[-0.01,-0.0025,-0.001], colors='hotpink', linewidths=1, linestyles='-')
    ax.scatter(x_median, z_median, s=20, color='dimgray', marker='.')
    ax.scatter(x_median[::20], z_median[::20], s=50, color='k', marker='.')
    for i in range(len(x_median[::20])):
        if np.mod(i,2) == 0:
            plt.text(x_median[::20][i], z_median[::20][i]+0.1, f"{ptime[ti0:ti+1][::20][i]/60:.0f} min", fontsize=10, fontweight='bold')
        else:
            plt.text(x_median[::20][i], z_median[::20][i]-0.2, f"{ptime[ti0:ti+1][::20][i]/60:.0f} min", fontsize=10, fontweight='bold')
    ax.set_xlabel('x (km)')
    ax.set_ylabel('z (km)')
    ax.set_title(f"Mid-level parcels in the MV at {mv_time} min")
    
    if figsave:
        plt.savefig(f"/Users/morgan.schneider/Documents/merger/merger-125m/MV1_w_cross_{mv_time}min.png", dpi=300)
    

# CROSS SECTION, 1 PANEL - shade thrpert; contour p'; parcels colored by w
# ***PAPER FIG*** - runs both times separately, need to combine into two panels in powerpoint
if False: # this one is True
    # parcel_cm = 'ChaseSpectral'
    import matplotlib as mpl
    # parcel_cm = mpl.colors.LinearSegmentedColormap.from_list('parcel_cm',
    #                     np.vstack((plt.cm.Spectral_r(np.linspace(0,0.5,128)),
    #                     pyart.graph.cmweather.cm_colorblind.ChaseSpectral(np.linspace(0.5,1,128)))))
    parcel_cm = mpl.colors.LinearSegmentedColormap.from_list('parcel_cm',
                        np.vstack((pyart.graph.cmweather.cm_colorblind.ChaseSpectral(np.linspace(0.1,0.5,154)),
                        pyart.graph.cmweather.cm_colorblind.ChaseSpectral(np.linspace(0.5,1,102)))))
    wl = [-15,10]
    
    fig,ax = plt.subplots(1,1,figsize=(10,4), layout='constrained')
    plot_contourf(xh[ix], z[iz], thrpert, 'thrpert', ax, levels=thr_levs, datalims=thr_lims, cmap='PuOr_r', cbfs=13, xlims=xlims, ylims=zlims, extend='both')
    ax.contour(xh[ix], z[iz], prspert, levels=[-5,-4,-3,-2,-1,1,2,3,4,5], colors=['b','b','b','b','b','r','r','r','r','r'], linewidths=1, linestyles='-')
    # ax.contour(xh[ix], z[iz], OW, levels=[-0.01,-0.0025,-0.001], colors='hotpink', linewidths=1, linestyles='-')
    # ax.quiver(xh[ix][::qix], z[iz][::qiz], 0*winterp[::qiz,::qix], winterp[::qiz,::qix], scale=250, width=0.002, pivot='tail')
    # ax.scatter(x_median, z_median, s=20, color='dimgray', marker='.')
    p = ax.scatter(x_median, z_median, s=30, c=w_median, cmap=parcel_cm, edgecolor='none', linewidth=0.25, vmin=wl[0], vmax=wl[1])
    ax.scatter(x_median[::20], z_median[::20], s=50, c=w_median[::20], cmap=parcel_cm, edgecolor='k', linewidth=1.5, vmin=wl[0], vmax=wl[1])
    cb = plt.colorbar(p, ax=ax, extend='both') #, ticks=np.linspace(wl[0],wl[1],11))
    cb.set_label("Parcel $w$ (m s$^{-1}$)", fontsize=13)
    if mv_time == 210:
        for i in range(len(x_median[::20])):
            t = ptime[ti0:ti+1][::20][i]/60
            if (t < 200) & (np.mod(i,2) == 0):
                ax.text(x_median[::20][i]-0.5, z_median[::20][i]+0.1, f"{t:.0f}", fontsize=12, fontweight='bold')
            elif (t < 200) & (np.mod(i,2) != 0):
                ax.text(x_median[::20][i]-0.5, z_median[::20][i]-0.25, f"{t:.0f}", fontsize=12, fontweight='bold')
            elif t == 200:
                ax.text(x_median[::20][i]+0.4, z_median[::20][i], f"{t:.0f}", fontsize=12, fontweight='bold')
            elif t == 205:
                ax.text(x_median[::20][i]+0.4, z_median[::20][i], f"{t:.0f}", fontsize=12, fontweight='bold')
            elif t == 210:
                ax.text(x_median[::20][i]-1.0, z_median[::20][i]+0.1, f"{t:.0f}", fontsize=12, fontweight='bold')
    if mv_time == 220:
        ax.text(x_median[::20][0]+0.3, z_median[::20][0]-0.1, f"180", fontsize=12, fontweight='bold')
        # ax.text(x_median[::20][1]-0.1, z_median[::20][1]+0.1, f"185", fontsize=12, fontweight='bold')
        # ax.text(x_median[::20][2]-0.35, z_median[::20][2]-0.25, f"190", fontsize=12, fontweight='bold')
        ax.text(x_median[::20][2]-0.5, z_median[::20][2]-0.25, f"190", fontsize=12, fontweight='bold')
        # ax.text(x_median[::20][3]-1.0, z_median[::20][3]-0.2, f"195", fontsize=12, fontweight='bold')
        # ax.text(x_median[::20][3]-1.4, z_median[::20][3]-0.15, f"195", fontsize=12, fontweight='bold')
        # ax.text(x_median[::20][4]-0.1, z_median[::20][4]+0.1, f"200", fontsize=12, fontweight='bold')
        ax.text(x_median[::20][4]+0.05, z_median[::20][4]+0.06, f"200", fontsize=12, fontweight='bold')
        # ax.text(x_median[::20][5]-1.3, z_median[::20][5]-0.05, f"205", fontsize=12, fontweight='bold')
        ax.text(x_median[::20][5]-1.6, z_median[::20][5]-0.05, f"205", fontsize=12, fontweight='bold')
        ax.text(x_median[::20][6]+0.3, z_median[::20][6]-0.05, f"210", fontsize=12, fontweight='bold')
        # ax.text(x_median[::20][7]-1.3, z_median[::20][7]-0.05, f"215", fontsize=12, fontweight='bold')
        ax.text(x_median[::20][7]-1.6, z_median[::20][7]-0.05, f"215", fontsize=12, fontweight='bold')
        ax.text(x_median[::20][8], z_median[::20][8]+0.1, f"220", fontsize=12, fontweight='bold')
    if mv_time == 225:
        ax.text(x_median[::20][0]-0.5, z_median[::20][0]-0.25, f"180", fontsize=12, fontweight='bold')
        ax.text(x_median[::20][1]-0.5, z_median[::20][1]+0.1, f"185", fontsize=12, fontweight='bold')
        ax.text(x_median[::20][2]-0.25, z_median[::20][2]-0.25, f"190", fontsize=12, fontweight='bold')
        ax.text(x_median[::20][3]-1, z_median[::20][3]-0.25, f"195", fontsize=12, fontweight='bold')
        ax.text(x_median[::20][4]+0.3, z_median[::20][4]+0.05, f"200", fontsize=12, fontweight='bold')
        ax.text(x_median[::20][5]+0.4, z_median[::20][5], f"205", fontsize=12, fontweight='bold')
        ax.text(x_median[::20][6]-0.75, z_median[::20][6]+0.1, f"210", fontsize=12, fontweight='bold')
        ax.text(x_median[::20][7]-2.2, z_median[::20][7]-0.05, f"215", fontsize=12, fontweight='bold')
        ax.text(x_median[::20][8]+0.3, z_median[::20][8]+0.05, f"220", fontsize=12, fontweight='bold')
        ax.text(x_median[::20][9]+0.3, z_median[::20][9]+0.05, f"225", fontsize=12, fontweight='bold')
    ax.set_xlabel('x (km)', fontsize=14)
    ax.set_ylabel('z (km)', fontsize=14)
    # ax.set_title(f"Mid-level parcels in the MV at {mv_time} min", fontsize=14)
    
    if figsave:
        plt.savefig(f"/Users/morgan.schneider/Documents/merger/traj_thr_cross_{mv_time}min.png", dpi=300)




#% Plot plan views ***PAPER FIG***

dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/cross_sections/MV1_{mv_time}min/plan{stime/60:.0f}.pkl", 'rb')
pl = pickle.load(dbfile)
xlp = pl['xlims']
ylp = pl['ylims']
# w_plan = pl['w']
thr_plan = pl['thrpert']
# prs_plan = pl['prspert']
dbfile.close()

ix = slice(np.where(xh >= xlp[0])[0][0], np.where(xh >= xlp[1])[0][1])
iy = slice(np.where(yh >= ylp[0])[0][0], np.where(xh >= ylp[1])[0][1])


if mv_time == 210:
    xlp = [-40,10]
    # ylp = [-100,-50]
    ylp = [-95,-55]
elif mv_time == 220:
    xlp = [-40,10]
    # ylp = [-90,-40]
    ylp = [-85,-45]
elif mv_time == 225:
    xlp = [-40,10]
    ylp = [-85,-35]




# figsave = False

# PLAN VIEW, 1 PANEL - shade thrpert; parcels colored by w
# ***PAPER FIG*** - runs both times separately, need to combine into two panels in powerpoint
if False:
    import matplotlib as mpl
    # parcel_cm = mpl.colors.LinearSegmentedColormap.from_list('parcel_cm',
    #                     np.vstack((plt.cm.Spectral_r(np.linspace(0,0.5,128)),
    #                     pyart.graph.cmweather.cm_colorblind.ChaseSpectral(np.linspace(0.5,1,128)))))
    parcel_cm = mpl.colors.LinearSegmentedColormap.from_list('parcel_cm',
                        np.vstack((pyart.graph.cmweather.cm_colorblind.ChaseSpectral(np.linspace(0.1,0.5,154)),
                        pyart.graph.cmweather.cm_colorblind.ChaseSpectral(np.linspace(0.5,1,102)))))
    wl = [-15,10]
    
    fig,ax = plt.subplots(1,1,figsize=(9,6), subplot_kw=dict(box_aspect=1), layout='constrained')
    plot_contourf(xh[ix], yh[iy], thr_plan, 'thrpert', ax, levels=thr_levs, datalims=thr_lims, cmap='PuOr_r', cbfs=14, xlims=xlp, ylims=ylp, extend='both', cbar=True)
    # ax.contour(xh[ix], yh[iy], w_plan, levels=[-15,-10,-5,5,10,15], colors=['b','b','b','r','r','r'], linewidths=1, linestyles='-')
    p = ax.scatter(x_median, y_median, s=30, c=w_median, cmap=parcel_cm, edgecolor='none', linewidth=0.25, vmin=wl[0], vmax=wl[1])
    ax.scatter(x_median[::20], y_median[::20], s=50, c=w_median[::20], cmap=parcel_cm, edgecolor='k', linewidth=1.5, vmin=wl[0], vmax=wl[1])
    cb = plt.colorbar(p, ax=ax, extend='both') #, ticks=np.linspace(wl[0],wl[1],11))
    cb.set_label(label="Parcel $w$ (m s$^{-1}$)", fontsize=14)
    for i in range(len(x_median[::20])):
        if y_median[::20][i] >= ylp[0]:
            t = ptime[ti0:ti+1][::20][i]/60
            ax.text(x_median[::20][i]+1, y_median[::20][i], f"{t:.0f} min", fontsize=12, fontweight='bold')
    ax.set_xlabel('x (km)', fontsize=14)
    ax.set_ylabel('y (km)', fontsize=14)
    # ax.set_title(f"Mid-level parcels in the MV at {mv_time} min", fontsize=14)
    
    if figsave:
        plt.savefig(f"/Users/morgan.schneider/Documents/merger/traj_thr_plan_{mv_time}min.png", dpi=300)



#% Plot trajectories with vorticity and SR velocity vectors 
# ***PAPER FIG*** - runs both times separately, need to combine into two panels in powerpoint

from matplotlib.patches import FancyArrowPatch

if mv_time == 210:
    xlp = [-40,10]
    ylp = [-100,-50]
elif mv_time == 220:
    xlp = [-40,10]
    ylp = [-90,-40]

ix = slice(np.where(xh >= xlp[0])[0][0], np.where(xh >= xlp[1])[0][1])
iy = slice(np.where(yh >= ylp[0])[0][0], np.where(xh >= ylp[1])[0][1])

if True:
    ds = nc.Dataset(fp+f"cm1out_{fn+2:06d}.nc")
    dbz = ds.variables['dbz'][:].data[0,0,iy,ix]
    ds.close()
    
dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/hvort_traj_{mv_time:.0f}min.pkl", 'rb')
vort_traj = pickle.load(dbfile)
xvort_ml = vort_traj['xvort_ml']
yvort_ml = vort_traj['yvort_ml']
# hvort_ml = vort_traj['hvort_ml']
# vort_ml = vort_traj['vort_ml']
# vort_sw = vort_traj['vort_sw_ml']
# vort_cw = vort_traj['vort_cw_ml_signed']
dbfile.close()

xvort_median = np.median(xvort_ml[ti0:ti+1,:], axis=1)
yvort_median = np.median(yvort_ml[ti0:ti+1,:], axis=1)


dbfile = open('/Users/morgan.schneider/Documents/merger/merger-125m/storm_motion.pkl', 'rb')
sm = pickle.load(dbfile)
u_storm = sm['u_storm']
v_storm = sm['v_storm']
dbfile.close()

qit = 4

u_sr = u_median[::qit] - u_storm[:len(x_median[::qit])]
v_sr = v_median[::qit] - v_storm[:len(x_median[::qit])]

if mv_time == 210:
    xlp = [-25,0]
elif mv_time == 220:
    xlp = [-20,5]
elif mv_time == 225:
    xlp = [-16.5,8.5]


# figsave = False


def add_vectors(ax, x, y, dx, dy, z, lengthscale=10, *args, **kwargs):
    # SMALLER lengthscale for SHORTER arrows -- this is REVERSED from scale kwarg in quiver
    zinds = np.argsort(z)
    for j in range(len(x)):
        i = zinds[j]
        x1 = x[i]
        y1 = y[i]
        dx1 = dx[i]*lengthscale
        dy1 = dy[i]*lengthscale
        x2 = x1 + dx1
        y2 = y1 + dy1
        arrow = FancyArrowPatch((x1, y1), (x2, y2), *args, **kwargs)
        ax.add_patch(arrow)
        
    return arrow

from matplotlib.colors import LinearSegmentedColormap
cmap = LinearSegmentedColormap.from_list('cmap', pyart.graph.cmweather.cm_colorblind.ChaseSpectral(np.linspace(0.1,0.9,205)))



fig,ax = plt.subplots(1,1, figsize=(5,6), subplot_kw=dict(box_aspect=2), layout='constrained')

l, = ax.plot([xh[0], xh[1]], [yh[0], yh[1]], 'gray', linewidth=1)
ax.contour(xh[ix], yh[iy], dbz, levels=[30], colors='gray', linewidths=1, zorder=0)
p = ax.scatter(x_median, y_median, s=40, c=z_median, cmap=cmap, vmin=0, vmax=3)
cb = plt.colorbar(p, ax=ax, extend='max')
cb.set_label("Parcel height (km)", fontsize=13)
# p = ax.scatter(x_median, y_median, s=30, c=zvort_median, cmap='ChaseSpectral', vmin=-0.04, vmax=0.04)
# cb = plt.colorbar(p, ax=ax, label="Parcel \u03B6 (s$^{-1}$)", extend='both')
a_wind = add_vectors(ax, x_median[::qit], y_median[::qit], u_sr, v_sr, z_median[::qit],
            lengthscale=0.2, arrowstyle='simple', mutation_scale=5, ec='k', fc='k', lw=0.5)
a_vort = add_vectors(ax, x_median[::qit]-0.1, y_median[::qit], xvort_median[::qit], yvort_median[::qit], z_median[::qit],
            lengthscale=275, arrowstyle='simple', mutation_scale=8, ec='k', fc='lightgray', lw=0.75)
ax.set_xlabel('x (km)', fontsize=12)
ax.set_ylabel('y (km)', fontsize=12)
# ax.set_title(f"Parcels in the MV at {mv_time} min", fontsize=13)
ax.set_xlim(xlp) # [-25,0] for 210 min, [-17,8] for 225 min
ax.set_ylim(ylp)
ax.set_yticks(np.arange(ylp[0], ylp[1]+5, 5))
plt.legend(handles=[l,a_vort,a_wind], labels=['30 dBZ',"\u03c9$_H$","SR wind"], loc=1, fontsize=10)

if figsave:
    plt.savefig(f"/Users/morgan.schneider/Documents/merger/traj_vort2d_parcelheight_{mv_time:.0f}min.png", dpi=300)



#%% Same figures but with storm-relative trajectories ***PAPER FIGS HERE***

from scipy.interpolate import interp1d

fp = '/Volumes/Promise_Pegasus_70TB/merger/merger-125m/'
ip = '/Users/morgan.schneider/Documents/merger/merger-125m/'

mv_time = 210

if mv_time == 210:
    fn = 41
elif mv_time == 220:
    fn = 51
elif mv_time == 225:
    fn = 56




# Load parcel time
ds = nc.Dataset(fp+'cm1out_pdata.nc')
ptime = ds.variables['time'][:].data
ds.close()

ti = np.where(ptime == mv_time*60)[0][0]
ti0 = np.where(ptime == 180*60)[0][0]

# Load MV parcel data
dbfile = open('/Users/morgan.schneider/Documents/merger/traj_MV1.pkl', 'rb')
traj = pickle.load(dbfile)
pids_mv = traj[f"{mv_time}min"]['pids']
x_mv = traj[f"{mv_time}min"]['x']
y_mv = traj[f"{mv_time}min"]['y']
z_mv = traj[f"{mv_time}min"]['z']
u_mv = traj[f"{mv_time}min"]['u']
v_mv = traj[f"{mv_time}min"]['v']
w_mv = traj[f"{mv_time}min"]['w']
dbfile.close()

# Load parcel source region clusters
dbfile = open(ip+f"traj_clusters/traj_clusters_{mv_time}min_v2.pkl", 'rb')
c = pickle.load(dbfile)
cc = c['mv1']
dbfile.close()


# Mid-level source (0 -> LL inflow, 1 -> ML inflow, 2 -> supercell outflow, 3 -> QLCS outflow)
x_ml = x_mv[:,(cc==1)]
y_ml = y_mv[:,(cc==1)]
z_ml = z_mv[:,(cc==1)]
u_ml = u_mv[:,(cc==1)]
v_ml = v_mv[:,(cc==1)]
w_ml = w_mv[:,(cc==1)]


# Load saved storm motion estimates
dbfile = open('/Users/morgan.schneider/Documents/merger/merger-125m/storm_motion.pkl', 'rb')
sm = pickle.load(dbfile)
u_storm = sm['u_storm']
v_storm = sm['v_storm']
dbfile.close()

stimes = np.linspace(10800, 14400, 61)


# Interpolate storm motion from model output timesto parcel output times (60 s to 15 s)
fu = interp1d(stimes, u_storm)
u_storm_interp = fu(ptime)
fv = interp1d(stimes, v_storm)
v_storm_interp = fv(ptime)


# Get storm-relative trajectories starting from MV
x_sr = np.zeros(shape=x_ml[:ti+1,:].shape, dtype=float)
y_sr = np.zeros(shape=y_ml[:ti+1,:].shape, dtype=float)
x_sr[ti,:] = x_ml[ti,:]
y_sr[ti,:] = y_ml[ti,:]
inds = np.linspace(ti-1, 0, ti)
for i in inds:
    i = int(i)
    dt = ptime[i+1] - ptime[i]
    delta_x = np.sum(u_storm_interp[i:ti] * dt)
    delta_y = np.sum(v_storm_interp[i:ti] * dt)
    # delta_x = u_storm_interp[i] * dt
    # delta_y = v_storm_interp[i] * dt
    x_sr[i,:] = x_ml[i] + delta_x
    y_sr[i,:] = y_ml[i] + delta_y


# Median parcel trajectory
x_median = np.median(x_sr[ti0:ti+1,:], axis=1)/1000
y_median = np.median(y_sr[ti0:ti+1,:], axis=1)/1000
z_median = np.median(z_ml[ti0:ti+1,:], axis=1)/1000
u_median = np.median(u_ml[ti0:ti+1,:], axis=1)
v_median = np.median(v_ml[ti0:ti+1,:], axis=1)
w_median = np.median(w_ml[ti0:ti+1,:], axis=1)

# Load saved parcel trajectory vorticity (interpolated from model fields)
dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/hvort_traj_{mv_time}min.pkl", 'rb')
vort_traj = pickle.load(dbfile)
xvort_ml = vort_traj['xvort_ml']
yvort_ml = vort_traj['yvort_ml']
dbfile.close()

xvort_median = np.median(xvort_ml[ti0:ti+1,:], axis=1)
yvort_median = np.median(yvort_ml[ti0:ti+1,:], axis=1)

# only plot every qit SR wind vectors (not actually used here oops)
qit = 4
u_sr = u_median[::qit] - u_storm[:len(x_median[::qit])]
v_sr = v_median[::qit] - v_storm[:len(x_median[::qit])]



ds = nc.Dataset(fp+f"cm1out_{fn:06d}.nc")
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
z = ds.variables['z'][:].data
ds.close()




figsave = False



### TRAJECTORIES WITH VORTICITY VECTORS
if False:
    from matplotlib.patches import FancyArrowPatch
    
    def add_vectors(ax, x, y, dx, dy, z, lengthscale=10, *args, **kwargs):
        # SMALLER lengthscale for SHORTER arrows -- this is REVERSED from scale kwarg in quiver
        zinds = np.argsort(z)
        for j in range(len(x)):
            i = zinds[j]
            x1 = x[i]
            y1 = y[i]
            dx1 = dx[i]*lengthscale
            dy1 = dy[i]*lengthscale
            x2 = x1 + dx1
            y2 = y1 + dy1
            arrow = FancyArrowPatch((x1, y1), (x2, y2), *args, **kwargs)
            ax.add_patch(arrow)
            
        return arrow
    
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list('cmap', pyart.graph.cmweather.cm_colorblind.ChaseSpectral(np.linspace(0.1,0.9,205)))
    
    if mv_time == 210:
        ix = slice(np.where(xh >= -30)[0][0], np.where(xh >= 20)[0][1])
        iy = slice(np.where(yh >= -85)[0][0], np.where(xh >= -45)[0][1])
    elif mv_time == 220:
        ix = slice(np.where(xh >= -25)[0][0], np.where(xh >= 25)[0][1])
        iy = slice(np.where(yh >= -80)[0][0], np.where(xh >= -40)[0][1])
    
    ds = nc.Dataset(fp+f"cm1out_{fn+2:06d}.nc")
    dbz = ds.variables['dbz'][:].data[0,0,iy,ix]
    ds.close()
    
    if mv_time == 210:
        xlp = [-25,15]
    elif mv_time == 220:
        xlp = [-20,20]
    elif mv_time == 225:
        xlp = [-16.5,8.5]
    
    
    fig,ax = plt.subplots(1,1, figsize=(7,5), subplot_kw=dict(box_aspect=1), layout='constrained')
    
    l, = ax.plot([xh[0], xh[1]], [yh[0], yh[1]], 'gray', linewidth=1)
    ax.contour(xh[ix], yh[iy], dbz, levels=[30], colors='gray', linewidths=1, zorder=0)
    p = ax.scatter(x_median, y_median, s=80, c=z_median, cmap=cmap, vmin=0, vmax=3)
    cb = plt.colorbar(p, ax=ax, extend='max')
    cb.set_label("Parcel height (km)", fontsize=13)
    # p = ax.scatter(x_median, y_median, s=30, c=zvort_median, cmap='ChaseSpectral', vmin=-0.04, vmax=0.04)
    # cb = plt.colorbar(p, ax=ax, label="Parcel \u03B6 (s$^{-1}$)", extend='both')
    # a_wind = add_vectors(ax, x_median[::qit], y_median[::qit], u_sr, v_sr, z_median[::qit],
    #             lengthscale=0.2, arrowstyle='simple', mutation_scale=5, ec='k', fc='k', lw=0.5)
    a_vort = add_vectors(ax, x_median[::qit], y_median[::qit]+0.5, xvort_median[::qit], yvort_median[::qit], z_median[::qit],
                lengthscale=275, arrowstyle='simple', mutation_scale=9, ec='k', fc='lightgray', lw=0.75)
    ax.set_xlabel('x (km)', fontsize=12)
    ax.set_ylabel('y (km)', fontsize=12)
    # ax.set_title(f"Parcels in the MV at {mv_time} min", fontsize=13)
    ax.set_xlim(xlp)
    ax.set_ylim(ylp)
    ax.set_yticks(np.arange(ylp[0], ylp[1]+5, 5))
    # plt.legend(handles=[l,a_vort,a_wind], labels=['30 dBZ',"\u03c9$_H$","SR wind"], loc=1, fontsize=10)
    plt.legend(handles=[l,a_vort], labels=['30 dBZ',"\u03c9$_H$"], loc=1, fontsize=10)
    
    if figsave:
        plt.savefig(f"/Users/morgan.schneider/Documents/merger/traj_vort2d_parcelheight_{mv_time:.0f}min_SR.png", dpi=300)




import matplotlib as mpl
parcel_cm = mpl.colors.LinearSegmentedColormap.from_list('parcel_cm',
                    np.vstack((pyart.graph.cmweather.cm_colorblind.ChaseSpectral(np.linspace(0.1,0.5,154)),
                    pyart.graph.cmweather.cm_colorblind.ChaseSpectral(np.linspace(0.5,1,102)))))


### TRAJECTORY PLAN VIEWS
if False:
    if mv_time == 210:
        xl = [-30,20]
        # yl = [-100,-50]
        yl = [-90,-40]
    elif mv_time == 220:
        xl = [-30,20]
        # yl = [-90,-40]
        yl = [-85,-35]
    
    ix = slice(np.where(xh >= xl[0])[0][0], np.where(xh >= xl[1])[0][1])
    iy = slice(np.where(yh >= yl[0])[0][0], np.where(xh >= yl[1])[0][1])
    iz1 = slice(0, np.where(z >= 1)[0][1])
    
    ds = nc.Dataset(fp+f"cm1out_{fn:06d}.nc")
    thr_plan = ds.variables['thrpert'][:].data[0,iy,ix]
    # w_plan = np.max(ds.variables['winterp'][:].data[0,iz1,iy,ix], axis=0)
    ds.close()
    
    wl = [-15,10]
    
    fig,ax = plt.subplots(1,1,figsize=(9,6), subplot_kw=dict(box_aspect=1), layout='constrained')
    plot_contourf(xh[ix], yh[iy], thr_plan, 'thrpert', ax, levels=np.linspace(-16,16,17), datalims=[-16,16], cmap='PuOr_r', cbfs=14, xlims=xl, ylims=yl, extend='both', cbar=True)
    # ax.contour(xh[ix], yh[iy], w_plan, levels=[5,10,15], colors=['gray','k','k'], linewidths=1, linestyles='-')
    p = ax.scatter(x_median, y_median, s=80, c=w_median, cmap=parcel_cm, edgecolor='none', linewidth=0.25, vmin=wl[0], vmax=wl[1])
    ax.scatter(x_median[::20], y_median[::20], s=80, c=w_median[::20], cmap=parcel_cm, edgecolor='k', linewidth=1.5, vmin=wl[0], vmax=wl[1])
    cb = plt.colorbar(p, ax=ax, extend='both') #, ticks=np.linspace(wl[0],wl[1],11))
    cb.set_label(label="Parcel $w$ (m s$^{-1}$)", fontsize=14)
    for i in range(len(x_median[::20])):
        t = ptime[ti0:ti+1][::20][i]/60
        if mv_time == 210:
            ax.text(x_median[::20][i]-2, y_median[::20][i]+1.5, f"{t:.0f}", fontsize=14, fontweight='bold')
        elif mv_time == 220:
            ax.text(x_median[::20][i]+0.75, y_median[::20][i]+0.5, f"{t:.0f}", fontsize=14, fontweight='bold')
    ax.set_xlabel('x (km)', fontsize=14)
    ax.set_ylabel('y (km)', fontsize=14)
    # ax.set_title(f"Mid-level parcels in the MV at {mv_time} min", fontsize=14)
    
    if figsave:
        plt.savefig(f"/Users/morgan.schneider/Documents/merger/traj_thr_plan_{mv_time}min_SR.png", dpi=300)





# TRAJECTORY CROSS SECTIONS
if True:
    xl = [-30,15] #[-55,25]
    zl = [0,3.5]
    iz = slice(0, np.where(z >= zl[1])[0][1])
    ix = slice(np.where(xh >= xl[0])[0][0], np.where(xh >= xl[1])[0][1])
    iy = (np.abs(yh - np.median(y_ml[ti-9,:]/1000))).argmin()
    
    ds = nc.Dataset(fp+'cm1out_000013.nc')
    prs0 = ds.variables['prs0'][:].data[0,iz,iy,ix]
    ds.close()
    
    ds = nc.Dataset(fp+f"cm1out_{fn:06d}.nc")
    thrpert = ds.variables['thrpert'][:].data[iz,iy,ix]
    prspert = (ds.variables['prs'][:].data[0,iz,iy,ix] - prs0)/100
    # winterp = ds.variables['winterp'][:].data[0,iz,iy,ix]
    ds.close()
    
    wl = [-15,10]
    
    fig,ax = plt.subplots(1,1,figsize=(10,4), layout='constrained')
    plot_contourf(xh[ix], z[iz], thrpert, 'thrpert', ax, levels=np.linspace(-16,16,17), datalims=[-16,16], cmap='PuOr_r', cbfs=13, xlims=xl, ylims=zl, extend='both')
    # ax.contour(xh[ix], z[iz], prspert, levels=[-5,-4,-3,-2,-1,1,2,3,4,5], colors=['b','b','b','b','b','r','r','r','r','r'], linewidths=1, linestyles='-')
    ax.contour(xh[ix], z[iz], prspert, levels=[-3,-2,-1,1,2,3], colors=['b','b','b','r','r','r'], linewidths=0.75, linestyles='-')
    # ax.contour(xh[ix], z[iz], OW, levels=[-0.01,-0.0025,-0.001], colors='hotpink', linewidths=1, linestyles='-')
    # ax.quiver(xh[ix][::qix], z[iz][::qiz], 0*winterp[::qiz,::qix], winterp[::qiz,::qix], scale=250, width=0.002, pivot='tail')
    # ax.scatter(x_median, z_median, s=20, color='dimgray', marker='.')
    p = ax.scatter(x_median, z_median, s=60, c=w_median, cmap=parcel_cm, edgecolor='none', linewidth=0.25, vmin=wl[0], vmax=wl[1])
    ax.scatter(x_median[::20], z_median[::20], s=60, c=w_median[::20], cmap=parcel_cm, edgecolor='k', linewidth=1.5, vmin=wl[0], vmax=wl[1])
    cb = plt.colorbar(p, ax=ax, extend='both') #, ticks=np.linspace(wl[0],wl[1],11))
    cb.set_label("Parcel $w$ (m s$^{-1}$)", fontsize=13)
    for i in range(len(x_median[::20])):
        t = ptime[ti0:ti+1][::20][i]/60
        if mv_time == 210:
            ax.text(x_median[::20][i]+0.3, z_median[::20][i]+0.05, f"{t:.0f}", fontsize=12, fontweight='bold')
        elif mv_time == 220:
            if (i==1) | (i==3):
                x_offset = -1; z_offset = -0.25
            elif (i == 6):
                x_offset = -0.2; z_offset = 0.1
            elif i >= 7:
                x_offset = -2.25; z_offset = 0.1
            else:
                x_offset = -1; z_offset = 0.1
            ax.text(x_median[::20][i]+x_offset, z_median[::20][i]+z_offset, f"{t:.0f}", fontsize=12, fontweight='bold')
    ax.set_xlabel('x (km)', fontsize=14)
    ax.set_ylabel('z (km)', fontsize=14)
    # ax.set_title(f"Mid-level parcels in the MV at {mv_time} min", fontsize=14)
    
    if figsave:
        plt.savefig(f"/Users/morgan.schneider/Documents/merger/traj_thr_cross_{mv_time}min_SR.png", dpi=300)
    
    
    
    
    
    


    
#%% Animate w, thrpert, prspert cross sections

from matplotlib.animation import FuncAnimation

ip = '/Users/morgan.schneider/Documents/merger/merger-125m/'
mv_time = 220

# Read parcel data
ds = nc.Dataset(fp+'cm1out_pdata.nc')
ptime = ds.variables['time'][:].data
ds.close()

w_lims = [-20,20]
thr_lims = [-16,16]
pp_lims = [-4,4]
# OW_lims = [-0.01,0]
zlims = [0,3.5]

ds = nc.Dataset('/Volumes/Promise_Pegasus_70TB/merger/merger-125m/cm1out_000014.nc')
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
z = ds.variables['z'][:].data
ds.close()

iz5 = np.where(z >= 5)[0][0]
iz = slice(0,iz5+1)

w_levs = np.linspace(w_lims[0], w_lims[1], (w_lims[1]-w_lims[0])+1)
thr_levs = np.linspace(thr_lims[0], thr_lims[1], (thr_lims[1]-thr_lims[0])+1)
pp_levs = np.linspace(pp_lims[0], pp_lims[1], 2*(pp_lims[1]-pp_lims[0])+1)


xlims = [-30,10] #[-55,25]

ix = slice(np.where(xh >= xlims[0])[0][0], np.where(xh >= xlims[1])[0][1])

qix = int(np.round(len(xh[(xh>=xlims[0]) & (xh<=xlims[1])])/60)*2) # scale wind vector spacing to number of grid pts
qiz = 3

xl = [-25,5]


if False:
    dbfile = open(ip+f"cross_sections/MV1_{mv_time}min/cross180.pkl", 'rb')
    cs = pickle.load(dbfile)
    stime = cs['time']
    w_cross = cs['w']
    prs_cross = cs['prspert']
    thr_cross = cs['thrpert']
    # OW_cross = cs['OW']
    # zvort_cross = cs['zvort']
    # B_cross = cs['B']
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
        # OW_cross = cs['OW']
        # zvort_cross = cs['zvort']
        # B_cross = cs['B']
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
    
    anim = FuncAnimation(fig, animate_w_cross, frames=73-12, interval=100, repeat=False, blit=False)
    if figsave:
        anim.save(ip+f"MV1_w_cross_{mv_time}min.gif", dpi=300)
    plt.show()



if True:
    dbfile = open(ip+f"cross_sections/MV1_{mv_time}min/cross205.pkl", 'rb')
    cs = pickle.load(dbfile)
    stime = cs['time']
    # w_cross = cs['w']
    prs_cross = cs['prspert']
    thr_cross = cs['thrpert']
    # OW_cross = cs['OW']
    # zvort_cross = cs['zvort']
    # B_cross = cs['B']
    dbfile.close()
    
    ti = np.where(ptime == stime)[0][0]
    
    fig,ax = plt.subplots(1,1,figsize=(9,4))
    plot_contourf(xh[ix], z[iz], thr_cross, 'thrpert', ax, levels=thr_levs, datalims=thr_lims, xlims=xl, ylims=zlims, cmap='PuOr_r')
    ax.contour(xh[ix], z[iz], prs_cross, levels=[1,2,3,4,5], colors='r', linewidths=1)
    ax.contour(xh[ix], z[iz], prs_cross, levels=[-5,-4,-3,-2,-1], colors='b', linewidths=1, linestyles='-')
    # ax.contour(xh[ix], z[iz], OW, levels=[-0.01,-0.0025,-0.001], colors='hotpink', linewidths=1, linestyles='-')
    # ax.quiver(xh[ix][::qix], z[iz][::qiz], 0*w_cross[::qiz,::qix], w_cross[::qiz,::qix], scale=250, width=0.002, pivot='tail')
    ax.scatter(x_ml[ti,:], z_ml[ti,:], s=20, color='k', marker='.')
    ax.set_xlabel('x (km)')
    ax.set_ylabel('z (km)')
    ax.set_title(f"Parcels in the mesocyclone at {mv_time} min \n t={stime/60:.0f} min")
    plt.tight_layout()
    
    def animate_thr_cross(i):
        global ax
        ax.clear()
        
        dbfile = open(ip+f"cross_sections/MV1_{mv_time}min/cross{i+205}.pkl", 'rb')
        cs = pickle.load(dbfile)
        stime = cs['time']
        # w_cross = cs['w']
        prs_cross = cs['prspert']
        thr_cross = cs['thrpert']
        # OW_cross = cs['OW']
        # zvort_cross = cs['zvort']
        # B_cross = cs['B']
        dbfile.close()
        ti = np.where(ptime == stime)[0][0]
        
        plot_contourf(xh[ix], z[iz], thr_cross, 'thrpert', ax, levels=thr_levs, datalims=thr_lims, xlims=xl, ylims=zlims, cmap='PuOr_r', cbar=False)
        ax.contour(xh[ix], z[iz], prs_cross, levels=[1,2,3,4,5], colors='r', linewidths=1)
        ax.contour(xh[ix], z[iz], prs_cross, levels=[-5,-4,-3,-2,-1], colors='b', linewidths=1, linestyles='-')
        # ax.contour(xh[ix], z[iz], OW, levels=[-0.01,-0.0025,-0.001], colors='hotpink', linewidths=1, linestyles='-')
        # ax.quiver(xh[ix][::qix], z[iz][::qiz], 0*w_cross[::qiz,::qix], w_cross[::qiz,::qix], scale=250, width=0.002, pivot='tail')
        ax.scatter(x_ml[ti,:], z_ml[ti,:], s=20, color='k', marker='.')
        ax.set_xlabel('x (km)')
        ax.set_ylabel('z (km)')
        ax.set_title(f"Parcels in the mesocyclone at {mv_time} min \n t={stime/60:.0f} min")
        plt.tight_layout()
    
    figsave = True
    
    anim = FuncAnimation(fig, animate_thr_cross, frames=16, interval=500, repeat=False, blit=False)
    if figsave:
        anim.save(ip+f"MV1_thr_cross_{mv_time}min_slow.gif", dpi=300)
    plt.show()


if False:
    dbfile = open(ip+f"cross_sections/MV1_{mv_time}min/cross180.pkl", 'rb')
    cs = pickle.load(dbfile)
    stime = cs['time']
    # w_cross = cs['w']
    prs_cross = cs['prspert']
    # thr_cross = cs['thrpert']
    OW_cross = cs['OW']
    # zvort_cross = cs['zvort']
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
        # w_cross = cs['w']
        prs_cross = cs['prspert']
        # thr_cross = cs['thrpert']
        OW_cross = cs['OW']
        # zvort_cross = cs['zvort']
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
    
    anim = FuncAnimation(fig, animate_prs_cross, frames=73-12, interval=100, repeat=False, blit=False)
    if figsave:
        anim.save(ip+f"MV1_prs_cross_{mv_time}min_extended.gif", dpi=300)
    plt.show()


#%% Animate plan views

from matplotlib.animation import FuncAnimation

ip = '/Users/morgan.schneider/Documents/merger/merger-125m/'
mv_time = 210

w_lims = [-20,20]
thr_lims = [-16,16]
pp_lims = [-4,4]
# OW_lims = [-0.01,0]
zlims = [0,3.5]

ds = nc.Dataset('/Volumes/Promise_Pegasus_70TB/merger/merger-125m/cm1out_000014.nc')
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
z = ds.variables['z'][:].data
ds.close()

iz5 = np.where(z >= 5)[0][0]
iz = slice(0,iz5+1)

w_levs = np.linspace(w_lims[0], w_lims[1], (w_lims[1]-w_lims[0])+1)
thr_levs = np.linspace(thr_lims[0], thr_lims[1], (thr_lims[1]-thr_lims[0])+1)
pp_levs = np.linspace(pp_lims[0], pp_lims[1], 2*(pp_lims[1]-pp_lims[0])+1)



if False:
    dbfile = open(ip+f"cross_sections/MV1_{mv_time}min/plan180.pkl", 'rb')
    cs = pickle.load(dbfile)
    stime = cs['time']
    xlims = cs['xlims']
    ylims = cs['ylims']
    w_plan = cs['w'] # 500m
    # prs_plan = cs['prspert'] # mean 0-500m
    # thr_plan = cs['thrpert'] # sfc
    # OW_plan = cs['OW'] # min 0-500m
    # zvort_plan = cs['zvort'] # max 0-500m
    # B_plan = cs['B'] # mean 0-500m
    dbfile.close()
    
    ix = slice(np.where(xh >= xlims[0])[0][0], np.where(xh >= xlims[1])[0][1])
    iy = slice(np.where(yh >= ylims[0])[0][0], np.where(yh >= ylims[1])[0][1])
    
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    plot_contourf(xh[ix], yh[iy], w_plan, 'w', ax, levels=w_levs, datalims=w_lims, xlims=xlims, ylims=ylims)
    # ax.contour(xh[ix], yh[iy], thr_plan, levels=[-5,-1], colors='purple', linewidths=1, linestyles=['-','--'])
    # ax.contour(xh[ix], yh[iy], OW_plan, levels=[-0.01,-0.0025,-0.001], colors='hotpink', linewidths=1, linestyles='-')
    # ax.contour(xh[ix], yh[iy], prs_plan, levels=[-4,-3,-2,-1,1,2,3,4], colors='purple', linewidths=1)
    ax.plot([-30,10], [np.median(y_ml[0,:]),np.median(y_ml[0,:])], 'k', linewidth=0.75)
    # ax.scatter(x_ll[0,:], y_ll[0,:], s=20, color='gray', marker='.')
    # ax.scatter(x_sc[0,:], y_sc[0,:], s=20, color='gray', marker='.')
    # ax.scatter(x_qc[0,:], y_qc[0,:], s=20, color='gray', marker='.')
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
        # prs_plan = cs['prspert']
        # thr_plan = cs['thrpert']
        # OW_plan = cs['OW']
        # zvort_plan = cs['zvort']
        # B_plan = cs['B']
        dbfile.close()
        
        ti = np.where(ptime == stime)[0][0]
        ix = slice(np.where(xh >= xlims[0])[0][0], np.where(xh >= xlims[1])[0][1])
        iy = slice(np.where(yh >= ylims[0])[0][0], np.where(yh >= ylims[1])[0][1])
        
        plot_contourf(xh[ix], yh[iy], w_plan, 'w', ax, levels=w_levs, datalims=w_lims, xlims=xlims, ylims=ylims, cbar=False)
        # ax.contour(xh[ix], yh[iy], thr_plan, levels=[-5,-1], colors='purple', linewidths=1, linestyles=['-','--'])
        # ax.contour(xh[ix], yh[iy], OW_plan, levels=[-0.01,-0.0025,-0.001], colors='hotpink', linewidths=1, linestyles='-')
        # ax.contour(xh[ix], yh[iy], prs_plan, levels=[-4,-3,-2,-1,1,2,3,4], colors='purple', linewidths=1)
        ax.plot([-30,10], [np.median(y_ml[ti,:]),np.median(y_ml[ti,:])], 'k', linewidth=0.75)
        # ax.scatter(x_ll[ti,:], y_ll[ti,:], s=20, color='gray', marker='.')
        # ax.scatter(x_sc[ti,:], y_sc[ti,:], s=20, color='gray', marker='.')
        # ax.scatter(x_qc[ti,:], y_qc[ti,:], s=20, color='gray', marker='.')
        ax.scatter(x_ml[ti,:], y_ml[ti,:], s=50, color='k', marker='.')
        ax.set_xlabel('x (km)')
        ax.set_ylabel('y (km)')
        ax.set_title(f"Mid-level parcels in the MV at {mv_time} min \n t={stime/60:.0f} min")
        plt.tight_layout()
    
    figsave = True
    
    anim = FuncAnimation(fig, animate_w_plan, frames=73-12, interval=100, repeat=False, blit=False)
    if figsave:
        anim.save(ip+f"MV1_w_plan_{mv_time}min_extended.gif", dpi=300)
    plt.show()
    
    
if True:
    dbfile = open(ip+f"cross_sections/MV1_{mv_time}min/plan180.pkl", 'rb')
    cs = pickle.load(dbfile)
    stime = cs['time']
    xlims = cs['xlims']
    ylims = cs['ylims']
    # w_plan = cs['w'] # 500m
    # prs_plan = cs['prspert'] # mean 0-500m
    thr_plan = cs['thrpert'] # sfc
    # OW_plan = cs['OW'] # min 0-500m
    # zvort_plan = cs['zvort'] # max 0-500m
    # B_plan = cs['B'] # mean 0-500m
    dbfile.close()
    
    ti = np.where(ptime == stime)[0][0]
    
    ix = slice(np.where(xh >= xlims[0])[0][0], np.where(xh >= xlims[1])[0][1])
    iy = slice(np.where(yh >= ylims[0])[0][0], np.where(yh >= ylims[1])[0][1])
    
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    plot_contourf(xh[ix], yh[iy], thr_plan, 'thrpert', ax, levels=thr_levs, datalims=thr_lims, xlims=xlims, ylims=ylims, cmap='PuOr_r')
    # ax.contour(xh[ix], yh[iy], B_plan, levels=[1,2,3,4,5], colors='r', linewidths=1)
    # ax.contour(xh[ix], yh[iy], OW_plan, levels=[-0.01,-0.0025,-0.001], colors='hotpink', linewidths=1, linestyles='-')
    # ax.contour(xh[ix], yh[iy], w_plan, levels=[-20,-15,-10,-5,5,10,15,20], colors='purple', linewidths=1)
    ax.plot(xl, [np.median(y_ml[0,:]),np.median(y_ml[0,:])], 'k', linewidth=0.75)
    # ax.scatter(x_ll[0,:], y_ll[0,:], s=20, color='gray', marker='.')
    # ax.scatter(x_sc[0,:], y_sc[0,:], s=20, color='gray', marker='.')
    # ax.scatter(x_qc[0,:], y_qc[0,:], s=20, color='gray', marker='.')
    ax.scatter(x_ml[ti,:], y_ml[ti,:], s=10, edgecolor='k', facecolor='w', marker='o')
    ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')
    ax.set_title(f"Parcels in the mesocyclone at {mv_time} min \n t={stime/60:.0f} min")
    plt.tight_layout()
    
    def animate_thr_plan(i):
        global ax
        ax.clear()
        
        dbfile = open(ip+f"cross_sections/MV1_{mv_time}min/plan{i+180}.pkl", 'rb')
        cs = pickle.load(dbfile)
        stime = cs['time']
        xlims = cs['xlims']
        ylims = cs['ylims']
        # w_plan = cs['w']
        # prs_plan = cs['prspert']
        thr_plan = cs['thrpert']
        # OW_plan = cs['OW']
        # zvort_plan = cs['zvort']
        # B_plan = cs['B']
        dbfile.close()
        
        ti = np.where(ptime == stime)[0][0]
        ix = slice(np.where(xh >= xlims[0])[0][0], np.where(xh >= xlims[1])[0][1])
        iy = slice(np.where(yh >= ylims[0])[0][0], np.where(yh >= ylims[1])[0][1])
        
        plot_contourf(xh[ix], yh[iy], thr_plan, 'thrpert', ax, levels=thr_levs, datalims=thr_lims, xlims=xlims, ylims=ylims, cmap='PuOr_r', cbar=False)
        # ax.contour(xh[ix], yh[iy], B_plan, levels=[1,2,3,4,5], colors='r', linewidths=1)
        # ax.contour(xh[ix], yh[iy], OW_plan, levels=[-0.01,-0.0025,-0.001], colors='hotpink', linewidths=1, linestyles='-')
        # ax.contour(xh[ix], yh[iy], w_plan, levels=[-20,-15,-10,-5,5,10,15,20], colors='purple', linewidths=1)
        ax.plot(xl, [np.median(y_ml[ti,:]),np.median(y_ml[ti,:])], 'k', linewidth=0.75)
        # ax.scatter(x_ll[ti,:], y_ll[ti,:], s=20, color='gray', marker='.')
        # ax.scatter(x_sc[ti,:], y_sc[ti,:], s=20, color='gray', marker='.')
        # ax.scatter(x_qc[ti,:], y_qc[ti,:], s=20, color='gray', marker='.')
        ax.scatter(x_ml[ti,:], y_ml[ti,:], s=10, edgecolor='k', facecolor='w', marker='o')
        ax.set_xlabel('x (km)')
        ax.set_ylabel('y (km)')
        ax.set_title(f"Parcels in the mesocyclone at {mv_time} min \n t={stime/60:.0f} min")
        plt.tight_layout()
    
    figsave = True
    
    anim = FuncAnimation(fig, animate_thr_plan, frames=46, interval=100, repeat=False, blit=False)
    if figsave:
        anim.save(ip+f"MV1_thr_plan_{mv_time}min.gif", dpi=300)
    plt.show()
    
    
    
    
if False:
    dbfile = open(ip+f"cross_sections/MV1_{mv_time}min/plan180.pkl", 'rb')
    cs = pickle.load(dbfile)
    stime = cs['time']
    xlims = cs['xlims']
    ylims = cs['ylims']
    w_plan = cs['w'] # 500m
    prs_plan = cs['prspert'] # mean 0-500m
    # thr_plan = cs['thrpert'] # sfc
    # OW_plan = cs['OW'] # min 0-500m
    # zvort_plan = cs['zvort'] # max 0-500m
    # B_plan = cs['B'] # mean 0-500m
    dbfile.close()
    
    ix = slice(np.where(xh >= xlims[0])[0][0], np.where(xh >= xlims[1])[0][1])
    iy = slice(np.where(yh >= ylims[0])[0][0], np.where(yh >= ylims[1])[0][1])
    
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    plot_contourf(xh[ix], yh[iy], prs_plan, 'prspert', ax, levels=pp_levs, datalims=pp_lims, xlims=xlims, ylims=ylims)
    # ax.contour(xh[ix], yh[iy], B_plan, levels=[1,2,3,4,5], colors='r', linewidths=1)
    # ax.contour(xh[ix], yh[iy], OW_plan, levels=[-0.01,-0.0025,-0.001], colors='hotpink', linewidths=1, linestyles='-')
    ax.contour(xh[ix], yh[iy], w_plan, levels=[-20,-15,-10,-5,5,10,15,20], colors='purple', linewidths=1)
    ax.plot([-30,10], [np.median(y_ml[0,:]),np.median(y_ml[0,:])], 'k', linewidth=0.75)
    # ax.scatter(x_ll[0,:], y_ll[0,:], s=20, color='gray', marker='.')
    # ax.scatter(x_sc[0,:], y_sc[0,:], s=20, color='gray', marker='.')
    # ax.scatter(x_qc[0,:], y_qc[0,:], s=20, color='gray', marker='.')
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
        # thr_plan = cs['thrpert']
        # OW_plan = cs['OW']
        # zvort_plan = cs['zvort']
        # B_plan = cs['B']
        dbfile.close()
        
        ti = np.where(ptime == stime)[0][0]
        ix = slice(np.where(xh >= xlims[0])[0][0], np.where(xh >= xlims[1])[0][1])
        iy = slice(np.where(yh >= ylims[0])[0][0], np.where(yh >= ylims[1])[0][1])
        
        plot_contourf(xh[ix], yh[iy], prs_plan, 'prspert', ax, levels=pp_levs, datalims=pp_lims, xlims=xlims, ylims=ylims, cbar=False)
        # ax.contour(xh[ix], yh[iy], B_plan, levels=[1,2,3,4,5], colors='r', linewidths=1)
        # ax.contour(xh[ix], yh[iy], OW_plan, levels=[-0.01,-0.0025,-0.001], colors='hotpink', linewidths=1, linestyles='-')
        ax.contour(xh[ix], yh[iy], w_plan, levels=[-20,-15,-10,-5,5,10,15,20], colors='purple', linewidths=1)
        ax.plot([-30,10], [np.median(y_ml[ti,:]),np.median(y_ml[ti,:])], 'k', linewidth=0.75)
        # ax.scatter(x_ll[ti,:], y_ll[ti,:], s=20, color='gray', marker='.')
        # ax.scatter(x_sc[ti,:], y_sc[ti,:], s=20, color='gray', marker='.')
        # ax.scatter(x_qc[ti,:], y_qc[ti,:], s=20, color='gray', marker='.')
        ax.scatter(x_ml[ti,:], y_ml[ti,:], s=50, color='k', marker='.')
        ax.set_xlabel('x (km)')
        ax.set_ylabel('y (km)')
        ax.set_title(f"Mid-level parcels in the MV at {mv_time} min \n t={stime/60:.0f} min")
        plt.tight_layout()
    
    figsave = True
    
    anim = FuncAnimation(fig, animate_prs_plan, frames=73-12, interval=100, repeat=False, blit=False)
    if figsave:
        anim.save(ip+f"MV1_prs_plan_{mv_time}min_extended.gif", dpi=300)
    plt.show()

