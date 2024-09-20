#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 16:32:28 2024

@author: morgan.schneider
"""

from CM1utils import *
import wrf

#%% Calculate time heights

fp = '/Volumes/Promise_Pegasus_70TB/merger/supercell-125m/'
ip = '/Users/morgan.schneider/Documents/merger/supercell-125m/'

ds = nc.Dataset(fp+'cm1out_000014.nc')
xh = ds.variables['xh'][:].data
xf = ds.variables['xf'][:].data
yh = ds.variables['yh'][:].data
yf = ds.variables['yf'][:].data
z = ds.variables['z'][:].data
zf = ds.variables['zf'][:].data
ds.close()

dbfile = open(ip+'boxes_s2.pkl', 'rb')
box = pickle.load(dbfile)
x1 = box['x1_pp']
x2 = box['x2_pp']
y1 = box['y1_pp']
y2 = box['y2_pp']
dbfile.close()

iz5 = np.where(z >= 5.0)[0][0]
iz = slice(0,iz5)
iz1 = np.where(z >= 1.0)[0][0]
iz200 = np.where(z >= 0.2)[0][0]

# times = np.zeros(shape=(60,), dtype=float)
# w_max = np.zeros(shape=(60,len(z[iz])), dtype=float)
# zvort_max = np.zeros(shape=(60,len(z[iz])), dtype=float)
# OW_min = np.zeros(shape=(60,len(z[iz])), dtype=float)
# wspd_max = np.zeros(shape=(60,len(z[iz])), dtype=float)
# uh25_max = np.zeros(shape=(60,), dtype=float)
# uh02_max = np.zeros(shape=(60,), dtype=float)
pp_min = np.zeros(shape=(60,len(z[iz])), dtype=float)
vppga_1km_max = np.zeros(shape=(60,), dtype=float)
# p_dn_min = np.zeros(shape=(60,len(z[iz])), dtype=float)
# vppga_dn_max = np.zeros(shape=(60,), dtype=float)
# vppga_b_max = np.zeros(shape=(60,), dtype=float)
# vppga_dl_max = np.zeros(shape=(60,), dtype=float)

i = 0
for fn in np.arange(14,74):
    print(f"cm1out_{fn:06d}")
    ix1 = np.where(xh >= x1[i])[0][0]
    ix2 = np.where(xh >= x2[i])[0][0]
    iy1 = np.where(yh >= y1[i])[0][0]
    iy2 = np.where(yh >= y2[i])[0][0]
    
    ix = slice(ix1,ix2)
    iy = slice(iy1,iy2)
    
    # ds = nc.Dataset(fp+f"cm1out_{fn:06d}.nc")
    # time = ds.variables['time'][:].data[0]/60
    # w = ds.variables['winterp'][:].data[0,iz,iy,ix]
    # zvort = ds.variables['zvort'][:].data[0,iz,iy,ix]
    # u = ds.variables['uinterp'][:].data[0,iz,iy,ix]
    # v = ds.variables['vinterp'][:].data[0,iz,iy,ix]
    # wspd = np.sqrt(u**2 + v**2)
    # S_N = np.gradient(u, xh[ix]*1000, axis=2) - np.gradient(v, yh[iy]*1000, axis=1)
    # S_S = np.gradient(v, xh[ix]*1000, axis=2) + np.gradient(u, yh[iy]*1000, axis=1)
    # OW = S_N**2 + S_S**2 - zvort**2
    # del u,v,S_N,S_S
    # ds.close()
    
    # times[i] = time
    # w_max[i,:] = np.max(w, axis=(1,2))
    # zvort_max[i,:] = np.max(zvort, axis=(1,2))
    # OW_min[i,:] = np.min(OW, axis=(1,2))
    # wspd_max[i,:] = np.max(wspd, axis=(1,2))
    # del w,zvort,OW,wspd
    
    
    # ds = nc.Dataset(fp+f"cm1out_{fn:06d}.nc")
    # u = ds.variables['uinterp'][:].data[0,:,iy,ix]
    # v = ds.variables['vinterp'][:].data[0,:,iy,ix]
    # w = ds.variables['w'][:].data[0,:,iy,ix]
    # ds.close()
    
    # zf3d = np.tile(np.expand_dims(zf,(1,2)), (1,len(yh[iy]),len(xh[ix])))
    # mapfct = np.ones(shape=(len(yh[iy]),len(xh[ix])), dtype=float)
    # uh_25km = wrf.to_np(wrf.udhel(zf3d*1000, mapfct, u, v, w, 125, 125, bottom=2000., top=5000.))
    # uh_02km = wrf.to_np(wrf.udhel(zf3d*1000, mapfct, u, v, w, 125, 125, bottom=0., top=2000.))
    # del u,v,w,zf3d,mapfct
    # uh_25km = np.asarray(uh_25km.data)
    # uh_02km = np.asarray(uh_02km.data)
    # uh25_max[i] = np.max(uh_25km)
    # uh02_max[i] = np.max(uh_02km)
    # del uh_25km,uh_02km
    
    
    ds = nc.Dataset(fp+'base/cm1out_000013.nc')
    prs0 = ds.variables['prs0'][:].data[0,iz,iy,ix]
    ds.close()
    ds = nc.Dataset(fp+f"cm1out_{fn:06d}.nc")
    prspert = ds.variables['prs'][:].data[0,iz,iy,ix] - prs0
    del prs0
    ds.close()
    
    dp_1km = prspert[iz1,:,:] - prspert[0,:,:]
    dz = (z[iz1] - z[0])*1000
    vppga_1km = -1/1.1 * dp_1km / dz
    pp_min[i,:] = np.min(prspert, axis=(1,2))
    vppga_1km_max[i] = np.max(vppga_1km)
    del prspert,dp_1km,vppga_1km
    
    
    # ds = nc.Dataset(fp+f"pp/dyn_{fn:06d}.nc")
    # dz = (z[iz1] - z[iz200])*1000
    # p_dn = ds.variables['p_dn'][:].data[iz,iy,ix]
    # p_b = ds.variables['p_b'][:].data[iz,iy,ix]
    # p_dl = ds.variables['p_dl'][:].data[iz,iy,ix]
    # ds.close()
    
    # dp_dn = p_dn[iz1,:,:] - p_dn[iz200,:,:]
    # pgfz_dn = -1/1.1 * dp_dn / dz
    # dp_b = p_b[iz1,:,:] - p_b[iz200,:,:]
    # pgfz_b = -1/1.1 * dp_b / dz
    # dp_dl = p_dl[iz1,:,:] - p_dl[iz200,:,:]
    # pgfz_dl = -1/1.1 * dp_dl / dz
    # del dp_dn,dp_b,dp_dl
    
    # p_dn_min[i,:] = np.min(p_dn, axis=(1,2))
    # vppga_dn_max[i] = np.max(pgfz_dn)
    # vppga_b_max[i] = np.max(pgfz_b)
    # vppga_dl_max[i] = np.max(pgfz_dl)
    # del pgfz_dn,pgfz_b,pgfz_dl,p_dn,p_b,p_dl
    
    i = i + 1


#%% Save or update pickle file

new_pkl = False
update_pkl = False
fn = 'timeheight_s2.pkl'

if new_pkl:
    timeheight = {'times':times, 'z':z[iz], 'w_max':w_max.transpose(), 'zvort_max':zvort_max.transpose(), 
                   'OW_min':OW_min.transpose(), 'wspd_max':wspd_max.transpose(),
                   'uh25_max':uh25_max, 'uh02_max':uh02_max}
    dbfile = open(ip+fn, 'wb')
    pickle.dump(timeheight, dbfile)
    dbfile.close()
    
if update_pkl:
    dbfile = open(ip+fn, 'rb')
    timeheight = pickle.load(dbfile)
    dbfile.close()
    
    new_vars = {'pp_min':pp_min.transpose(), 'vppga_1km_max':vppga_1km_max}
    timeheight.update(new_vars)
    dbfile = open(ip+fn, 'wb')
    pickle.dump(timeheight, dbfile)
    dbfile.close()


#%% Load time heights from saved pickle files

fp = '/Volumes/Promise_Pegasus_70TB/merger/supercell-125m/'
ip = '/Users/morgan.schneider/Documents/merger/supercell-125m/'
sim = 'SUPERCELL (S2)'
fn = 's2'

dbfile = open(ip+f"timeheight_{fn}.pkl", 'rb')
tmp = pickle.load(dbfile)
times = tmp['times']
zz = tmp['z']
w_max = tmp['w_max']
zvort_max = tmp['zvort_max']
OW_min = tmp['OW_min']
wspd_max = tmp['wspd_max']
uh25_max = tmp['uh25_max']
uh02_max = tmp['uh02_max']
pp_min = tmp['pp_min']
vppga_1km_max = tmp['vppga_1km_max']

# p_dn_min = tmp['p_dn_min']
# vppga_dn_max = tmp['vppga_dn_max']
# vppga_b_max = tmp['vppga_b_max']
# vppga_dl_max = tmp['vppga_dl_max']
dbfile.close()


#%% Plot timeheights from pickle files

ztop = 4
figsave = False


# w, zvort, OW, wspd
if True:
    fig,ax = plt.subplots(1,1,figsize=(12,5))
    plot_cfill(times, zz, w_max, 'w', ax, datalims=[0,20], cmap='Reds')
    ax.contour(times, zz, w_max, levels=[10,15,20], colors='k', linewidths=1)
    ax.axvline(210, color='w', linewidth=1.5, linestyle='--')
    ax.set_xlabel('Time (min)')
    ax.set_ylabel('Height (km)')
    ax.set_ylim([0,ztop])
    ax.set_title(f"{sim} - Maximum w")
    if figsave:
        plt.savefig(ip+f"imgs_timeheights/timeheight_{fn}_w.png", dpi=300)
    
    fig,ax = plt.subplots(1,1,figsize=(12,5))
    plot_cfill(times, zz, zvort_max, 'zvort', ax, datalims=[0,0.1], cmap='Reds')
    ax.contour(times, zz, zvort_max, levels=[0.1], colors='k', linewidths=1)
    ax.axvline(210, color='w', linewidth=1.5, linestyle='--')
    ax.set_xlabel('Time (min)')
    ax.set_ylabel('Height (km)')
    ax.set_ylim([0,ztop])
    ax.set_title(f"{sim} - Maximum \u03B6")
    if figsave:
        plt.savefig(ip+f"imgs_timeheights/timeheight_{fn}_zvort.png", dpi=300)
    
    fig,ax = plt.subplots(1,1,figsize=(12,5))
    plot_cfill(times, zz, OW_min, 'OW', ax, datalims=[-0.01,0], cmap='Blues_r')
    ax.contour(times, zz, OW_min, levels=[-0.01,-0.005,-0.001], colors='k', linewidths=1, linestyles='-')
    ax.axvline(210, color='k', linewidth=1.5, linestyle='--')
    ax.set_xlabel('Time (min)')
    ax.set_ylabel('Height (km)')
    ax.set_ylim([0,ztop])
    ax.set_title(f"{sim} - Minimum OW")
    if figsave:
        plt.savefig(ip+f"imgs_timeheights/timeheight_{fn}_OW.png", dpi=300)
    
    fig,ax = plt.subplots(1,1,figsize=(12,5))
    plot_cfill(times, zz, wspd_max, 'wspd', ax, datalims=[25,50], cmap='pyart_HomeyerRainbow')
    # ax.contour(times, zz, wspd_max, levels=[25], colors='k', linewidths=1, linestyles='-')
    ax.axvline(210, color='k', linewidth=1.5, linestyle='--')
    ax.set_xlabel('Time (min)')
    ax.set_ylabel('Height (km)')
    ax.set_ylim([0,3])
    ax.set_title(f"{sim} - Maximum wind speed")
    if figsave:
        plt.savefig(ip+f"imgs_timeheights/timeheight_{fn}_wspd.png", dpi=300)
    
    fig,ax = plt.subplots(1,1,figsize=(12,5))
    plot_cfill(times, zz, w_max, 'w', ax, datalims=[0,20], cmap='Reds')
    ax.contour(times, zz, OW_min, levels=[-0.01,-0.005,-0.001], colors='k', linestyles=['-','-','--'], linewidths=[2,1,1])
    ax.axvline(210, color='w', linewidth=1.5, linestyle='--')
    ax.set_xlabel('Time (min)')
    ax.set_ylabel('Height (km)')
    ax.set_ylim([0,ztop])
    ax.set_title(f"{sim} - Maximum w, minimum OW")
    if figsave:
        plt.savefig(ip+f"imgs_timeheights/timeheight_{fn}_w+OW.png", dpi=300)
    
    fig,ax = plt.subplots(1,1,figsize=(12,5))
    plot_cfill(times, zz, OW_min, 'OW', ax, datalims=[-0.01,0], cmap='Blues_r')
    ax.contour(times, zz, w_max, levels=[10,15,20], colors=['k','k','k'], linestyles=['--','-','-'], linewidths=[1,1,2])
    ax.axvline(210, color='k', linewidth=1.5, linestyle='--')
    ax.set_xlabel('Time (min)')
    ax.set_ylabel('Height (km)')
    ax.set_ylim([0,ztop])
    ax.set_title(f"{sim} - Minimum OW, maximum w")
    if figsave:
        plt.savefig(ip+f"imgs_timeheights/timeheight_{fn}_OW+w.png", dpi=300)

# Pressure perturbation
if True:
    fig,ax = plt.subplots(1,1,figsize=(12,5))
    plot_cfill(times, zz, pp_min/100, 'prspert', ax, datalims=[-4,0], cmap='Blues_r')
    # ax.contour(times, zz, pp_min/100, levels=[-5,-3,-1], colors='k', linewidths=1, linestyles='-')
    ax.axvline(210, color='k', linewidth=1.5, linestyle='--')
    ax.set_xlabel('Time (min)')
    ax.set_ylabel('Height (km)')
    ax.set_ylim([0,ztop])
    ax.set_title(f"{sim} - Minimum p'")
    if figsave:
        plt.savefig(ip+f"imgs_timeheights/timeheight_{fn}_prspert.png", dpi=300)
    
    fig,ax = plt.subplots(1,1,figsize=(12,5))
    ax.plot(times, vppga_1km_max, linewidth=1.5)
    ax.axvline(210, color='k', linewidth=1.5, linestyle='--')
    ax.set_xlabel('Time (min)')
    ax.set_xlim([180,240])
    ax.set_ylim([0,0.5])
    ax.set_title(f"{sim} - Maximum 0-1 km VPPGA")
    if figsave:
        plt.savefig(ip+f"imgs_timeheights/timeseries_{fn}_vppga.png", dpi=300)

# Updraft helicity
if False:
    fig,ax = plt.subplots(1,1,figsize=(12,5))
    ax.plot(times, uh25_max, linewidth=1.5)
    ax.plot(times, uh02_max, linewidth=1.5)
    ax.axvline(210, color='k', linewidth=1.5, linestyle='--')
    ax.set_xlabel('Time (min)')
    ax.set_xlim([180,240])
    ax.set_ylim([0,2000])
    ax.set_title(f"{sim} - Maximum updraft helicity")
    plt.legend(['2-5 km', '0-2 km'])
    if figsave:
        plt.savefig(ip+f"imgs_timeheights/timeseries_{fn}_UH.png", dpi=300)

# all simulations combined UH
if False:
    dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/timeheight_{fn}.pkl", 'rb')
    tmp = pickle.load(dbfile)
    times = tmp['times']
    uh25_max_m = tmp['uh25_max']
    uh02_max_m = tmp['uh02_max']
    dbfile.close()
    dbfile = open(f"/Users/morgan.schneider/Documents/merger/qlcs-125m/timeheight_{fn}.pkl", 'rb')
    tmp = pickle.load(dbfile)
    times = tmp['times']
    uh25_max_q = tmp['uh25_max']
    uh02_max_q = tmp['uh02_max']
    dbfile.close()
    if 's' in fn:
        dbfile = open(f"/Users/morgan.schneider/Documents/merger/supercell-125m/timeheight_{fn}.pkl", 'rb')
        tmp = pickle.load(dbfile)
        uh25_max_s = tmp['uh25_max']
        uh02_max_s = tmp['uh02_max']
        dbfile.close()
    
    
    fig,((ax1),(ax2)) = plt.subplots(2,1,figsize=(12,9))
    l1, = ax1.plot(times, uh25_max_q) # k
    l2, = ax1.plot(times, uh25_max_m) # r
    if 's' in fn:
        l3, = ax1.plot(times, uh25_max_s) # b
        h = [l1,l2,l3]
        l = ['QLCS','Merger','Supercell']
    else:
        h = [l1,l2]
        l = ['QLCS','Merger']
    ax1.axvline(210, color='k', linewidth=1.5, linestyle='--')
    ax1.set_xlim([180,240])
    ax1.set_ylim([0,2000])
    ax1.set_title(f"Maximum 2-5 km UH")
    ax1.legend(handles=h, labels=l)
    ax2.plot(times, uh02_max_q) # k
    ax2.plot(times, uh02_max_m) # r
    if 's' in fn:
        ax2.plot(times, uh02_max_s) # b
    ax2.axvline(210, color='k', linewidth=1.5, linestyle='--')
    ax2.set_xlabel('Time (min)')
    ax2.set_xlim([180,240])
    ax2.set_ylim([0,1500])
    ax2.set_title(f"Maximum 0-2 km UH")
    ax2.set_yticks([0,300,600,900,1200,1500])
    if figsave:
        plt.savefig(ip+f"imgs_timeheights/timeseries_{fn}_UH_all.png", dpi=300)


# # Pressure perturbation components
# if False:
#     fig,ax = plt.subplots(1,1,figsize=(12,5))
#     plot_cfill(times, zz, p_dn_min/100, 'prspert', ax, datalims=[-6,0], cmap='pyart_HomeyerRainbow')
#     ax.contour(times, zz, w_max, levels=[15,20], colors='white', linewidths=1.5)
#     # ax.contour(times, zz, zvort_max, levels=[0.09,0.1], colors='k', linewidths=1.5, linestyles='-')
#     ax.contour(times, zz, OW_min, levels=[-0.01,-0.001], colors='k', linewidths=1.5, linestyles='-')
#     ax.set_xlabel('Time (min)')
#     ax.set_ylabel('Height (km)')
#     ax.set_ylim([0,ztop])
#     ax.set_title(f"{sim} - Minimum $p_{{DN}}'$")
#     if figsave:
#         plt.savefig(ip+f"imgs_timeheights/timeheight_{fn}_pdn_min.png", dpi=300)
    
    
#     fig,ax = plt.subplots(1,1,figsize=(12,5))
#     ax.plot(times, vppga_dn_max, 'k', linewidth=1.5)
#     ax.plot(times, vppga_dl_max, 'purple', linewidth=1.5)
#     ax.plot(times, vppga_b_max, 'b', linewidth=1.5)
#     ax.set_xlabel('Time (min)')
#     ax.set_xlim([180,240])
#     ax.set_ylim([0,0.5])
#     ax.set_title(f"{sim} - Maximum 200-1000 m VPPGA")
#     plt.legend(['Dynamic (nonlinear)', 'Dynamic (linear)', 'Buoyant'])
#     if figsave:
#         plt.savefig(ip+f"imgs_timeheights/timeseries_{fn}_vppga_max.png", dpi=300)


#%% Multi panel timeheight figures

figsave = False
ztop = 4


### SUPERCELL 1 ###

dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/timeheight_s1.pkl", 'rb')
tmp = pickle.load(dbfile)
times = tmp['times']
zz = tmp['z']
w_max_m = tmp['w_max']
OW_min_m = tmp['OW_min']
wspd_max_m = tmp['wspd_max']
uh25_max_m = tmp['uh25_max']
uh02_max_m = tmp['uh02_max']
# pp_min_m = tmp['pp_min']
# vppga_1km_max_m = tmp['vppga_1km_max']
dbfile.close()

dbfile = open(f"/Users/morgan.schneider/Documents/merger/qlcs-125m/timeheight_s1.pkl", 'rb')
tmp = pickle.load(dbfile)
w_max_q = tmp['w_max']
OW_min_q = tmp['OW_min']
wspd_max_q = tmp['wspd_max']
uh25_max_q = tmp['uh25_max']
uh02_max_q = tmp['uh02_max']
# pp_min_q = tmp['pp_min']
# vppga_1km_max_q = tmp['vppga_1km_max']
dbfile.close()

dbfile = open(f"/Users/morgan.schneider/Documents/merger/supercell-125m/timeheight_s1.pkl", 'rb')
tmp = pickle.load(dbfile)
w_max_s = tmp['w_max']
OW_min_s = tmp['OW_min']
wspd_max_s = tmp['wspd_max']
uh25_max_s = tmp['uh25_max']
uh02_max_s = tmp['uh02_max']
# pp_min_s = tmp['pp_min']
# vppga_1km_max_s = tmp['vppga_1km_max']
dbfile.close()



fig,((ax1),(ax2),(ax3)) = plt.subplots(3,1,figsize=(12,15))

plot_cfill(times, zz, OW_min_m, 'OW', ax1, datalims=[-0.01,0], cmap='Blues_r')
ax1.contour(times, zz, w_max_m, levels=[10,15,20], colors=['k','k','k'], linestyles=['--','-','-'], linewidths=[1,1,2])
ax1.axvline(210, color='k', linewidth=1.5, linestyle='--')
# ax1.set_xlabel('Time (min)')
ax1.set_ylabel('Height (km)')
ax1.set_ylim([0,ztop])
ax1.set_title('MERGER')

plot_cfill(times, zz, OW_min_q, 'OW', ax2, datalims=[-0.01,0], cmap='Blues_r')
ax2.contour(times, zz, w_max_q, levels=[10,15,20], colors=['k','k','k'], linestyles=['--','-','-'], linewidths=[1,1,2])
ax2.axvline(210, color='k', linewidth=1.5, linestyle='--')
# ax2.set_xlabel('Time (min)')
ax2.set_ylabel('Height (km)')
ax2.set_ylim([0,ztop])
ax2.set_title('QLCS')

plot_cfill(times, zz, OW_min_s, 'OW', ax3, datalims=[-0.01,0], cmap='Blues_r')
ax3.contour(times, zz, w_max_s, levels=[10,15,20], colors=['k','k','k'], linestyles=['--','-','-'], linewidths=[1,1,2])
ax3.axvline(210, color='k', linewidth=1.5, linestyle='--')
ax3.set_xlabel('Time (min)')
ax3.set_ylabel('Height (km)')
ax3.set_ylim([0,ztop])
ax3.set_title('SUPERCELL')
if figsave:
    plt.savefig(f"/Users/morgan.schneider/Documents/merger/timeheights_s1_OW+w.png", dpi=300)



fig,((ax1),(ax2)) = plt.subplots(2,1,figsize=(12,9))

ax1.plot(times, uh02_max_q, linewidth=1.5)
ax1.plot(times, uh02_max_m, linewidth=1.5)
ax1.plot(times, uh02_max_s, linewidth=1.5)
ax1.axvline(210, color='k', linewidth=1.5, linestyle='--')
ax1.set_xlabel('Time (min)')
ax1.set_xlim([180,240])
ax1.set_ylim([0,1500])
ax1.set_title(f"Maximum 0-2 km updraft helicity")
plt.legend(['No merger', 'Merger', 'Supercell'])

ax2.plot(times, uh25_max_q, linewidth=1.5)
ax2.plot(times, uh25_max_m, linewidth=1.5)
ax2.plot(times, uh25_max_s, linewidth=1.5)
ax2.axvline(210, color='k', linewidth=1.5, linestyle='--')
ax2.set_xlabel('Time (min)')
ax2.set_xlim([180,240])
ax2.set_ylim([0,2000])
ax2.set_title(f"Maximum 2-5 km updraft helicity")
if figsave:
    plt.savefig(f"/Users/morgan.schneider/Documents/merger/timeseries_s1_UH.png", dpi=300)




### SUPERCELL 2 ###

dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/timeheight_s2.pkl", 'rb')
tmp = pickle.load(dbfile)
w_max_m2 = tmp['w_max']
OW_min_m2 = tmp['OW_min']
wspd_max_m2 = tmp['wspd_max']
uh25_max_m2 = tmp['uh25_max']
uh02_max_m2 = tmp['uh02_max']
# pp_min_m2 = tmp['pp_min']
# vppga_1km_max_m2 = tmp['vppga_1km_max']
dbfile.close()

dbfile = open(f"/Users/morgan.schneider/Documents/merger/supercell-125m/timeheight_s2.pkl", 'rb')
tmp = pickle.load(dbfile)
w_max_s2 = tmp['w_max']
OW_min_s2 = tmp['OW_min']
wspd_max_s2 = tmp['wspd_max']
uh25_max_s2 = tmp['uh25_max']
uh02_max_s2 = tmp['uh02_max']
# pp_min_s2 = tmp['pp_min']
# vppga_1km_max_s2 = tmp['vppga_1km_max']
dbfile.close()


fig,((ax1),(ax2),(ax3)) = plt.subplots(3,1,figsize=(12,15))

plot_cfill(times, zz, OW_min_m2, 'OW', ax1, datalims=[-0.01,0], cmap='Blues_r')
ax1.contour(times, zz, w_max_m2, levels=[10,15,20], colors=['k','k','k'], linestyles=['--','-','-'], linewidths=[1,1,2])
ax1.axvline(210, color='k', linewidth=1.5, linestyle='--')
# ax1.set_xlabel('Time (min)')
ax1.set_ylabel('Height (km)')
ax1.set_ylim([0,ztop])
ax1.set_title('MERGER')

plot_cfill(times, zz, OW_min_s2, 'OW', ax2, datalims=[-0.01,0], cmap='Blues_r')
ax2.contour(times, zz, w_max_s2, levels=[10,15,20], colors=['k','k','k'], linestyles=['--','-','-'], linewidths=[1,1,2])
ax2.axvline(210, color='k', linewidth=1.5, linestyle='--')
# ax2.set_xlabel('Time (min)')
ax2.set_ylabel('Height (km)')
ax2.set_ylim([0,ztop])
ax2.set_title('SUPERCELL')

if figsave:
    plt.savefig(f"/Users/morgan.schneider/Documents/merger/timeheights_s2_OW+w.png", dpi=300)



fig,((ax1),(ax2)) = plt.subplots(2,1,figsize=(12,9))

ax1.plot(times, uh02_max_m2, linewidth=1.5)
ax1.plot(times, uh02_max_s2, linewidth=1.5)
ax1.axvline(210, color='k', linewidth=1.5, linestyle='--')
ax1.set_xlabel('Time (min)')
ax1.set_xlim([180,240])
ax1.set_ylim([0,1500])
ax1.set_title(f"Maximum 0-2 km updraft helicity")
plt.legend(['Merger', 'Supercell'])

ax2.plot(times, uh25_max_m2, linewidth=1.5)
ax2.plot(times, uh25_max_s2, linewidth=1.5)
ax2.axvline(210, color='k', linewidth=1.5, linestyle='--')
ax2.set_xlabel('Time (min)')
ax2.set_xlim([180,240])
ax2.set_ylim([0,2000])
ax2.set_title(f"Maximum 2-5 km updraft helicity")
if figsave:
    plt.savefig(f"/Users/morgan.schneider/Documents/merger/timeseries_s2_UH.png", dpi=300)



