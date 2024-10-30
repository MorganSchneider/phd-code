#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 16:39:33 2024

@author: morgan.schneider
"""

from CM1utils import *
import wrf
from scipy.ndimage import gaussian_filter


#%% Calculate swaths

fp = '/Volumes/Promise_Pegasus_70TB/merger/qlcs-125m/'
ip = '/Users/morgan.schneider/Documents/merger/qlcs-125m/'


ds = nc.Dataset(fp+'cm1out_000014.nc')
xh = ds.variables['xh'][:].data
xf = ds.variables['xf'][:].data
yh = ds.variables['yh'][:].data
yf = ds.variables['yf'][:].data
z = ds.variables['z'][:].data
zf = ds.variables['zf'][:].data
ds.close()

ix1 = np.where(xh >= -60)[0][0]
ix2 = np.where(xh >= 30)[0][0]
iy1 = np.where(yh >= -110)[0][0]
iy2 = np.where(yh >= -20)[0][0]
iz1 = np.where(z >= 1.0)[0][0]
iz15 = np.where(z >= 1.5)[0][0]
iz2 = np.where(z >= 2.0)[0][0]
iz3 = np.where(z >= 3.0)[0][0]

ix = slice(ix1,ix2)
iy = slice(iy1,iy2)

# w_2km_swath = np.zeros(shape=(len(yh[iy]),len(xh[ix])), dtype=float) # max 2 km w
# zvort_sfc_swath = np.zeros(shape=(len(yh[iy]),len(xh[ix])), dtype=float) # max sfc zvort
# OW_sfc_swath = np.zeros(shape=(len(yh[iy]),len(xh[ix])), dtype=float) # min sfc OW
# wspd_sfc_swath = np.zeros(shape=(len(yh[iy]),len(xh[ix])), dtype=float) # max sfc wind speed
# conv_sfc_swath = np.zeros(shape=(len(yh[iy]),len(xh[ix])), dtype=float) # min sfc convergence
# w_1km_swath = np.zeros(shape=(len(yh[iy]),len(xh[ix])), dtype=float) # max 1 km w
# zvort_1km_swath = np.zeros(shape=(len(yh[iy]),len(xh[ix])), dtype=float) # max 1 km zvort
# OW_1km_swath = np.zeros(shape=(len(yh[iy]),len(xh[ix])), dtype=float) # min 1 km OW
# OW_2km_swath = np.zeros(shape=(len(yh[iy]),len(xh[ix])), dtype=float) # min 2 km OW
OW_3km_swath = np.zeros(shape=(len(yh[iy]),len(xh[ix])), dtype=float) # min 3 km OW
# wspd_1km_swath = np.zeros(shape=(len(yh[iy]),len(xh[ix])), dtype=float) # max 1 km wind speed
# UH_25_swath = np.zeros(shape=(len(yh[iy]),len(xh[ix])), dtype=float) # max 2-5 km UH
# UH_02_swath = np.zeros(shape=(len(yh[iy]),len(xh[ix])), dtype=float) # max 0-2 km UH
# pp_sfc_swath = np.zeros(shape=(len(yh[iy]),len(xh[ix])), dtype=float) # min sfc total prspert
# pp_1km_swath = np.zeros(shape=(len(yh[iy]),len(xh[ix])), dtype=float) # min 1 km total prspert
# vppga_max_swath = np.zeros(shape=(len(yh[iy]),len(xh[ix])), dtype=float) # max upward total VPPGA
# vppga_min_swath = np.zeros(shape=(len(yh[iy]),len(xh[ix])), dtype=float) # min downward total VPPGA


# p_dn_sfc_swath = np.zeros(shape=(len(yh[iy]),len(xh[ix])), dtype=float) # min sfc p' (dn)
# p_dn_1km_swath = np.zeros(shape=(len(yh[iy]),len(xh[ix])), dtype=float) # min 0-1 km mean p' (dn)
# vppga_dn_max_swath = np.zeros(shape=(len(yh[iy]),len(xh[ix])), dtype=float) # max upward VPPGA (dn)
# vppga_dn_min_swath = np.zeros(shape=(len(yh[iy]),len(xh[ix])), dtype=float) # min downward VPPGA (dn)
# p_b_sfc_swath = np.zeros(shape=(len(yh[iy]),len(xh[ix])), dtype=float) # min sfc p' (b)
# p_b_1km_swath = np.zeros(shape=(len(yh[iy]),len(xh[ix])), dtype=float) # min 0-1 km mean p' (b)
# vppga_b_swath = np.zeros(shape=(len(yh[iy]),len(xh[ix])), dtype=float) # min downward VPPGA (b)

for fn in np.arange(14,74):
    print(f"{fn:06d}")
    
    ### Surface TLV variables ###
    # ds = nc.Dataset(fp+f"cm1out_{fn:06d}.nc")
    # w_2km = ds.variables['winterp'][:].data[0,iz2,iy,ix]
    # zvort_sfc = ds.variables['zvort'][:].data[0,0,iy,ix]
    # u_sfc = ds.variables['uinterp'][:].data[0,0,iy,ix]
    # v_sfc = ds.variables['vinterp'][:].data[0,0,iy,ix]
    # S_N = np.gradient(u_sfc, xh[ix]*1000, axis=1) - np.gradient(v_sfc, yh[iy]*1000, axis=0)
    # S_S = np.gradient(v_sfc, xh[ix]*1000, axis=1) + np.gradient(u_sfc, yh[iy]*1000, axis=0)
    # OW_sfc = S_N**2 + S_S**2 - zvort_sfc**2
    # wspd_sfc = (u_sfc**2 + v_sfc**2)**0.5
    # conv_sfc = np.gradient(u_sfc, xh[ix]*1000, axis=1) + np.gradient(v_sfc, yh[iy]*1000, axis=0)
    # del u_sfc,v_sfc,S_N,S_S
    # 
    # ds.close()
    
    # w_2km_swath = np.maximum(w_2km_swath, w_2km)
    # zvort_sfc_swath = np.maximum(zvort_sfc_swath, zvort_sfc)
    # OW_sfc_swath = np.minimum(OW_sfc_swath, OW_sfc)
    # wspd_sfc_swath = np.maximum(wspd_sfc_swath, wspd_sfc)
    # conv_sfc_swath = np.minimum(conv_sfc_swath, conv_sfc)
    # del w_2km,zvort_sfc,OW_sfc,wspd_sfc,conv_sfc
    
    ### 1 km MV variables ###
    # ds = nc.Dataset(fp+f"cm1out_{fn:06d}.nc")
    # w_1km = ds.variables['winterp'][:].data[0,iz1,iy,ix]
    # zvort_1km = ds.variables['zvort'][:].data[0,iz1,iy,ix]
    # u_1km = ds.variables['uinterp'][:].data[0,iz1,iy,ix]
    # v_1km = ds.variables['vinterp'][:].data[0,iz1,iy,ix]
    # S_N = np.gradient(u_1km, xh[ix]*1000, axis=1) - np.gradient(v_1km, yh[iy]*1000, axis=0)
    # S_S = np.gradient(v_1km, xh[ix]*1000, axis=1) + np.gradient(u_1km, yh[iy]*1000, axis=0)
    # OW_1km = S_N**2 + S_S**2 - zvort_1km**2
    # wspd_1km = (u_1km**2 + v_1km**2)**0.5
    # del u_1km,v_1km,S_N,S_S
    # ds.close()
    
    # w_1km_swath = np.maximum(w_1km_swath, w_1km)
    # zvort_1km_swath = np.maximum(zvort_1km_swath, zvort_1km)
    # OW_1km_swath = np.minimum(OW_1km_swath, OW_1km)
    # wspd_1km_swath = np.maximum(wspd_1km_swath, wspd_1km)
    # del w_1km,zvort_1km,OW_1km,wspd_1km
    
    
    # ds = nc.Dataset(fp+f"cm1out_{fn:06d}.nc")
    # zvort_2km = ds.variables['zvort'][:].data[0,iz2,iy,ix]
    # u_2km = ds.variables['uinterp'][:].data[0,iz2,iy,ix]
    # v_2km = ds.variables['vinterp'][:].data[0,iz2,iy,ix]
    # S_N = np.gradient(u_2km, xh[ix]*1000, axis=1) - np.gradient(v_2km, yh[iy]*1000, axis=0)
    # S_S = np.gradient(v_2km, xh[ix]*1000, axis=1) + np.gradient(u_2km, yh[iy]*1000, axis=0)
    # OW_2km = S_N**2 + S_S**2 - zvort_2km**2
    # del u_2km,v_2km,S_N,S_S,zvort_2km
    # ds.close()
    
    # OW_2km_swath = np.minimum(OW_2km_swath, OW_2km)
    # del OW_2km
    
    ds = nc.Dataset(fp+f"cm1out_{fn:06d}.nc")
    zvort_3km = ds.variables['zvort'][:].data[0,iz3,iy,ix]
    u_3km = ds.variables['uinterp'][:].data[0,iz3,iy,ix]
    v_3km = ds.variables['vinterp'][:].data[0,iz3,iy,ix]
    S_N = np.gradient(u_3km, xh[ix]*1000, axis=1) - np.gradient(v_3km, yh[iy]*1000, axis=0)
    S_S = np.gradient(v_3km, xh[ix]*1000, axis=1) + np.gradient(u_3km, yh[iy]*1000, axis=0)
    OW_3km = S_N**2 + S_S**2 - zvort_3km**2
    del u_3km,v_3km,S_N,S_S,zvort_3km
    ds.close()
    
    OW_3km_swath = np.minimum(OW_3km_swath, OW_3km)
    del OW_3km
    
    ### Updraft helicity ###
    # ds = nc.Dataset(fp+f"cm1out_{fn:06d}.nc")
    # u = ds.variables['uinterp'][:].data[0,:,iy,ix]
    # v = ds.variables['vinterp'][:].data[0,:,iy,ix]
    # w = ds.variables['w'][:].data[0,:,iy,ix]
    # ds.close()
    
    # zf3d = np.tile(np.expand_dims(zf,(1,2)), (1,len(yh[iy]),len(xh[ix])))
    # mapfct = np.ones(shape=(len(yh[iy]),len(xh[ix])), dtype=float)
    # UH_25km = wrf.to_np(wrf.udhel(zf3d*1000, mapfct, u, v, w, 125, 125, bottom=2000., top=5000.))
    # UH_02km = wrf.to_np(wrf.udhel(zf3d*1000, mapfct, u, v, w, 125, 125, bottom=0., top=2000.))
    # del u,v,w,zf3d,mapfct
    # UH_25km = np.asarray(UH_25km.data)
    # UH_02km = np.asarray(UH_02km.data)
    # UH_25_swath = np.maximum(UH_25_swath, UH_25km)
    # UH_02_swath = np.maximum(UH_02_swath, UH_02km)
    # del UH_25km,UH_02km
    
    ### Pressure perturbation and VPPGA ###
    # ds = nc.Dataset(fp+f"base/cm1out_000013.nc")
    # prs0 = ds.variables['prs0'][:].data[0,slice(0,iz1+1),iy,ix]
    # ds.close()
    # ds = nc.Dataset(fp+f"cm1out_{fn:06d}.nc")
    # prspert = ds.variables['prs'][:].data[0,slice(0,iz1+1),iy,ix] - prs0
    # ds.close()
    
    # pp_sfc = prspert[0,:,:]
    # pp_1km = prspert[iz1,:,:]
    # pp_1km_mean = np.mean(prspert, axis=0)
    # dp = pp_1km - pp_sfc
    # dz = (z[iz1] - z[0])*1000
    # vppga = -1/1.1 * dp / dz
    # pp_sfc_swath = np.minimum(pp_sfc_swath, pp_sfc)
    # pp_1km_swath = np.minimum(pp_1km_swath, pp_1km_mean)
    # vppga_max_swath = np.maximum(vppga_max_swath, vppga)
    # vppga_min_swath = np.minimum(vppga_min_swath, vppga)
    # del prs0,prspert,pp_sfc,pp_1km,pp_1km_mean,dp,vppga
    
    
    # ds = nc.Dataset(fp+f"pp/dyn_{fn:06d}.nc")
    # dz = (z[iz1] - z[0])*1000
    # p_dn_sfc = ds.variables['p_dn'][:].data[0,iy,ix]
    # p_dn_1km = np.mean(ds.variables['p_dn'][:].data[0:iz1,iy,ix], axis=0)
    # dp_dn = ds.variables['p_dn'][:].data[iz1,iy,ix] - ds.variables['p_dn'][:].data[0,iy,ix]
    # # p_b_sfc = ds.variables['p_b'][:].data[0,iy,ix]
    # # p_b_1km = np.mean(ds.variables['p_b'][:].data[0:iz1,iy,ix], axis=0)
    # # dp_b = ds.variables['p_b'][:].data[iz1,iy,ix] - ds.variables['p_b'][:].data[0,iy,ix]
    # ds.close()
    
    # p_dn_sfc_swath = np.minimum(p_dn_sfc_swath, p_dn_sfc)
    # p_dn_1km_swath = np.minimum(p_dn_1km_swath, p_dn_1km)
    # vppga_dn = -1/1.1 * dp_dn / dz
    # vppga_dn_max_swath = np.maximum(vppga_dn_max_swath, vppga_dn)
    # vppga_dn_min_swath = np.minimum(vppga_dn_min_swath, vppga_dn)
    # del p_dn_sfc,p_dn_1km,dp_dn,vppga_dn
    
    # # p_b_sfc_swath = np.minimum(p_b_sfc_swath, p_b_sfc)
    # # p_b_1km_swath = np.minimum(p_b_1km_swath, p_b_1km)
    # # vppga_b = -1/1.1 * dp_b / dz
    # # vppga_b_swath = np.minimum(vppga_b_swath, vppga_b)
    # # del p_b_sfc,p_b_1km,dp_b,vppga_b


#%% Save or update pickle file
save_pkl = False
update_pkl = False

if save_pkl:
    swaths = {'ix':[ix1,ix2], 'iy':[iy1,iy2], 'w_2km':w_2km_swath, 'zvort_sfc':zvort_sfc_swath, 
              'OW_sfc':OW_sfc_swath, 'wspd_sfc':wspd_sfc_swath, 'conv_sfc':conv_sfc_swath}
    dbfile = open(ip+'swaths.pkl', 'wb')
    pickle.dump(swaths, dbfile)
    dbfile.close()
    
if update_pkl:
    dbfile = open(ip+'swaths.pkl', 'rb')
    swaths = pickle.load(dbfile)
    dbfile.close()
    
    new_vars = {'OW_3km':OW_3km_swath}
    swaths.update(new_vars)
    dbfile = open(ip+'swaths.pkl', 'wb')
    pickle.dump(swaths, dbfile)
    dbfile.close()



#%% Load swaths from saved pickle file

fp = '/Volumes/Promise_Pegasus_70TB/merger/supercell-125m/'
ip = '/Users/morgan.schneider/Documents/merger/supercell-125m/'
sim = 'SUPERCELL'

ds = nc.Dataset(fp+'cm1out_000014.nc')
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
ds.close()

dbfile = open(ip+'swaths.pkl', 'rb')
swaths = pickle.load(dbfile)
ix = slice(swaths['ix'][0], swaths['ix'][1])
iy = slice(swaths['iy'][0], swaths['iy'][1])
w_2km_swath = swaths['w_2km']
zvort_sfc_swath = swaths['zvort_sfc']
OW_sfc_swath = swaths['OW_sfc']
wspd_sfc_swath = swaths['wspd_sfc']
conv_sfc_swath = swaths['conv_sfc']
w_1km_swath = swaths['w_1km']
zvort_1km_swath = swaths['zvort_1km']
OW_1km_swath = swaths['OW_1km']
OW_2km_swath = swaths['OW_2km']
OW_3km_swath = swaths['OW_3km']
wspd_1km_swath = swaths['wspd_1km']
UH_25_swath = swaths['UH_25km']
UH_02_swath= swaths['UH_02km']
pp_sfc_swath = swaths['pp_sfc']
pp_1km_swath = swaths['pp_1km']
vppga_max_swath = swaths['vppga_max']
vppga_min_swath = swaths['vppga_min']

# p_dn_sfc_swath = swaths['pdn_sfc']
# p_dn_1km_swath = swaths['pdn_1km_mean']
# vppga_dn_max_swath = swaths['vppga_dn_max']
# vppga_dn_min_swath = swaths['vppga_dn_min']
# p_b_sfc_swath = swaths['pb_sfc']
# p_b_1km_swath = swaths['pb_1km_mean']
# vppga_b_swath = swaths['vppga_b_min']
dbfile.close()

#%% Plot swaths from saved pickle files


figsave = False

# xl = [-60,20]
# yl = [-100,-20]
xl = [-60,30]
yl = [-110,-20]


if True:
    cm = plt.get_cmap('Blues_r')
    cm.set_over('white')
    
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    plot_cfill(xh[ix], yh[iy], np.ma.masked_array(OW_2km_swath, OW_2km_swath>-0.001), 'OW', ax, datalims=[-0.003,0], cmap='Blues_r', xlims=xl, ylims=yl)
    # ax.contour(xh[ix], yh[iy], gaussian_filter(w_2km_swath,2.5), levels=[5,10], colors='k', linewidths=[1,1.5])
    ax.set_title(f"{sim} - Minimum 2 km OW")
    
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    plot_contourf(xh[ix], yh[iy], OW_2km_swath, 'OW', ax, levels=np.linspace(-0.001,0,21), datalims=[-0.001,-0.0001], cmap=cm, xlims=xl, ylims=yl, extend='min')
    ax.set_title(f"{sim} - Minimum 2 km OW")
    
    
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    plot_cfill(xh[ix], yh[iy], np.ma.masked_array(OW_3km_swath, OW_3km_swath>-0.001), 'OW', ax, datalims=[-0.003,0], cmap='Blues_r', xlims=xl, ylims=yl)
    # ax.contour(xh[ix], yh[iy], gaussian_filter(w_2km_swath,2.5), levels=[5,10], colors='k', linewidths=[1,1.5])
    ax.set_title(f"{sim} - Minimum 3 km OW")
    
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    plot_contourf(xh[ix], yh[iy], OW_3km_swath, 'OW', ax, levels=np.linspace(-0.001,0,21), datalims=[-0.001,-0.0001], cmap=cm, xlims=xl, ylims=yl, extend='min')
    ax.set_title(f"{sim} - Minimum 3 km OW")
#%%

if False:
    # 2 km max w
    fig,ax = plt.subplots(1,1,figsize=(9,6))
    plot_cfill(xh[ix], yh[iy], np.ma.masked_array(w_2km_swath, w_2km_swath<5), 'w', ax, datalims=[5,20], cmap='Reds', xlims=xl, ylims=yl)
    # ax.contour(xh[ix], yh[iy], gaussian_filter(w_2km_swath,2.5), levels=[5,10], colors='k', linewidths=[1,1.5])
    # plot_cfill(xh[ix], yh[iy], np.ma.masked_array(OW_sfc_swath, OW_sfc_swath>-0.001), 'OW', ax, datalims=[-0.001,-0.003], cmap=cmocean.cm.ice)
    ax.set_title(f"{sim} - Maximum 2 km w")
    if figsave:
        plt.savefig(ip+'imgs_swaths/swaths_w_2km.png', dpi=300)
    
    # 1 km max w
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    plot_cfill(xh[ix], yh[iy], np.ma.masked_array(w_1km_swath, w_1km_swath<5), 'w', ax, datalims=[5,20], cmap='Reds', xlims=xl, ylims=yl)
    # ax.contour(xh[ix], yh[iy], gaussian_filter(w_1km_swath,2.5), levels=[5,10], colors='k', linewidths=[1,1.5])
    # plot_cfill(xh[ix], yh[iy], np.ma.masked_array(OW_1km_swath, OW_1km_swath>-0.001), 'OW', ax, datalims=[-0.001,-0.003], cmap=cmocean.cm.ice)
    ax.set_title(f"{sim} - Maximum 1 km w")
    if figsave:
        plt.savefig(ip+'imgs_swaths/swaths_w_1km.png', dpi=300)

    # sfc max zvort
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    plot_cfill(xh[ix], yh[iy], np.ma.masked_array(zvort_sfc_swath, zvort_sfc_swath<0.02), 'zvort', ax, datalims=[0.02,0.08], cmap='Reds', xlims=xl, ylims=yl)
    # ax.contour(xh[ix], yh[iy], gaussian_filter(w_1km_swath,2.5), levels=[5,10], colors='k', linewidths=[1,1.5])
    ax.set_title(f"{sim} - Maximum surface \u03B6")
    if figsave:
        plt.savefig(ip+'imgs_swaths/swaths_zvort_sfc.png', dpi=300)
    
    # 1 km max zvort
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    plot_cfill(xh[ix], yh[iy], gaussian_filter(np.ma.masked_array(zvort_1km_swath, zvort_1km_swath<0.02),1), 'zvort', ax, datalims=[0.02,0.08], cmap='Reds', xlims=xl, ylims=yl)
    # ax.contour(xh[ix], yh[iy], gaussian_filter(w_1km_swath,2.5), levels=[5,10], colors='k', linewidths=[1,1.5])
    ax.set_title(f"{sim} - Maximum 1 km \u03B6")
    if figsave:
        plt.savefig(ip+'imgs_swaths/swaths_zvort_1km.png', dpi=300)

    # sfc min OW
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    plot_cfill(xh[ix], yh[iy], np.ma.masked_array(OW_sfc_swath, OW_sfc_swath>-0.001), 'OW', ax, datalims=[-0.003,0], cmap='Blues_r', xlims=xl, ylims=yl)
    # ax.contour(xh[ix], yh[iy], gaussian_filter(w_1km_swath,2.5), levels=[5,10], colors='k', linewidths=[1,1.5])
    ax.set_title(f"{sim} - Minimum surface OW")
    if figsave:
        plt.savefig(ip+'imgs_swaths/swaths_OW_sfc.png', dpi=300)
    
    # 1 km min OW
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    plot_cfill(xh[ix], yh[iy], gaussian_filter(np.ma.masked_array(OW_1km_swath, OW_1km_swath>-0.001),1), 'OW', ax, datalims=[-0.003,0], cmap='Blues_r', xlims=xl, ylims=yl)
    # ax.contour(xh[ix], yh[iy], gaussian_filter(w_1km_swath,2.5), levels=[5,10], colors='k', linewidths=[1,1.5])
    ax.set_title(f"{sim} - Minimum 1 km OW")
    if figsave:
        plt.savefig(ip+'imgs_swaths/swaths_OW_1km.png', dpi=300)
    
    # 2 km min OW
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    plot_cfill(xh[ix], yh[iy], gaussian_filter(np.ma.masked_array(OW_2km_swath, OW_2km_swath>-0.001),1), 'OW', ax, datalims=[-0.003,0], cmap='Blues_r', xlims=xl, ylims=yl)
    # ax.contour(xh[ix], yh[iy], gaussian_filter(w_1km_swath,2.5), levels=[5,10], colors='k', linewidths=[1,1.5])
    ax.set_title(f"{sim} - Minimum 2 km OW")
    if figsave:
        plt.savefig(ip+'imgs_swaths/swaths_OW_2km.png', dpi=300)
    
    # sfc max wind speed
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    plot_cfill(xh[ix], yh[iy], np.ma.masked_array(wspd_sfc_swath, wspd_sfc_swath<15), 'wspd', ax, datalims=[15,35], xlims=xl, ylims=yl)
    # ax.contour(xh[ix], yh[iy], gaussian_filter(w_1km_swath,2.5), levels=[5,10], colors='k', linewidths=[1,1.5])
    ax.set_title(f"{sim} - Maximum surface wind speed")
    if figsave:
        plt.savefig(ip+'imgs_swaths/swaths_wspd_sfc.png', dpi=300)
    
    # 1 km max wind speed
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    plot_cfill(xh[ix], yh[iy], np.ma.masked_array(wspd_1km_swath, wspd_1km_swath<15), 'wspd', ax, datalims=[15,35], xlims=xl, ylims=yl)
    # ax.contour(xh[ix], yh[iy], gaussian_filter(w_1km_swath,2.5), levels=[5,10], colors='k', linewidths=[1,1.5])
    ax.set_title(f"{sim} - Maximum 1 km wind speed")
    if figsave:
        plt.savefig(ip+'imgs_swaths/swaths_wspd_1km.png', dpi=300)
    
    # sfc min convergence
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    plot_cfill(xh[ix], yh[iy], gaussian_filter(np.ma.masked_array(conv_sfc_swath, conv_sfc_swath>-0.02),1), 'divh', ax, datalims=[-0.08,0], cmap='Blues_r', xlims=xl, ylims=yl)
    # ax.contour(xh[ix], yh[iy], gaussian_filter(w_1km_swath,2.5), levels=[5,10], colors='k', linewidths=[1,1.5])
    # ax.contour(xh[ix], yh[iy], p_dn_1km_swath/100, levels=[-2], colors='k', linewidths=1, linestyles='-')
    ax.set_title(f"{sim} - Minimum surface convergence")
    if figsave:
        plt.savefig(ip+'imgs_swaths/swaths_convergence_sfc.png', dpi=300)
    
    # 2-5 km max UH
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    plot_cfill(xh[ix], yh[iy], np.ma.masked_array(UH_25_swath, UH_25_swath<50), 'uh', ax, datalims=[50,500], xlims=xl, ylims=yl)
    # ax.contour(xh[ix], yh[iy], gaussian_filter(w_2km_swath,3), levels=[5,10], colors='k', linewidths=[1,1.5])
    ax.set_title(f"{sim} - Maximum 2-5 km UH")
    if figsave:
        plt.savefig(ip+'imgs_swaths/swaths_UH_2-5.png', dpi=300)
    
    # 0-2 km max UH
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    plot_cfill(xh[ix], yh[iy], np.ma.masked_array(UH_02_swath, UH_02_swath<50), 'uh', ax, datalims=[50,250], xlims=xl, ylims=yl)
    # ax.contour(xh[ix], yh[iy], gaussian_filter(w_1km_swath,2.5), levels=[5,10], colors='k', linewidths=[1,1.5])
    ax.set_title(f"{sim} - Maximum 0-2 km UH")
    if figsave:
        plt.savefig(ip+'imgs_swaths/swaths_UH_0-2.png', dpi=300)
        

if False:
    # sfc min pressure perturbation
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    plot_cfill(xh[ix], yh[iy], np.ma.masked_array(pp_sfc_swath/100, pp_sfc_swath>=0), 'prspert', ax, datalims=[-2,0], cmap='Blues_r', xlims=xl, ylims=yl)
    # ax.contour(xh[ix], yh[iy], gaussian_filter(w_1km_swath,2.5), levels=[5,10], colors='k', linewidths=[1,1.5])
    ax.set_title(f"{sim} - Minimum surface p'")
    if figsave:
        plt.savefig(ip+'imgs_swaths/swaths_prspert_sfc.png', dpi=300)
    
    # 1 km min pressure perturbation
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    plot_cfill(xh[ix], yh[iy], np.ma.masked_array(pp_1km_swath/100, pp_1km_swath>=0), 'prspert', ax, datalims=[-2,0], cmap='Blues_r', xlims=xl, ylims=yl)
    # ax.contour(xh[ix], yh[iy], gaussian_filter(w_1km_swath,2.5), levels=[5,10], colors='k', linewidths=[1,1.5])
    ax.set_title(f"{sim} - Minimum 0-1 km mean p'")
    if figsave:
        plt.savefig(ip+'imgs_swaths/swaths_prspert_1km.png', dpi=300)
    
    # 0-1 km max upward VPPGA
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    plot_cfill(xh[ix], yh[iy], np.ma.masked_array(vppga_max_swath, vppga_max_swath<0.2), 'pgfz', ax, datalims=[0.2,0.4], cmap='Reds', xlims=xl, ylims=yl)
    # ax.contour(xh[ix], yh[iy], gaussian_filter(w_1km_swath,2.5), levels=[5,10], colors='k', linewidths=[1,1.5])
    ax.set_title(f"{sim} - Maximum upward 0-1 km VPPGA")
    if figsave:
        plt.savefig(ip+'imgs_swaths/swaths_vppga_max.png', dpi=300)
    
    # # 0-1 km min downward VPPGA
    # fig,ax = plt.subplots(1,1,figsize=(8,6))
    # plot_cfill(xh[ix], yh[iy], np.ma.masked_array(vppga_min_swath, vppga_min_swath>0), 'pgfz', ax, datalims=[-0.05,0], cmap='Blues_r', xlims=xl, ylims=yl)
    # # ax.contour(xh[ix], yh[iy], gaussian_filter(w_1km_swath,2.5), levels=[5,10], colors='k', linewidths=[1,1.5])
    # ax.set_title(f"{sim} - Minimum downward 0-1 km VPPGA")
    # if figsave:
    #     plt.savefig(ip+'imgs_swaths/swaths_vppga_min.png', dpi=300)

if False:
    # Sfc p_dn
    fig,ax = plt.subplots(1,1,figsize=(9,6))
    plot_cfill(xh[ix], yh[iy], p_dn_sfc_swath/100, 'prspert', ax, datalims=[-5,0], cmap='Blues_r', xlims=xl, ylims=yl)
    ax.contour(xh[ix], yh[iy], w_2km_swath, levels=[10], colors='grey', linewidths=1)
    plot_cfill(xh[ix], yh[iy], np.ma.masked_array(OW_sfc_swath, OW_sfc_swath>-0.001), 'OW', ax, datalims=[-0.001,-0.003], cmap='Reds_r')
    ax.set_title(f"{sim} - Sfc p' (2 km w, sfc OW)")
    if figsave:
        plt.savefig(ip+'imgs_swaths/swaths_pdn_sfc.png', dpi=300)
    
    # 0-1 km mean p_dn
    fig,ax = plt.subplots(1,1,figsize=(9,6))
    plot_cfill(xh[ix], yh[iy], p_dn_1km_swath/100, 'prspert', ax, datalims=[-5,0], cmap='Blues_r', xlims=xl, ylims=yl)
    ax.contour(xh[ix], yh[iy], w_2km_swath, levels=[10], colors='grey', linewidths=1)
    plot_cfill(xh[ix], yh[iy], np.ma.masked_array(OW_sfc_swath, OW_sfc_swath>-0.001), 'OW', ax, datalims=[-0.001,-0.003], cmap='Reds_r')
    ax.set_title(f"{sim} - 0-1 km mean p' (2 km w, sfc OW)")
    if figsave:
        plt.savefig(ip+'imgs_swaths/swaths_pdn_1km.png', dpi=300)

    # 0-1 km upward VPPGA dn
    fig,ax = plt.subplots(1,1,figsize=(9,6))
    plot_cfill(xh[ix], yh[iy], np.ma.masked_array(vppga_dn_max_swath, vppga_dn_max_swath<0.02), 'pgfz', ax, datalims=[0,0.2], cmap='Reds', xlims=xl, ylims=yl)
    # plot_cfill(xh[ix], yh[iy], vppga_dn_max_swath, 'pgfz', ax, datalims=[0,0.2], cmap='Reds', xlims=xl, ylims=yl)
    # ax.contour(xh[ix], yh[iy], w_2km_swath, levels=[10], colors='k', linewidths=1)
    plot_cfill(xh[ix], yh[iy], np.ma.masked_array(OW_sfc_swath, OW_sfc_swath>-0.001), 'OW', ax, datalims=[-0.001,-0.003], cmap=cmocean.cm.ice)
    ax.set_title(f"{sim} - 0-1 km upward VPPGA (nonlinear dynamic), sfc OW")
    if figsave:
        plt.savefig(ip+'imgs_swaths/swaths_vppga_dn_max.png', dpi=300)

    # 0-1 km downward VPPGA dn
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    plot_cfill(xh[ix], yh[iy], np.ma.masked_array(vppga_dn_min_swath, vppga_dn_min_swath>-0.01), 'pgfz', ax, datalims=[-0.1,0], cmap='Blues_r', xlims=xl, ylims=yl)
    # plot_cfill(xh[ix], yh[iy], vppga_dn_min_swath, 'pgfz', ax, datalims=[-0.1,0], cmap='Blues_r', xlims=xl, ylims=yl)
    # ax.contour(xh[ix], yh[iy], wspd_sfc_swath, levels=[22], colors='k', linewidths=1)
    # ax.contour(xh[ix], yh[iy], w_2km_swath, levels=[10], colors='k', linewidths=1)
    ax.set_title(f"{sim} - 0-1 km downward VPPGA (nonlinear dynamic)")
    if figsave:
        plt.savefig(ip+'imgs_swaths/swaths_vppga_dn_min.png', dpi=300)



#%% Multi panel swath figure

from CM1utils import *

ds = nc.Dataset('/Volumes/Promise_Pegasus_70TB/merger/merger-125m/cm1out_000014.nc')
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
ds.close()

dbfile = open('/Users/morgan.schneider/Documents/merger/merger-125m/swaths.pkl', 'rb')
swaths = pickle.load(dbfile)
ix = slice(swaths['ix'][0], swaths['ix'][1])
iy = slice(swaths['iy'][0], swaths['iy'][1])
# OW_sfc_m = swaths['OW_sfc']
wspd_sfc_m = swaths['wspd_sfc']
w_1km_m = swaths['w_1km']
OW_1km_m = swaths['OW_1km']
UH_25_m = swaths['UH_25km']
UH_02_m = swaths['UH_02km']
pp_1km_m = swaths['pp_1km']
vppga_max_m = swaths['vppga_max']
dbfile.close()

dbfile = open('/Users/morgan.schneider/Documents/merger/qlcs-125m/swaths.pkl', 'rb')
swaths = pickle.load(dbfile)
ix = slice(swaths['ix'][0], swaths['ix'][1])
iy = slice(swaths['iy'][0], swaths['iy'][1])
# OW_sfc_q = swaths['OW_sfc']
wspd_sfc_q = swaths['wspd_sfc']
w_1km_q = swaths['w_1km']
OW_1km_q = swaths['OW_1km']
UH_25_q = swaths['UH_25km']
UH_02_q = swaths['UH_02km']
pp_1km_q = swaths['pp_1km']
vppga_max_q = swaths['vppga_max']
dbfile.close()

dbfile = open('/Users/morgan.schneider/Documents/merger/supercell-125m/swaths.pkl', 'rb')
swaths = pickle.load(dbfile)
ix = slice(swaths['ix'][0], swaths['ix'][1])
iy = slice(swaths['iy'][0], swaths['iy'][1])
# OW_sfc_s = swaths['OW_sfc']
wspd_sfc_s = swaths['wspd_sfc']
w_1km_s = swaths['w_1km']
OW_1km_s = swaths['OW_1km']
UH_25_s = swaths['UH_25km']
UH_02_s = swaths['UH_02km']
pp_1km_s = swaths['pp_1km']
vppga_max_s = swaths['vppga_max']
dbfile.close()




figsave = False

xl = [-60,30]
yl = [-110,-20]


fig,axs = plt.subplots(3, 3, figsize=(8,7), subplot_kw=dict(box_aspect=1), sharex=True, sharey=True, layout='constrained')

plot_cfill(xh[ix], yh[iy], np.ma.masked_array(w_1km_m, w_1km_m<5), 'w', axs[0,0], datalims=[5,15], cmap='Reds', xlims=xl, ylims=yl, cbar=False)
# axs[0,0].set_title("MERGER")
plot_cfill(xh[ix], yh[iy], np.ma.masked_array(w_1km_q, w_1km_q<5), 'w', axs[0,1], datalims=[5,15], cmap='Reds', xlims=xl, ylims=yl, cbar=False)
# axs[0,1].set_title("QLCS")
plot_cfill(xh[ix], yh[iy], np.ma.masked_array(w_1km_s, w_1km_s<5), 'w', axs[0,2], datalims=[5,15], cmap='Reds', xlims=xl, ylims=yl)
# axs[0,2].set_title("SUPERCELL")

plot_cfill(xh[ix], yh[iy], gaussian_filter(np.ma.masked_array(OW_1km_m, OW_1km_m>-0.001),1), 'OW', axs[1,0], datalims=[-0.003,0], cmap='Blues_r', xlims=xl, ylims=yl, cbar=False)
plot_cfill(xh[ix], yh[iy], gaussian_filter(np.ma.masked_array(OW_1km_q, OW_1km_q>-0.001),1), 'OW', axs[1,1], datalims=[-0.003,0], cmap='Blues_r', xlims=xl, ylims=yl, cbar=False)
plot_cfill(xh[ix], yh[iy], gaussian_filter(np.ma.masked_array(OW_1km_s, OW_1km_s>-0.001),1), 'OW', axs[1,2], datalims=[-0.003,0], cmap='Blues_r', xlims=xl, ylims=yl)
axs[1,0].set_ylabel('y distance (km)')

plot_cfill(xh[ix], yh[iy], np.ma.masked_array(wspd_sfc_m, wspd_sfc_m<15), 'wspd', axs[2,0], datalims=[15,35], xlims=xl, ylims=yl, cbar=False)
plot_cfill(xh[ix], yh[iy], np.ma.masked_array(wspd_sfc_q, wspd_sfc_q<15), 'wspd', axs[2,1], datalims=[15,35], xlims=xl, ylims=yl, cbar=False)
plot_cfill(xh[ix], yh[iy], np.ma.masked_array(wspd_sfc_s, wspd_sfc_s<15), 'wspd', axs[2,2], datalims=[15,35], xlims=xl, ylims=yl)
axs[2,1].set_xlabel('x distance (km)')

if figsave:
    plt.savefig('/Users/morgan.schneider/Documents/merger/swaths_w-OW-wspd.png', dpi=300)



fig,axs = plt.subplots(2, 3, figsize=(9,5.5), subplot_kw=dict(box_aspect=1), sharex=True, sharey=True, layout='constrained')

plot_cfill(xh[ix], yh[iy], np.ma.masked_array(UH_02_m, UH_02_m<50), 'uh', axs[0,0], datalims=[50,250], xlims=xl, ylims=yl, cbar=False)
plot_cfill(xh[ix], yh[iy], np.ma.masked_array(UH_02_q, UH_02_q<50), 'uh', axs[0,1], datalims=[50,250], xlims=xl, ylims=yl, cbar=False)
plot_cfill(xh[ix], yh[iy], np.ma.masked_array(UH_02_s, UH_02_s<50), 'uh', axs[0,2], datalims=[50,250], xlims=xl, ylims=yl)
axs[0,0].set_ylabel('y distance (km)')

plot_cfill(xh[ix], yh[iy], np.ma.masked_array(UH_25_m, UH_25_m<50), 'uh', axs[1,0], datalims=[50,500], xlims=xl, ylims=yl, cbar=False)
plot_cfill(xh[ix], yh[iy], np.ma.masked_array(UH_25_q, UH_25_q<50), 'uh', axs[1,1], datalims=[50,500], xlims=xl, ylims=yl, cbar=False)
plot_cfill(xh[ix], yh[iy], np.ma.masked_array(UH_25_s, UH_25_s<50), 'uh', axs[1,2], datalims=[50,500], xlims=xl, ylims=yl)
axs[1,0].set_ylabel('y distance (km)')
axs[1,1].set_xlabel('x distance (km)')

if figsave:
    plt.savefig('/Users/morgan.schneider/Documents/merger/swaths_UH.png', dpi=300)



fig,axs = plt.subplots(2, 3, figsize=(9,5.5), subplot_kw=dict(box_aspect=1), sharex=True, sharey=True, layout='constrained')

plot_cfill(xh[ix], yh[iy], pp_1km_m/100, 'prspert', axs[0,0], datalims=[-3,0], cmap='Blues_r', xlims=xl, ylims=yl, cbar=False)
plot_cfill(xh[ix], yh[iy], pp_1km_q/100, 'prspert', axs[0,1], datalims=[-3,0], cmap='Blues_r', xlims=xl, ylims=yl, cbar=False)
plot_cfill(xh[ix], yh[iy], pp_1km_s/100, 'prspert', axs[0,2], datalims=[-3,0], cmap='Blues_r', xlims=xl, ylims=yl)
axs[0,0].set_ylabel('y distance (km)')

plot_cfill(xh[ix], yh[iy], vppga_max_m, 'pgfz', axs[1,0], datalims=[0,0.5], cmap='Reds', xlims=xl, ylims=yl, cbar=False)
plot_cfill(xh[ix], yh[iy], vppga_max_q, 'pgfz', axs[1,1], datalims=[0,0.5], cmap='Reds', xlims=xl, ylims=yl, cbar=False)
plot_cfill(xh[ix], yh[iy], vppga_max_s, 'pgfz', axs[1,2], datalims=[0,0.5], cmap='Reds', xlims=xl, ylims=yl)
axs[1,1].set_xlabel('x distance (km)')
axs[1,0].set_ylabel('y distance (km)')

if figsave:
    plt.savefig('/Users/morgan.schneider/Documents/merger/swaths_prspert+vppga.png', dpi=300)



