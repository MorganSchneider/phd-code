#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 16:24:06 2024

@author: morgan.schneider
"""

from CM1utils import *
import wrf
import metpy.calc as mc
from metpy.units import units


#%% Calculate environmental fields

fp = '/Volumes/Promise_Pegasus_70TB/merger/merger-125m/' 
ip = '/Users/morgan.schneider/Documents/merger/merger-125m/'   

fn = 14
# 181 min -> 14 | 195 min -> 28 | 210 min -> 43 | 225 min -> 58 | 240 min -> 73                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     

ds = nc.Dataset(fp+'cm1out_000014.nc')
time = ds.variables['time'][:].data[0]/60
xh = ds.variables['xh'][:].data
xf = ds.variables['xf'][:].data
yh = ds.variables['yh'][:].data
yf = ds.variables['yf'][:].data
z = ds.variables['z'][:].data
zf = ds.variables['zf'][:].data
ds.close()



xlims = [-70,10]
ylims = [-140,-60]

# ix1 = np.where(xh >= xlims[0])[0][0]
# ix2 = np.where(xh >= xlims[1])[0][0]
# iy1 = np.where(yh >= ylims[0])[0][0]
# iy2 = np.where(yh >= ylims[1])[0][0]
ix1 = 0; ix2 = len(xh)
iy1 = 0; iy2 = len(yh)
ix = slice(ix1,ix2)
iy = slice(iy1,iy2)


# Calculate environmental fields using wrf package
calc_CAPE = True
calc_SRH = True
calc_UH = True

z3d = np.tile(np.expand_dims(z,(1,2)), (1,len(yh),len(xh)))
ter = np.zeros(shape=(len(yh),len(xh)), dtype=float)

if calc_CAPE:
    if 'prs' not in locals():
        print('Loading data...')
        ds = nc.Dataset(fp+f"cm1out_{fn:06d}.nc")
        th = ds.variables['th'][:].data[0,:,:,:]
        prs = ds.variables['prs'][:].data[0,:,:,:]
        qv = ds.variables['qv'][:].data[0,:,:,:]
        T = th * (prs/100000.)**0.286
        del th
        ds.close()
    
    print('Calculating CAPE...')
    tmp = wrf.to_np(wrf.cape_2d(prs/100, T, qv, z3d*1000, ter, prs[0,:,:]/100, False))
    del prs,qv,T
    MUcape = np.asarray(tmp[0,:,:].data)
    MUcin = np.asarray(tmp[1,:,:].data)
    z_LCL = np.asarray(tmp[2,:,:].data)
    # z_LFC = tmp[3,:,:]
    del tmp
    
    dat = {'MUcape':MUcape, 'MUcin':MUcin, 'z_LCL':z_LCL}
    # del MUcape,MUcin,z_LCL

if calc_SRH:
    if 'u' not in locals():
        print('Loading data...')
        ds = nc.Dataset(fp+f"cm1out_{fn:06d}.nc")
        u = ds.variables['uinterp'][:].data[0,:,:,:] - 4
        v = ds.variables['vinterp'][:].data[0,:,:,:] - 23
        ds.close()
    
    print('Calculating 0-3 km SRH...')
    SRH_3km = wrf.to_np(wrf.srhel(u, v, z3d*1000, ter, top=3000.))
    print('Calculating 0-1 km SRH...')
    SRH_1km = wrf.to_np(wrf.srhel(u, v, z3d*1000, ter, top=1000.))
    
    SRH_3km = np.asarray(SRH_3km.data)
    SRH_1km = np.asarray(SRH_1km.data)
    
    print('Calculating bulk shear...')
    iz6 = np.where(z >= 6)[0][0]
    u_shr = u[iz6,:,:] - u[0,:,:]
    v_shr = v[iz6,:,:] - v[0,:,:]
    shear_6km = np.sqrt(u_shr**2 + v_shr**2)
    del u,v
    
    dat = {'SRH_3km':SRH_3km, 'SRH_1km':SRH_1km, 'shear_6km':shear_6km, 'u_shr':u_shr, 'v_shr':v_shr}
    # del SRH_3km,SRH_1km,shear_6km,u_shr,v_shr

if calc_UH:
    if 'w' not in locals():
        print('Loading data...')
        ds = nc.Dataset(fp+f"cm1out_{fn:06d}.nc")
        u = ds.variables['uinterp'][:].data[0,:,:,:]
        v = ds.variables['vinterp'][:].data[0,:,:,:]
        w = ds.variables['w'][:].data[0,:,:,:]
        ds.close()
    zf3d = np.tile(np.expand_dims(zf,(1,2)), (1,len(yh),len(xh)))
    mapfct = np.ones(shape=(len(yh),len(xh)), dtype=float)
    
    print('Calculating 2-5 km UH...')
    UH_25km = wrf.to_np(wrf.udhel(zf3d*1000, mapfct, u, v, w, 125, 125, bottom=2000., top=5000.))
    print('Calculating 0-2 km UH...')
    UH_02km = wrf.to_np(wrf.udhel(zf3d*1000, mapfct, u, v, w, 125, 125, bottom=0., top=2000.))
    del u,v,w,mapfct
    UH_25km = np.asarray(UH_25km.data)
    UH_02km = np.asarray(UH_02km.data)
    
    dat = {'UH_25km':UH_25km, 'UH_02km':UH_02km}
    # del UH_25km,UH_02km



#%% Save or update pickle file

new_pkl = False
update_pkl = False

if new_pkl:
    dbfile = open(ip+f"environment/env_{fn:06d}.pkl", 'wb')
    env = {'time':time, 'MUcape':MUcape, 'MUcin':MUcin, 'z_LCL':z_LCL, 
            'SRH_3km':SRH_3km, 'SRH_1km':SRH_1km, 
            'shear_6km':shear_6km, 'u_shr':u_shr, 'v_shr':v_shr, 
            'UH_25km':UH_25km, 'UH_02km':UH_02km}
    # env = {'time':time}
    pickle.dump(env, dbfile)
    dbfile.close()

if update_pkl:
    if 'fn' not in locals():
        fn = input('File number:')
    dbfile = open(ip+f"environment/env_{fn:06d}.pkl", 'rb')
    env = pickle.load(dbfile)
    dbfile.close()
    
    # new_vars = {}
    new_vars = dat
    env.update(new_vars)
    dbfile = open(ip+f"environment/env_{fn:06d}.pkl", 'wb')
    pickle.dump(env, dbfile)
    dbfile.close()
    
#%% Load variables from saved pickle files

fp = '/Volumes/Promise_Pegasus_70TB/merger/merger-125m/' 
ip = '/Users/morgan.schneider/Documents/merger/merger-125m/'

fn = 14


dbfile = open(ip+f"environment/env_{fn:06d}.pkl", 'rb')
env = pickle.load(dbfile)
time = env['time']
MUcape = env['MUcape']
MUcin = env['MUcin']
z_LCL = env['z_LCL']
SRH_3km = env['SRH_3km']
SRH_1km = env['SRH_1km']
shear_6km = env['shear_6km']
u_shr = env['u_shr']
v_shr = env['v_shr']
UH_25km = env['UH_25km']
UH_02km = env['UH_02km']
dbfile.close()


ds = nc.Dataset(fp+f"cm1out_{fn:06d}.nc")
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
ix1 = 0 # ix1 = np.where(xh >= xlims[0])[0][0]
ix2 = len(xh) # ix2 = np.where(xh >= xlims[1])[0][0]
iy1 = 0 # iy1 = np.where(yh >= ylims[0])[0][0]
iy2 = len(yh) # iy2 = np.where(yh >= ylims[1])[0][0]
ix = slice(ix1,ix2)
iy = slice(iy1,iy2)
dbz = ds.variables['dbz'][:].data[0,0,iy,ix]
ds.close()


#%% Plot environmental fields from saved pickle files

figsave = False
plot_CAPE = True
plot_SRH = True
plot_UH = False
plot_SCP = False

if plot_SCP:
    SCP = mc.supercell_composite(MUcape*units('J/kg'), SRH_3km*units('m^2/s^2'), shear_6km*units('m/s'))
    SCP = SCP.magnitude
    STP = mc.significant_tornado(MUcape*units('J/kg'), z_LCL*units.m, SRH_1km*units('m^2/s^2'), shear_6km*units('m/s'))
    STP = STP.magnitude
    

xlims = [-60,-10]
ylims = [-25,25]


if plot_CAPE:
    # MUCAPE
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    plot_cfill(xh[ix], yh[iy], MUcape, 'cape', ax, datalims=[0,300], xlims=xlims, ylims=ylims, cmap='pyart_HomeyerRainbow')
    ax.contour(xh[ix], yh[iy], dbz, levels=[30], colors='k', linewidths=1, linestyles='-')
    ax.contourf(xh[ix], yh[iy], MUcin, levels=[50,300], hatches=['///'], colors=['gray'], alpha=0.4)
    # ax.contour(xh[ix], yh[iy], w_2km, levels=[-15,-10,10,15], colors='silver', linewidths=1)
    ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')
    ax.set_title(f"MUCAPE, MUCIN > 50 J/kg, 30 dBZ contour\n t={time:.0f} min", fontsize=14)
    if figsave:
        plt.savefig(ip+f"MUCAPE_{time:03.0f}min.png", dpi=400)
    
    # LCL height
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    plot_cfill(xh[ix], yh[iy], np.ma.masked_array(z_LCL, z_LCL is np.nan), 'z', ax, datalims=[500,2500], xlims=xlims, ylims=ylims, cmap='pyart_HomeyerRainbow')
    ax.contour(xh[ix], yh[iy], dbz, levels=[30], colors='k', linewidths=1, linestyles='-')
    # ax.contour(xh[ix], yh[iy], w_2km, levels=[-15,-10,10,15], colors='silver', linewidths=1)
    ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')
    ax.set_title(f"LCL height, 30 dBZ contour\n t={time:.0f} min", fontsize=14)
    if figsave:
        plt.savefig(ip+f"LCL_{time:03.0f}min.png", dpi=400)

if plot_SRH:
    # 0-3 km SRH
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    plot_cfill(xh[ix], yh[iy], env['SRH_3km'], 'srh', ax, datalims=[0,1000], xlims=xlims, ylims=ylims, cmap='pyart_HomeyerRainbow')
    ax.contour(xh[ix], yh[iy], dbz, levels=[30], colors='k', linewidths=1, linestyles='-')
    # ax.contour(xh[ix], yh[iy], w_2km, levels=[-15,-10,10,15], colors='silver', linewidths=1)
    ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')
    ax.set_title(f"0-3 km SRH, 30 dBZ contour\n t={time:.0f} min", fontsize=14)
    if figsave:
        plt.savefig(ip+f"SRH3km_{time:03.0f}min.png", dpi=400)
    
    # 0-1 km SRH
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    plot_cfill(xh[ix], yh[iy], env['SRH_1km'], 'srh', ax, datalims=[0,500], xlims=xlims, ylims=ylims, cmap='pyart_HomeyerRainbow')
    ax.contour(xh[ix], yh[iy], dbz, levels=[30], colors='k', linewidths=1, linestyles='-')
    # ax.contour(xh[ix], yh[iy], w_2km, levels=[-15,-10,10,15], colors='silver', linewidths=1)
    ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')
    ax.set_title(f"0-1 km SRH, 30 dBZ contour\n t={time:.0f} min", fontsize=14)
    if figsave:
        plt.savefig(ip+f"SRH1km_{time:03.0f}min.png", dpi=400)
    
    # # 0-6 km bulk shear
    # fig,ax = plt.subplots(1,1,figsize=(8,6))
    # plot_cfill(xh[ix], yh[iy], shear_6km, 'ws', ax, datalims=[10,50], xlims=xlims, ylims=ylims, cmap='pyart_HomeyerRainbow')
    # ax.contour(xh[ix], yh[iy], dbz, levels=[30], colors='k', linewidths=1, linestyles='-')
    # ax.quiver(xh[ix][::40], yh[iy][::40], u_shr[::40,::40], v_shr[::40,::40], scale=550, width=0.003, pivot='middle')
    # ax.quiver(xh[ix], yh[iy])
    # # ax.contour(xh[ix], yh[iy], w_2km, levels=[-15,-10,10,15], colors='silver', linewidths=1)
    # ax.set_xlabel('x (km)')
    # ax.set_ylabel('y (km)')
    # ax.set_title(f"0-6 km bulk shear vectors, 30 dBZ contour\n t={time:.0f} min", fontsize=14)
    # if figsave:
    #     plt.savefig(ip+f"shear6km_{time:03.0f}min.png", dpi=400)

if plot_UH:
    # 2-5 km UH
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    plot_cfill(xh[ix], yh[iy], np.ma.masked_array(UH_25km, UH_25km<=25), 'uh', ax, datalims=[0,500], xlims=xlims, ylims=ylims, cmap='pyart_HomeyerRainbow')
    ax.contour(xh[ix], yh[iy], dbz, levels=[30], colors='k', linewidths=1, linestyles='-')
    # ax.contour(xh[ix], yh[iy], w_2km, levels=[-15,-10,10,15], colors='silver', linewidths=1)
    ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')
    ax.set_title(f"2-5 km UH, 30 dBZ contour\n t={time:.0f} min", fontsize=14)
    if figsave:
        plt.savefig(ip+f"UH2-5km_{time:03.0f}min.png", dpi=400)
    
    # 0-2 km UH
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    plot_cfill(xh[ix], yh[iy], np.ma.masked_array(UH_02km, UH_02km<=25), 'uh', ax, datalims=[0,250], xlims=xlims, ylims=ylims, cmap='pyart_HomeyerRainbow')
    ax.contour(xh[ix], yh[iy], dbz, levels=[30], colors='k', linewidths=1, linestyles='-')
    # ax.contour(xh[ix], yh[iy], w_2km, levels=[-15,-10,10,15], colors='silver', linewidths=1)
    ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')
    ax.set_title(f"0-2 km UH, 30 dBZ contour\n t={time:.0f} min", fontsize=14)
    if figsave:
        plt.savefig(ip+f"UH0-2km_{time:03.0f}min.png", dpi=400)

if plot_SCP:
    # SCP
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    plot_cfill(xh[ix], yh[iy], np.ma.masked_array(SCP, SCP==0), 'scp', ax, datalims=[0,3], xlims=xlims, ylims=ylims, cmap='pyart_HomeyerRainbow')
    ax.contour(xh[ix], yh[iy], dbz, levels=[30], colors='k', linewidths=1, linestyles='-')
    # ax.contour(xh[ix], yh[iy], w_2km, levels=[-15,-10,10,15], colors='silver', linewidths=1)
    ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')
    ax.set_title(f"SCP, 30 dBZ contour\n t={time:.0f} min", fontsize=14)
    if figsave:
        plt.savefig(ip+f"SCP_{time:03.0f}min.png", dpi=400)
    
    # STP
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    plot_cfill(xh[ix], yh[iy], np.ma.masked_array(STP, STP==0), 'stp', ax, datalims=[0,1], xlims=xlims, ylims=ylims, cmap='pyart_HomeyerRainbow')
    ax.contour(xh[ix], yh[iy], dbz, levels=[30], colors='k', linewidths=1, linestyles='-')
    # ax.contour(xh[ix], yh[iy], w_2km, levels=[-15,-10,10,15], colors='silver', linewidths=1)
    ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')
    ax.set_title(f"STP, 30 dBZ contour\n t={time:.0f} min", fontsize=14)
    if figsave:
        plt.savefig(ip+f"STP_{time:03.0f}min.png", dpi=400)








