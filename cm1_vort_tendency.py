#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 15:24:56 2025

@author: morgan.schneider
"""

####################
### Load modules ###
####################

from CM1utils import *
from scipy.interpolate import RegularGridInterpolator
from os.path import exists
import metpy.calc as mc

#%% Tendency along trajectories

ds = nc.Dataset("/Volumes/Promise_Pegasus_70TB/merger/merger-125m/base/cm1out_000001.nc")
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
zh = ds.variables['z'][:].data
iz = np.where(zh >= 4)[0][1]
ds.close()

xlims = [-25,5]
ylims = [-90,-40]
zlims = [0,4]

ix1 = np.where(xh >= xlims[0])[0][0]
ix2 = np.where(xh >= xlims[1])[0][0]
iy1 = np.where(yh >= ylims[0])[0][0]
iy2 = np.where(yh >= ylims[1])[0][0]
ix = slice(ix1,ix2+1)
iy = slice(iy1,iy2+1)
# xx,yy = np.meshgrid(xh[ix], yh[iy], indexing='xy')


dbfile = open(f"/Users/morgan.schneider/Documents/merger/traj_MV1.pkl", 'rb')
traj = pickle.load(dbfile)
dbfile.close()

ds = nc.Dataset('/Volumes/Promise_Pegasus_70TB/merger/merger-125m/cm1out_pdata.nc')
ptime = ds.variables['time'][:].data
ds.close()


# Calculate storm motion for streamwise/crosswise vorticity
if 'u_storm' not in locals():
    dbfile = open('/Users/morgan.schneider/Documents/merger/merger-125m/boxes_s1.pkl', 'rb') # interchange with boxes_q
    box = pickle.load(dbfile)
    x1_m = 1000 * np.array(box['x1_pp'])
    y1_m = 1000 * np.array(box['y1_pp'])
    dbfile.close()
    
    for i in range(len(x1_m)-1):
        if x1_m[i+1] == x1_m[i]:
            if (i != 0):
                x1_m[i] = (x1_m[i+1] + x1_m[i-1]) / 2
            else:
                x1_m[i+1] = (x1_m[i+2] + x1_m[i]) / 2
        
        if y1_m[i+1] == y1_m[i]:
            if (i != 0):
                y1_m[i] = (y1_m[i+1] + y1_m[i-1]) / 2
            else:
                y1_m[i+1] = (y1_m[i+2] + y1_m[i]) / 2
    
    u_s = np.zeros(shape=(61,), dtype=float); v_s = np.zeros(shape=(61,), dtype=float)
    u_s[1:] = np.gradient(x1_m, np.linspace(10860, 14400, 60))
    v_s[1:] = np.gradient(y1_m, np.linspace(10860, 14400, 60))
    
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


times = [210]
# n = np.where(ptime/60 == 225)[0][0]
fnums = [43]

calc_stretching = False
calc_tilting = False
calc_baroclinic = True
calc_friction = False


for i in range(len(times)):
    t = times[i]
    # it1 = np.where(ptime/60 >= t-15)[0][0]
    # it2 = np.where(ptime/60 > t)[0][0]
    # its = slice(it1,it2)
    
    dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/traj_clusters_{t:.0f}min_v2.pkl", 'rb')
    cc = pickle.load(dbfile)
    cc_mv1 = cc['mv1']
    dbfile.close()
    
    pids_ml = traj[f"{t}min"]['pids'][(cc_mv1 == 1)]
    x_ml = traj[f"{t}min"]['x'][:,(cc_mv1 == 1)]/1000
    y_ml = traj[f"{t}min"]['y'][:,(cc_mv1 == 1)]/1000
    z_ml = traj[f"{t}min"]['z'][:,(cc_mv1 == 1)]/1000
    u_ml = traj[f"{t}min"]['u'][:,(cc_mv1 == 1)]
    v_ml = traj[f"{t}min"]['v'][:,(cc_mv1 == 1)]
    w_ml = traj[f"{t}min"]['w'][:,(cc_mv1 == 1)]
    zvort_ml = traj[f"{t}min"]['zvort'][:,(cc_mv1 == 1)]
    
    dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/hvort_traj_{t:.0f}min.pkl", 'rb')
    vort_traj = pickle.load(dbfile)
    xvort_ml = vort_traj['xvort_ml']
    yvort_ml = vort_traj['yvort_ml']
    hvort_ml = vort_traj['hvort_ml']
    dbfile.close()
    
    # stimes = np.zeros(shape=(16,), dtype=float)
    xvort_term = np.zeros(shape=(16,len(pids_ml)), dtype=float)
    yvort_term = np.zeros(shape=(16,len(pids_ml)), dtype=float)
    zvort_term = np.zeros(shape=(16,len(pids_ml)), dtype=float)
    hvort_term = np.zeros(shape=(16,len(pids_ml)), dtype=float)
    svort_term = np.zeros(shape=(16,len(pids_ml)), dtype=float)
    cvort_term = np.zeros(shape=(16,len(pids_ml)), dtype=float)
    
    pkl_fn = f"/Users/morgan.schneider/Documents/merger/merger-125m/traj_vten_{t:.0f}min.pkl"
    
    for k in np.arange(fnums[i]-15,fnums[i]+1):
        m = k-13-(i+1)*15
        print(f"cm1out_{k:06d}")
        
        if k == 13:
            fp = '/Volumes/Promise_Pegasus_70TB/merger/merger-125m/base/'
        else:
            fp = '/Volumes/Promise_Pegasus_70TB/merger/merger-125m/'
        #end if k == 13
        
        us = u_ml - u_storm[k-13]
        vs = v_ml - v_storm[k-13]
        us_norm = us / np.sqrt(us**2 + vs**2)
        vs_norm = vs / np.sqrt(us**2 + vs**2)
        
        if calc_stretching | calc_tilting:
            ds = nc.Dataset(fp+f"cm1out_{k:06d}.nc")
            stime = ds.variables['time'][:].data[0]
            u = ds.variables['uinterp'][:].data[0,0:iz,:,:]
            v = ds.variables['vinterp'][:].data[0,0:iz,:,:]
            w = ds.variables['winterp'][:].data[0,0:iz,:,:]
            ds.close()
            
            if calc_stretching:
                print('...stretching...')
                
                dudx = np.gradient(u, xh*1000, axis=2)
                dvdy = np.gradient(v, yh*1000, axis=1)
                dwdz = np.gradient(w, zh[0:iz]*1000, axis=0)
                del u,v,w
                
                dudx_interp = RegularGridInterpolator((zh[0:iz], yh, xh), dudx)
                dvdy_interp = RegularGridInterpolator((zh[0:iz], yh, xh), dvdy)
                dwdz_interp = RegularGridInterpolator((zh[0:iz], yh, xh), dwdz)
                del dudx,dvdy,dwdz
                
                it = np.where(ptime == stime)[0][0]
                for p in range(len(pids_ml)):
                    dudx_ml = dudx_interp((z_ml[it,p], y_ml[it,p], x_ml[it,p]))
                    dvdy_ml = dvdy_interp((z_ml[it,p], y_ml[it,p], x_ml[it,p]))
                    dwdz_ml = dwdz_interp((z_ml[it,p], y_ml[it,p], x_ml[it,p]))
                    
                    xvort_term[m,p] = xvort_ml[it,p] * dudx_ml
                    yvort_term[m,p] = yvort_ml[it,p] * dvdy_ml
                    zvort_term[m,p] = zvort_ml[it,p] * dwdz_ml
                    hvort_term[m,p] = (xvort_ml[it,p] / hvort_ml[it,p]) * xvort_term[m,p] + (yvort_ml[it,p] / hvort_ml[it,p]) * yvort_term[m,p]
                    svort_term[m,p] = us_norm[it,p] * xvort_term[m,p] + vs_norm[it,p] * yvort_term[m,p]
                    cvort_term[m,p] = vs_norm[it,p] * xvort_term[m,p] - us_norm[it,p] * yvort_term[m,p]
                    
                    del dudx_ml,dvdy_ml,dwdz_ml
                #end for p in range(len(pids_ml))
            #end if calc_stretching
            
            if calc_tilting:
                print('...tilting...')
                
                dudy = np.gradient(u, yh*1000, axis=1)
                dudz = np.gradient(u, zh[0:iz]*1000, axis=0)
                dvdx = np.gradient(v, xh*1000, axis=2)
                dvdz = np.gradient(v, zh[0:iz]*1000, axis=0)
                dwdx = np.gradient(w, xh*1000, axis=2)
                dwdy = np.gradient(w, yh*1000, axis=1)
                del u,v,w
                
                dudy_interp = RegularGridInterpolator((zh[0:iz], yh, xh), dudy)
                dudz_interp = RegularGridInterpolator((zh[0:iz], yh, xh), dudz)
                dvdx_interp = RegularGridInterpolator((zh[0:iz], yh, xh), dvdx)
                dvdz_interp = RegularGridInterpolator((zh[0:iz], yh, xh), dvdz)
                dwdx_interp = RegularGridInterpolator((zh[0:iz], yh, xh), dwdx)
                dwdy_interp = RegularGridInterpolator((zh[0:iz], yh, xh), dwdy)
                del dudy,dudz,dvdx,dvdz,dwdx,dwdy
                
                it = np.where(ptime == stime)[0][0]
                for p in range(len(pids_ml)):
                    dudy_ml = dudy_interp((z_ml[it,p], y_ml[it,p], x_ml[it,p]))
                    dudz_ml = dudz_interp((z_ml[it,p], y_ml[it,p], x_ml[it,p]))
                    xvort_term[m,p] = yvort_ml[it,p] * dudy_ml + zvort_ml[it,p] * dudz_ml
                    
                    dvdx_ml = dvdx_interp((z_ml[it,p], y_ml[it,p], x_ml[it,p]))
                    dvdz_ml = dvdz_interp((z_ml[it,p], y_ml[it,p], x_ml[it,p]))
                    yvort_term[m,p] = xvort_ml[it,p] * dvdx_ml + zvort_ml[it,p] * dvdz_ml
                    
                    dwdx_ml = dwdx_interp((z_ml[it,p], y_ml[it,p], x_ml[it,p]))
                    dwdy_ml = dwdy_interp((z_ml[it,p], y_ml[it,p], x_ml[it,p]))
                    zvort_term[m,p] = xvort_ml[it,p] * dwdx_ml + yvort_ml[it,p] * dwdy_ml
                    
                    hvort_term[m,p] = (xvort_ml[it,p] / hvort_ml[it,p]) * xvort_term[m,p] + (yvort_ml[it,p] / hvort_ml[it,p]) * yvort_term[m,p]
                    svort_term[m,p] = us_norm[it,p] * xvort_term[m,p] + vs_norm[it,p] * yvort_term[m,p]
                    cvort_term[m,p] = vs_norm[it,p] * xvort_term[m,p] - us_norm[it,p] * yvort_term[m,p]
                    
                    del dudy_ml,dudz_ml,dvdx_ml,dvdz_ml,dwdx_ml,dwdy_ml
                #end for p in range(len(pids_ml))
            #end if calc_tilting
        #end if calc_stretching | calc_tilting
            
        if calc_baroclinic:
            print('...baroclinic...')
            
            ds = nc.Dataset(fp+f"cm1out_{k:06d}.nc")
            stime = ds.variables['time'][:].data[0]
            rho = ds.variables['rho'][:].data[0,0:iz,:,:]
            prs = ds.variables['prs'][:].data[0,0:iz,:,:]
            ds.close()
            
            drdx = np.gradient(rho, xh*1000, axis=2)
            drdy = np.gradient(rho, yh*1000, axis=1)
            drdz = np.gradient(rho, zh[0:iz]*1000, axis=0)
            dpdx = np.gradient(prs, xh*1000, axis=2)
            dpdy = np.gradient(prs, yh*1000, axis=1)
            dpdz = np.gradient(prs, zh[0:iz]*1000, axis=0)
            del rho,prs
            
            drdx_interp = RegularGridInterpolator((zh[0:iz], yh, xh), drdx)
            drdy_interp = RegularGridInterpolator((zh[0:iz], yh, xh), drdy)
            drdz_interp = RegularGridInterpolator((zh[0:iz], yh, xh), drdz)
            del drdx,drdy,drdz
            dpdx_interp = RegularGridInterpolator((zh[0:iz], yh, xh), dpdx)
            dpdy_interp = RegularGridInterpolator((zh[0:iz], yh, xh), dpdy)
            dpdz_interp = RegularGridInterpolator((zh[0:iz], yh, xh), dpdz)
            del dpdx,dpdy,dpdz
            
            it = np.where(ptime == stime)[0][0]
            for p in range(len(pids_ml)):
                drdx_ml = drdx_interp((z_ml[it,p], y_ml[it,p], x_ml[it,p]))
                drdy_ml = drdy_interp((z_ml[it,p], y_ml[it,p], x_ml[it,p]))
                drdz_ml = drdz_interp((z_ml[it,p], y_ml[it,p], x_ml[it,p]))
                dpdx_ml = dpdx_interp((z_ml[it,p], y_ml[it,p], x_ml[it,p]))
                dpdy_ml = dpdy_interp((z_ml[it,p], y_ml[it,p], x_ml[it,p]))
                dpdz_ml = dpdz_interp((z_ml[it,p], y_ml[it,p], x_ml[it,p]))
                
                xvort_term[m,p] = 1/(1.1**2) * (drdy_ml * dpdz_ml - drdz_ml * dpdy_ml)
                yvort_term[m,p] = 1/(1.1**2) * (drdz_ml * dpdx_ml - drdx_ml * dpdz_ml)
                zvort_term[m,p] = 1/(1.1**2) * (drdx_ml * dpdy_ml - drdy_ml * dpdx_ml)
                hvort_term[m,p] = (xvort_ml[it,p] / hvort_ml[it,p]) * xvort_term[m,p] + (yvort_ml[it,p] / hvort_ml[it,p]) * yvort_term[m,p]
                svort_term[m,p] = us_norm[it,p] * xvort_term[m,p] + vs_norm[it,p] * yvort_term[m,p]
                cvort_term[m,p] = vs_norm[it,p] * xvort_term[m,p] - us_norm[it,p] * yvort_term[m,p]
                
                del drdx_ml,drdy_ml,drdz_ml,dpdx_ml,dpdy_ml,dpdz_ml
            #end for p in range(len(pids_ml))
        #end if calc_baroclinic
        
        if calc_friction:
            print('...friction...')
            
            ds = nc.Dataset(fp+f"cm1out_{k:06d}.nc")
            stime = ds.variables['time'][:].data[0]
            xvort = ds.variables['xvort'][:].data[0,0:iz,:,:]
            yvort = ds.variables['yvort'][:].data[0,0:iz,:,:]
            zvort = ds.variables['zvort'][:].data[0,0:iz,:,:]
            ds.close()
            
            del2xvort = mc.laplacian(xvort, coordinates=(zh[0:iz]*1000, yh*1000, xh*1000))
            del2yvort = mc.laplacian(yvort, coordinates=(zh[0:iz]*1000, yh*1000, xh*1000))
            del2zvort = mc.laplacian(zvort, coordinates=(zh[0:iz]*1000, yh*1000, xh*1000))
            del xvort,yvort,zvort
            
            del2xvort_interp = RegularGridInterpolator((zh[0:iz], yh, xh), del2xvort)
            del2yvort_interp = RegularGridInterpolator((zh[0:iz], yh, xh), del2yvort)
            del2zvort_interp = RegularGridInterpolator((zh[0:iz], yh, xh), del2zvort)
            del del2xvort,del2yvort,del2zvort
            
            it = np.where(ptime == stime)[0][0]
            for p in range(len(pids_ml)):
                del2xvort_ml = del2xvort_interp((z_ml[it,p], y_ml[it,p], x_ml[it,p]))
                del2yvort_ml = del2yvort_interp((z_ml[it,p], y_ml[it,p], x_ml[it,p]))
                del2zvort_ml = del2zvort_interp((z_ml[it,p], y_ml[it,p], x_ml[it,p]))
                
                # kinematic viscosity??
                nu = 1.46e-5
                
                xvort_term[m,p] = nu * del2xvort_ml
                yvort_term[m,p] = nu * del2yvort_ml
                zvort_term[m,p] = nu * del2zvort_ml
                hvort_term[m,p] = (xvort_ml[it,p] / hvort_ml[it,p]) * xvort_term[m,p] + (yvort_ml[it,p] / hvort_ml[it,p]) * yvort_term[m,p]
                svort_term[m,p] = us_norm[it,p] * xvort_term[m,p] + vs_norm[it,p] * yvort_term[m,p]
                cvort_term[m,p] = vs_norm[it,p] * xvort_term[m,p] - us_norm[it,p] * yvort_term[m,p]
                
                del del2xvort_ml,del2yvort_ml,del2zvort_ml
            #end for p in range(len(pids_ml))
        #end if calc_friction
    #end for k in np.arange(fnums[i]-15,fnums[i]+1)
    
    if calc_stretching:
        data = {'stretch_x':xvort_term, 'stretch_y':yvort_term, 'stretch_z':zvort_term,
                'stretch_h':hvort_term, 'stretch_sw':svort_term, 'stretch_cw':cvort_term}
    elif calc_tilting:
        data = {'tilt_x':xvort_term, 'tilt_y':yvort_term, 'tilt_z':zvort_term,
                'tilt_h':hvort_term, 'tilt_sw':svort_term, 'tilt_cw':cvort_term}
    elif calc_baroclinic:
        data = {'bcl_x':xvort_term, 'bcl_y':yvort_term, 'bcl_z':zvort_term,
                'bcl_h':hvort_term, 'bcl_sw':svort_term, 'bcl_cw':cvort_term}
    elif calc_friction:
        data = {'fric_x':xvort_term, 'fric_y':yvort_term, 'fric_z':zvort_term,
                'fric_h':hvort_term, 'fric_sw':svort_term, 'fric_cw':cvort_term}
    #end if calc_stretching
    save_to_pickle(data, pkl_fn)
    
#end for i in range(len(times))


#%%

mvtime = 210
# times = np.arange(195,226)
t = np.linspace(mvtime-15, mvtime, 16)

dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/traj_vten_{mvtime}min.pkl", 'rb')
vten = pickle.load(dbfile)
dbfile.close()


tilt_z = vten['tilt_z']; stretch_z = vten['stretch_z']; bcl_z = vten['bcl_z']; fric_z = vten['fric_z']
tilt_h = vten['tilt_h']; stretch_h = vten['stretch_h']; bcl_h = vten['bcl_h']; fric_h = vten['fric_h']
tilt_sw = vten['tilt_sw']; stretch_sw = vten['stretch_sw']; bcl_sw = vten['bcl_sw']; fric_sw = vten['fric_sw']
tilt_cw = vten['tilt_cw']; stretch_cw = vten['stretch_cw']; bcl_cw = vten['bcl_cw']; fric_cw = vten['fric_cw']




fig,ax = plt.subplots(4, 1, figsize=(15,11), sharex=True, layout='constrained')

s1, = ax[0].plot(t, np.median(tilt_z, axis=1), 'mediumblue', linewidth=2)
s2, = ax[0].plot(t, np.median(stretch_z, axis=1), 'deepskyblue', linewidth=2)
s3, = ax[0].plot(t, np.median(bcl_z, axis=1), 'red', linewidth=2)
s4, = ax[0].plot(t, np.median(fric_z, axis=1), 'gold', linewidth=2)
ax[0].set_title("Parcel \u03B6 tendency")

ax[1].plot(t, np.median(tilt_h, axis=1), 'mediumblue', linewidth=2)
ax[1].plot(t, np.median(stretch_h, axis=1), 'deepskyblue', linewidth=2)
ax[1].plot(t, np.median(bcl_h, axis=1), 'red', linewidth=2)
ax[1].plot(t, np.median(fric_h, axis=1), 'gold', linewidth=2)
ax[1].set_title("Parcel \u03c9$_H$ tendency")

ax[2].plot(t, np.median(tilt_sw, axis=1), 'mediumblue', linewidth=2)
ax[2].plot(t, np.median(stretch_sw, axis=1), 'deepskyblue', linewidth=2)
ax[2].plot(t, np.median(bcl_sw, axis=1), 'red', linewidth=2)
ax[2].plot(t, np.median(fric_sw, axis=1), 'gold', linewidth=2)
ax[2].set_title("Parcel \u03c9$_{SW}$ tendency")

ax[3].plot(t, np.median(tilt_cw, axis=1), 'mediumblue', linewidth=2)
ax[3].plot(t, np.median(stretch_cw, axis=1), 'deepskyblue', linewidth=2)
ax[3].plot(t, np.median(bcl_cw, axis=1), 'red', linewidth=2)
ax[3].plot(t, np.median(fric_cw, axis=1), 'gold', linewidth=2)
ax[3].set_title("Parcel \u03c9$_{CW}$ tendency")

plt.show()


# fig,ax = plt.subplots(1, 1, figsize=(8,5), layout='constrained')

# ax.fill_between(times[15:], np.percentile(svort_term, 25, axis=1),
#                 np.percentile(svort_term, 75, axis=1), color='deepskyblue', alpha=0.2)
# ax.fill_between(times[15:], np.percentile(cvort_term, 25, axis=1),
#                 np.percentile(cvort_term, 75, axis=1), color='gold', alpha=0.3)
# ax.fill_between(times[15:], np.percentile(zvort_term, 25, axis=1),
#                 np.percentile(zvort_term, 75, axis=1), color='red', alpha=0.2)

# s1, = ax.plot(times[15:], np.median(svort_term, axis=1), 'deepskyblue', linewidth=2)
# s2, = ax.plot(times[15:], np.median(cvort_term, axis=1), 'gold', linewidth=2)
# s3, = ax.plot(times[15:], np.median(zvort_term, axis=1), 'red', linewidth=2)

# ax.set_xlabel('Time (min)')
# ax.set_ylabel('Friction term')
# ax.set_xlim([210,225])
# ax.legend(handles=[s1,s2,s3], labels=['sw','cw','z'], loc=1)

# plt.show()


#%% Calculate tendency plan views

# Calculate storm motion for streamwise/crosswise vorticity
if 'u_storm' not in locals():
    dbfile = open('/Users/morgan.schneider/Documents/merger/merger-125m/boxes_s1.pkl', 'rb') # interchange with boxes_q
    box = pickle.load(dbfile)
    x1_m = 1000 * np.array(box['x1_pp'])
    y1_m = 1000 * np.array(box['y1_pp'])
    dbfile.close()
    
    for i in range(len(x1_m)-1):
        if x1_m[i+1] == x1_m[i]:
            if (i != 0):
                x1_m[i] = (x1_m[i+1] + x1_m[i-1]) / 2
            else:
                x1_m[i+1] = (x1_m[i+2] + x1_m[i]) / 2
        
        if y1_m[i+1] == y1_m[i]:
            if (i != 0):
                y1_m[i] = (y1_m[i+1] + y1_m[i-1]) / 2
            else:
                y1_m[i+1] = (y1_m[i+2] + y1_m[i]) / 2
    
    u_s = np.zeros(shape=(61,), dtype=float); v_s = np.zeros(shape=(61,), dtype=float)
    u_s[1:] = np.gradient(x1_m, np.linspace(10860, 14400, 60))
    v_s[1:] = np.gradient(y1_m, np.linspace(10860, 14400, 60))
    
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



calc_stretching = False
calc_tilting = False
calc_baroclinic = False
calc_friction = False

ip = '/Users/morgan.schneider/Documents/merger/merger-125m/cross_sections/MV1_vten/'

for fn in np.arange(28,59):
    print(f"cm1out_{fn:06d}")
    
    # Read output file
    if np.any(calc_stretching, calc_tilting, calc_baroclinic, calc_friction):
        ds = nc.Dataset(f"/Volumes/Promise_Pegasus_70TB/merger/merger-125m/cm1out_{fn:06d}.nc")
        stime = ds.variables['time'][:].data[0]
        xh = ds.variables['xh'][:].data
        yh = ds.variables['yh'][:].data
        zh = ds.variables['z'][:].data
        
        iz4 = np.where(zh >= 4)[0][1]
        iz = slice(0,iz4)
        
        xlims = [-30,10] #[-55,25]
        ylims = [-90,-40]
        
        ix = slice(np.where(xh >= xlims[0])[0][0], np.where(xh >= xlims[1])[0][1])
        iy = slice(np.where(yh >= ylims[0])[0][0], np.where(yh >= ylims[1])[0][1])
        
        zz,yy,xx = np.meshgrid(zh[iz]*1000, yh[iy]*1000, xh[ix]*1000, indexing='ij')
        
        u = ds.variables['uinterp'][:].data[0,iz,iy,ix]
        v = ds.variables['vinterp'][:].data[0,iz,iy,ix]
        w = ds.variables['winterp'][:].data[0,iz,iy,ix]
        xvort = ds.variables['xvort'][:].data[0,iz,iy,ix]
        yvort = ds.variables['yvort'][:].data[0,iz,iy,ix]
        zvort = ds.variables['zvort'][:].data[0,iz,iy,ix]
        if calc_baroclinic:
            rho = ds.variables['rho'][:].data[0,iz,iy,ix]
            prs = ds.variables['prs'][:].data[0,iz,iy,ix]
        
        ds.close()
        
        u_sr = u - u_storm[fn-13]
        v_sr = v - v_storm[fn-13]
        ws_sr = np.sqrt(u_sr**2 + v_sr**2)
        hvort = np.sqrt(xvort**2 + yvort**2)
    #end if np.any(calc_stretching, calc_tilting, calc_baroclinic, calc_friction)
    
    # Stretching
    if calc_stretching:
        print('...stretching...')
        dudx = mc.gradient(u, coordinates=(zz, yy, xx))[2]
        dvdy = mc.gradient(v, coordinates=(zz, yy, xx))[1]
        dwdz = mc.gradient(w, coordinates=(zz, yy, xx))[0]
        
        stretch_x = xvort * dudx
        stretch_y = yvort * dvdy
        stretch_z = zvort * dwdz
        del dudx,dvdy,dwdz
        
        stretch_h = (xvort/hvort) * stretch_x + (yvort/hvort) * stretch_y
        stretch_sw = (u_sr/ws_sr) * stretch_x + (v_sr/ws_sr) * stretch_y
        stretch_cw = (v_sr/ws_sr) * stretch_x - (u_sr/ws_sr) * stretch_y
        
        dat = {'x':xh[ix], 'y':yh[iy], 'z':zh[iz],
               'stretch_x':stretch_x, 'stretch_y':stretch_y, 'stretch_z':stretch_z,
               'stretch_h':stretch_h, 'stretch_sw':stretch_sw, 'stretch_cw':stretch_cw}
        save_to_pickle(dat, ip+f"plan{stime/60:.0f}_stretch.pkl")
        del dat,stretch_x,stretch_y,stretch_z,stretch_h,stretch_sw,stretch_cw
    
    # Tilting
    if calc_tilting:
        print('...tilting...')
        dudy = mc.gradient(u, coordinates=(zz, yy, xx))[1]
        dudz = mc.gradient(u, coordinates=(zz, yy, xx))[0]
        tilt_x = yvort * dudy + zvort * dudz
        del dudy,dudz
        
        dvdx = mc.gradient(v, coordinates=(zz, yy, xx))[2]
        dvdz = mc.gradient(v, coordinates=(zz, yy, xx))[0]
        tilt_y = xvort * dvdx + zvort * dvdz
        del dvdx,dvdz
        
        dwdx = mc.gradient(w, coordinates=(zz, yy, xx))[2]
        dwdy = mc.gradient(w, coordinates=(zz, yy, xx))[1]
        tilt_z = xvort * dwdx + yvort * dwdy
        del dwdx,dwdy
        
        tilt_h = (xvort/hvort) * tilt_x + (yvort/hvort) * tilt_y
        tilt_sw = (u_sr/ws_sr) * tilt_x + (v_sr/ws_sr) * tilt_y
        tilt_cw = (v_sr/ws_sr) * tilt_x - (u_sr/ws_sr) * tilt_y
        
        dat = {'x':xh[ix], 'y':yh[iy], 'z':zh[iz],
               'tilt_x':tilt_x, 'tilt_y':tilt_y, 'tilt_z':tilt_z,
               'tilt_h':tilt_h, 'tilt_sw':tilt_sw, 'tilt_cw':tilt_cw}
        save_to_pickle(dat, ip+f"plan{stime/60:.0f}_tilt.pkl")
        del dat,tilt_x,tilt_y,tilt_z,tilt_h,tilt_sw,tilt_cw
    
    # Baroclinic
    if calc_baroclinic:
        print('...baroclinic...')
        drdx = mc.gradient(rho, coordinates=(zz, yy, xx))[2]
        drdy = mc.gradient(rho, coordinates=(zz, yy, xx))[1]
        drdz = mc.gradient(rho, coordinates=(zz, yy, xx))[0]
        dpdx = mc.gradient(prs, coordinates=(zz, yy, xx))[2]
        dpdy = mc.gradient(prs, coordinates=(zz, yy, xx))[1]
        dpdz = mc.gradient(prs, coordinates=(zz, yy, xx))[0]
        del rho,prs
        
        bcl_x = (1/1.1)**2 * (drdy * dpdz - drdz * dpdy)
        bcl_y = (1/1.1)**2 * (drdz * dpdx - drdx * dpdz)
        bcl_z = (1/1.1)**2 * (drdx * dpdy - drdy * dpdx)
        del drdx,drdy,drdz,dpdx,dpdy,dpdz
        
        bcl_h = (xvort/hvort) * bcl_x + (yvort/hvort) * bcl_y
        bcl_sw = (u_sr/ws_sr) * bcl_x + (v_sr/ws_sr) * bcl_y
        bcl_cw = (v_sr/ws_sr) * bcl_x - (u_sr/ws_sr) * bcl_y
        
        dat = {'x':xh[ix], 'y':yh[iy], 'z':zh[iz],
               'bcl_x':bcl_x, 'bcl_y':bcl_y, 'bcl_z':bcl_z,
               'bcl_h':bcl_h, 'bcl_sw':bcl_sw, 'bcl_cw':bcl_cw}
        save_to_pickle(dat, ip+f"plan{stime/60:.0f}_baroclinic.pkl")
        del dat,bcl_x,bcl_y,bcl_z,bcl_h,bcl_sw,bcl_cw
    
    # Friction
    if calc_friction:
        print('...friction...')
        nu = 1.46e-5
        
        del2xvort = mc.laplacian(xvort, coordinates=(zz, yy, xx))
        del2yvort = mc.laplacian(yvort, coordinates=(zz, yy, xx))
        del2zvort = mc.laplacian(zvort, coordinates=(zz, yy, xx))
        
        fric_x = nu * del2xvort
        fric_y = nu * del2yvort
        fric_z = nu * del2zvort
        del del2xvort,del2yvort,del2zvort
        
        fric_h = (xvort/hvort) * fric_x + (yvort/hvort) * fric_y
        fric_sw = (u_sr/ws_sr) * fric_x + (v_sr/ws_sr) * fric_y
        fric_cw = (v_sr/ws_sr) * fric_x - (u_sr/ws_sr) * fric_y
        
        dat = {'x':xh[ix], 'y':yh[iy], 'z':zh[iz],
               'fric_x':fric_x, 'fric_y':fric_y, 'fric_z':fric_z,
               'fric_h':fric_h, 'fric_sw':fric_sw, 'fric_cw':fric_cw}
        save_to_pickle(dat, ip+f"plan{stime/60:.0f}_friction.pkl")
        del dat,fric_x,fric_y,fric_z,fric_h,fric_sw,fric_cw
    
    del u,v,w,xvort,yvort,zvort,u_sr,v_sr,ws_sr,hvort


#%% Plot plan views of tendency

ip = '/Users/morgan.schneider/Documents/merger/merger-125m/cross_sections/MV1_vten/'

t = 210
mvtime = 210

dbfile = open(ip+f"plan{t}_tilt.pkl", 'rb')
tmp = pickle.load(dbfile)
xh = tmp['x']
yh = tmp['y']
zh = tmp['z']
tilt_h = tmp['tilt_h']
tilt_z = tmp['tilt_z']
dbfile.close()

dbfile = open(ip+f"plan{t}_stretch.pkl", 'rb')
tmp = pickle.load(dbfile)
stretch_h = tmp['stretch_h']
stretch_z = tmp['stretch_z']
dbfile.close()

dbfile = open(ip+f"plan{t}_baroclinic.pkl", 'rb')
tmp = pickle.load(dbfile)
bcl_h = tmp['bcl_h']
bcl_z = tmp['bcl_z']
bcl_sw = tmp['bcl_sw']
bcl_cw = tmp['bcl_cw']
dbfile.close()

dbfile = open(ip+f"plan{t}_friction.pkl", 'rb')
tmp = pickle.load(dbfile)
fric_h = tmp['fric_h']
fric_z = tmp['fric_z']
fric_sw = tmp['fric_sw']
fric_cw = tmp['fric_cw']
dbfile.close()


# Read parcel data
ds = nc.Dataset(fp+'cm1out_pdata.nc')
ptime = ds.variables['time'][:].data
# xp = ds.variables['x'][:].data
# yp = ds.variables['y'][:].data
# zp = ds.variables['z'][:].data
ds.close()

dbfile = open('/Users/morgan.schneider/Documents/merger/traj_MV1.pkl', 'rb')
traj = pickle.load(dbfile)
x_mv = traj[f"{mvtime}min"]['x']/1000
y_mv = traj[f"{mvtime}min"]['y']/1000
z_mv = traj[f"{mvtime}min"]['z']/1000
# w_mv = traj[f"{mvtime}min"]['w']
# zvort_mv = traj[f"{mvtime}min"]['zvort']
dbfile.close()

dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/traj_clusters_{mvtime}min_v2.pkl", 'rb')
tmp = pickle.load(dbfile)
cc = tmp['mv1']
dbfile.close()

# Mid-level source
x_ml = x_mv[:,(cc==1)]
y_ml = y_mv[:,(cc==1)]
z_ml = z_mv[:,(cc==1)]

it = np.where(ptime/60 == t)[0][0]
iz = np.where(zh >= np.median(z_ml[it,:]))[0][0]
iz1 = np.where(zh >= 1)[0][0]

ix = np.where(xh >= np.median(x_ml[it,:]))[0][0]
iy = np.where(yh >= np.median(y_ml[it,:]))[0][0]


### zvort lims (208 min)
# tilting max/min: 0.0023, -0.0024
# stretching max/min: 0.0026, -0.0031
# baroclinic max/min: 0.00001, -0.00002
# friction max/min: tiny, tiny

### hvort lims (223 min)
# tilting max/min: 0.0035, -0.0029
# stretching max/min: 0.0032, -0.0021
# baroclinic max/min: 0.00089, -0.00057
# friction max/min: tiny, tiny

print(f"...zvort limits, {t} min, {zh[iz]*1000:.0f} m...")
print(f"Tilting: {np.min(tilt_z[iz,:,:]):.6f}, {np.max(tilt_z[iz,:,:]):.6f}")
print(f"Stretching: {np.min(stretch_z[iz,:,:]):.6f}, {np.max(stretch_z[iz,:,:]):.6f}")
print(f"Baroclinic: {np.min(bcl_z[iz,:,:])}, {np.max(bcl_z[iz,:,:])}")
print(f"Friction: {np.min(fric_z[iz,:,:])}, {np.max(fric_z[iz,:,:])}\n")

print(f"...hvort limits, {t} min, {zh[iz]*1000:.0f} m...")
print(f"Tilting: {np.min(tilt_h[iz,:,:]):.6f}, {np.max(tilt_h[iz,:,:]):.6f}")
print(f"Stretching: {np.min(stretch_h[iz,:,:]):.6f}, {np.max(stretch_h[iz,:,:]):.6f}")
print(f"Baroclinic: {np.min(bcl_h[iz,:,:])}, {np.max(bcl_h[iz,:,:])}")
print(f"Friction: {np.min(fric_h[iz,:,:])}, {np.max(fric_h[iz,:,:])}")


#%%

# def find_power(x):
#     if x < 0:
#         logx = np.log10(np.abs(x))
#     else:
#         logx = np.log10(x)
#     ex = np.floor(logx)
#     c = x / 10**ex
    
#     return c,ex


xlims = [-20,0]
ylims = [-60,-40]

sz = 100

# 220 min: [-0.0020, 0.0026], [-0.0017, 0.0019], [-4.56e-6, 7.97e-6], [-2.07e-10, 2.23e-10]
# 221 min: [-0.0035, 0.0048], [-0.0030, 0.0019], [-1.27e-5, 1.99e-5], [-3.01e-10, 3.14e-10]
# 222 min: [-0.0029, 0.0046], [-0.0038, 0.0020], [-1.31e-5, 1.06e-5], [-4.21e-10, 3.99e-10]
# 223 min: [-0.0023, 0.0027], [-0.0023, 0.0020], [-9.02e-6, 1.03e-5], [-3.03e-10, 2.77e-10]
# 224 min: [-0.0035, 0.0063], [-0.0030, 0.0033], [-1.74e-5, 1.53e-5], [-5.34e-10, 7.52e-10]
# 225 min: [-0.0034, 0.0022], [-0.0023, 0.0020], [-1.45e-5, 1.38e-5], [-3.47e-10, 2.49e-10]

tlims = [-0.004, 0.004]; tlevs = np.linspace(tlims[0], tlims[1], 21)
slims = [-0.003, 0.003]; slevs = np.linspace(slims[0], slims[1], 21)
blims = [-2e-5, 2e-5]; blevs = np.linspace(blims[0], blims[1], 21)
flims = [-3e-10, 3e-10]; flevs = np.linspace(flims[0], flims[1], 21)


fig,ax = plt.subplots(2,2, figsize=(11,8), sharex=True, sharey=True, layout='constrained', subplot_kw=dict(box_aspect=1))

plot_contourf(xh, yh, tilt_z[iz,:,:], 'zvort', ax[0,0], levels=tlevs, datalims=tlims, xlims=xlims, ylims=ylims)
ax[0,0].set_title(f"Tilting")
ax[0,0].scatter(x_ml[it,:], y_ml[it,:], s=sz, color='k', marker='.')

plot_contourf(xh, yh, stretch_z[iz,:,:], 'zvort', ax[0,1], levels=slevs, datalims=slims, xlims=xlims, ylims=ylims)
ax[0,1].set_title(f"Stretching")
ax[0,1].scatter(x_ml[it,:], y_ml[it,:], s=sz, color='k', marker='.')

plot_contourf(xh, yh, bcl_z[iz,:,:], 'zvort', ax[1,0], levels=blevs, datalims=blims, xlims=xlims, ylims=ylims)
ax[1,0].set_title(f"Baroclinic")
ax[1,0].scatter(x_ml[it,:], y_ml[it,:], s=sz, color='k', marker='.')

plot_contourf(xh, yh, fric_z[iz,:,:], 'zvort', ax[1,1], levels=flevs, datalims=flims, xlims=xlims, ylims=ylims)
ax[1,1].set_title(f"Friction")
ax[1,1].scatter(x_ml[it,:], y_ml[it,:], s=sz, color='k', marker='.')

plt.suptitle(f"\u03B6 tendency - {t} min, {zh[iz]*1000:.0f} m")


# 220 min: [-0.0023, 0.0037], [-0.0023, 0.0024], [-4.60e-4, 5.81e-4], [-4.92e-10, 2.29e-10]
# 221 min: [-0.0019, 0.0030], [-0.0027, 0.0062], [-7.38e-4, 1.00e-3], [-1.04e-9, 4.89e-10]
# 222 min: [-0.0040, 0.0039], [-0.0022, 0.0055], [-6.17e-4, 1.02e-3], [-1.44e-9, 8.17e-10]
# 223 min: [-0.0029, 0.0035], [-0.0021, 0.0032], [-5.71e-4, 8.95e-4], [-7.56e-10, 4.57e-10]
# 224 min: [-0.0039, 0.0041], [-0.0027, 0.0052], [-9.75e-4, 1.03e-3], [-2.83e-9, 2.03e-9]
# 225 min: [-0.0030, 0.0034], [-0.0027, 0.0044], [-5.80e-4, 8.77e-4], [-7.13e-10, 3.98e-10]

tlims = [-0.003, 0.003]; tlevs = np.linspace(tlims[0], tlims[1], 21)
slims = [-0.006, 0.006]; slevs = np.linspace(slims[0], slims[1], 21)
blims = [-10e-4, 10e-4]; blevs = np.linspace(blims[0], blims[1], 21)
flims = [-10e-10, 10e-10]; flevs = np.linspace(flims[0], flims[1], 21)


fig,ax = plt.subplots(2,2, figsize=(11,8), sharex=True, sharey=True, layout='constrained', subplot_kw=dict(box_aspect=1))

plot_contourf(xh, yh, tilt_h[iz,:,:], 'hvort', ax[0,0], levels=tlevs, datalims=tlims, xlims=xlims, ylims=ylims, cmap='pyart_balance')
ax[0,0].set_title(f"Tilting")
ax[0,0].scatter(x_ml[it,:], y_ml[it,:], s=sz, color='k', marker='.')

plot_contourf(xh, yh, stretch_h[iz,:,:], 'hvort', ax[0,1], levels=slevs, datalims=slims, xlims=xlims, ylims=ylims, cmap='pyart_balance')
ax[0,1].set_title(f"Stretching")
ax[0,1].scatter(x_ml[it,:], y_ml[it,:], s=sz, color='k', marker='.')

plot_contourf(xh, yh, bcl_h[iz,:,:], 'hvort', ax[1,0], levels=blevs, datalims=blims, xlims=xlims, ylims=ylims, cmap='pyart_balance')
ax[1,0].set_title(f"Baroclinic")
ax[1,0].scatter(x_ml[it,:], y_ml[it,:], s=sz, color='k', marker='.')

plot_contourf(xh, yh, fric_h[iz,:,:], 'hvort', ax[1,1], levels=flevs, datalims=flims, xlims=xlims, ylims=ylims, cmap='pyart_balance')
ax[1,1].set_title(f"Friction")
ax[1,1].scatter(x_ml[it,:], y_ml[it,:], s=sz, color='k', marker='.')

plt.suptitle(f"\u03c9$_H$ tendency - {t} min, {zh[iz]*1000:.0f} m")

plt.show()






























