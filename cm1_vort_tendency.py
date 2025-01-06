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

#%%

ds = nc.Dataset(fp+"base/cm1out_000001.nc")
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
zh = ds.variables['z'][:].data
iz = np.where(zh >= 4)[0][1]
ds.close()

xlims = [-45,25]
ylims = [-120,-50]
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
pid = ds.variables['xh'][:].data
u_ml = ds.variables['u'][:].data[[np.nonzero(pid == pids_ml[i])[0][0] for i in range(len(pids_ml))]]
v_ml = ds.variables['v'][:].data[[np.nonzero(pid == pids_ml[i])[0][0] for i in range(len(pids_ml))]]
ds.close()

x_median = np.median(x_ml, axis=1)
y_median = np.median(y_ml, axis=1)
z_median = np.median(z_ml, axis=1)



# Calculate storm motion for streamwise/crosswise vorticity
# if 'u_storm' not in locals():
#     dbfile = open('/Users/morgan.schneider/Documents/merger/merger-125m/boxes_s1.pkl', 'rb') # interchange with boxes_q
#     box = pickle.load(dbfile)
#     x1_m = 1000 * np.array(box['x1_pp'])
#     y1_m = 1000 * np.array(box['y1_pp'])
#     dbfile.close()
    
#     for i in range(len(x1_m)-1):
#         if x1_m[i+1] == x1_m[i]:
#             if (i != 0):
#                 x1_m[i] = (x1_m[i+1] + x1_m[i-1]) / 2
#             else:
#                 x1_m[i+1] = (x1_m[i+2] + x1_m[i]) / 2
        
#         if y1_m[i+1] == y1_m[i]:
#             if (i != 0):
#                 y1_m[i] = (y1_m[i+1] + y1_m[i-1]) / 2
#             else:
#                 y1_m[i+1] = (y1_m[i+2] + y1_m[i]) / 2
    
#     u_s = np.zeros(shape=(61,), dtype=float); v_s = np.zeros(shape=(61,), dtype=float)
#     u_s[1:] = np.gradient(x1_m, np.linspace(10860, 14400, 60))
#     v_s[1:] = np.gradient(y1_m, np.linspace(10860, 14400, 60))
    
#     if u_s[1] == u_s[2]:
#         u_s[0] = u_s[1]
#     else:
#         u_s[0] = u_s[1] - np.diff(u_s)[1]
    
#     if v_s[1] == v_s[2]:
#         v_s[0] = v_s[1]
#     else:
#         v_s[0] = v_s[1] - np.diff(v_s)[1]
    
#     u_storm = u_s;  v_storm = v_s
#     u_storm[1:-1] = movmean(u_s,3)[1:-1]
#     v_storm[1:-1] = movmean(v_s,3)[1:-1]


times = [210, 225]
n = np.where(ptime/60 == 225)[0][0]


calc_stretching = False
calc_tilting = False
calc_baroclinic = False
calc_friction = False


for i in range(len(times)):
    t = times[i]
    it1 = np.where(ptime/60 >= t-15)[0][0]
    it2 = np.where(ptime/60 > t)[0][0]
    it = slice(it1,it2)
    
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
    dbfile.close()
    
    stimes = np.zeros(shape=(61,), dtype=float)
    xvort_term = np.zeros(shape=(61,len(pids_ml)), dtype=float)
    yvort_term = np.zeros(shape=(61,len(pids_ml)), dtype=float)
    zvort_term = np.zeros(shape=(61,len(pids_ml)), dtype=float)
    
    for k in np.arange(13,74):
        print(f"cm1out_{k:06d}")
        if k == 13:
            fp = '/Volumes/Promise_Pegasus_70TB/merger/merger-125m/base/'
        else:
            fp = '/Volumes/Promise_Pegasus_70TB/merger/merger-125m/'
        
        # u_sr = u_ml - u_storm[k-13]
        # v_sr = v_ml - v_storm[k-13]
        
        if calc_stretching | calc_tilting:
            ds = nc.Dataset(fp+f"cm1out_{k:06d}.nc")
            stime = ds.variables['time'][:].data[0]
            u = ds.variables['uinterp'][:].data[0,0:iz,:,:]
            v = ds.variables['vinterp'][:].data[0,0:iz,:,:]
            w = ds.variables['winterp'][:].data[0,0:iz,:,:]
            ds.close()
            
            if calc_stretching:
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
                    
                    xvort_term[k-13,p] = xvort_ml[it,p] * dudx_ml
                    yvort_term[k-13,p] = yvort_ml[it,p] * dvdy_ml
                    zvort_term[k-13,p] = zvort_ml[it,p] * dwdz_ml
            
            if calc_tilting:
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
                    xvort_term[k-13,p] = yvort_ml[it,p] * dudy_ml + zvort_ml[it,p] * dudz_ml
                    
                    dvdx_ml = dvdx_interp((z_ml[it,p], y_ml[it,p], x_ml[it,p]))
                    dvdz_ml = dvdz_interp((z_ml[it,p], y_ml[it,p], x_ml[it,p]))
                    yvort_term[k-13,p] = xvort_ml[it,p] * dvdx_ml + zvort_ml[it,p] * dvdz_ml
                    
                    dwdx_ml = dwdx_interp((z_ml[it,p], y_ml[it,p], x_ml[it,p]))
                    dwdy_ml = dwdy_interp((z_ml[it,p], y_ml[it,p], x_ml[it,p]))
                    zvort_term[k-13,p] = xvort_ml[it,p] * dwdx_ml + yvort_ml[it,p] * dwdy_ml
                    
            
        if calc_baroclinic:
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
                
                xvort_term[k-13,p] = 1/(1.1**2) * (drdy_ml * dpdz_ml - drdz_ml * dpdy_ml)
                yvort_term[k-13,p] = 1/(1.1**2) * (drdz_ml * dpdx_ml - drdx_ml * dpdz_ml)
                zvort_term[k-13,p] = 1/(1.1**2) * (drdx_ml * dpdy_ml - drdy_ml * dpdx_ml)
        
        
        if calc_friction:
            ds = nc.Dataset(fp+f"cm1out_{k:06d}.nc")
            stime = ds.variables['time'][:].data[0]
            xvort = ds.variables['xvort'][:].data[0,0:iz,:,:]
            yvort = ds.variables['yvort'][:].data[0,0:iz,:,:]
            zvort = ds.variables['zvort'][:].data[0,0:iz,:,:]
            ds.close()
            
            xvort_dx = np.gradient(xvort, xh*1000, axis=2)
            xvort_dy = np.gradient(xvort, yh*1000, axis=1)
            xvort_dz = np.gradient(xvort, zh[0:iz]*1000, axis=0)
            del2xvort = (np.gradient(xvort_dx, xh*1000, axis=2) + 
                         np.gradient(xvort_dy, yh*1000, axis=1) +
                         np.gradient(xvort_dz, zh[0:iz]*1000, axis=0))
            del xvort,xvort_dx,xvort_dy,xvort_dz
            
            yvort_dx = np.gradient(yvort, xh*1000, axis=2)
            yvort_dy = np.gradient(yvort, yh*1000, axis=1)
            yvort_dz = np.gradient(yvort, zh[0:iz]*1000, axis=0)
            del2yvort = (np.gradient(yvort_dx, xh*1000, axis=2) + 
                         np.gradient(yvort_dy, yh*1000, axis=1) +
                         np.gradient(yvort_dz, zh[0:iz]*1000, axis=0))
            del yvort,yvort_dx,yvort_dy,yvort_dz
            
            zvort_dx = np.gradient(zvort, xh*1000, axis=2)
            zvort_dy = np.gradient(zvort, yh*1000, axis=1)
            zvort_dz = np.gradient(zvort, zh[0:iz]*1000, axis=0)
            del2zvort = (np.gradient(zvort_dx, xh*1000, axis=2) + 
                         np.gradient(zvort_dy, yh*1000, axis=1) +
                         np.gradient(zvort_dz, zh[0:iz]*1000, axis=0))
            del zvort,zvort_dx,zvort_dy,zvort_dz
            
            del2xvort_interp = RegularGridInterpolator((zh[0:iz], yh, xh), del2xvort)
            del2yvort_interp = RegularGridInterpolator((zh[0:iz], yh, xh), del2yvort)
            del2zvort_interp = RegularGridInterpolator((zh[0:iz], yh, xh), del2zvort)
            del del2xvort,del2yvort,del2zvort
            
            it = np.where(ptime == stime)[0][0]
            
            for p in range(len(pids_ml)):
                del2xvort_ml = del2xvort_interp((z_ml[it,p], y_ml[it,p], x_ml[it,p]))
                del2yvort_ml = del2yvort_interp((z_ml[it,p], y_ml[it,p], x_ml[it,p]))
                del2zvort_ml = del2zvort_interp((z_ml[it,p], y_ml[it,p], x_ml[it,p]))
                
                # how do i calculate kinematic viscosity??
                nu = 1.46e-5
                
                xvort_term[k-13,p] = nu * del2xvort_ml
                yvort_term[k-13,p] = nu * del2yvort_ml
                zvort_term[k-13,p] = nu * del2zvort_ml
        
        
        
        
    
    
    

















































