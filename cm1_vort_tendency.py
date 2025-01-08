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


times = [210, 225]
n = np.where(ptime/60 == 225)[0][0]
fnums = [43, 58]

calc_stretching = False
calc_tilting = False
calc_baroclinic = False
calc_friction = True


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
                xvarname = 'stretch_x'
                yvarname = 'stretch_y'
                zvarname = 'stretch_z'
                hvarname = 'stretch_h'
                svarname = 'stretch_sw'
                cvarname = 'stretch_cw'
                
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
                #end for p in range(len(pids_ml))
            #end if calc_stretching
            
            if calc_tilting:
                xvarname = 'tilt_x'
                yvarname = 'tilt_y'
                zvarname = 'tilt_z'
                hvarname = 'tilt_h'
                svarname = 'tilt_sw'
                cvarname = 'tilt_cw'
                
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
                #end for p in range(len(pids_ml))
            #end if calc_tilting
        #end if calc_stretching | calc_tilting
            
        if calc_baroclinic:
            xvarname = 'baroclinic_x'
            yvarname = 'baroclinic_y'
            zvarname = 'baroclinic_z'
            hvarname = 'baroclinic_h'
            svarname = 'baroclinic_sw'
            cvarname = 'baroclinic_cw'
            
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
            #end for p in range(len(pids_ml))
        #end if calc_baroclinic
        
        if calc_friction:
            xvarname = 'fric_x'
            yvarname = 'fric_y'
            zvarname = 'fric_z'
            hvarname = 'fric_h'
            svarname = 'fric_sw'
            cvarname = 'fric_cw'
            
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
            #end for p in range(len(pids_ml))
        #end if calc_friction
    #end for k in np.arange(fnums[i]-15,fnums[i]+1)
    
    if ~exists(f"/Users/morgan.schneider/Documents/merger/merger-125m/traj_vten_{t:.0f}min.pkl"):
        dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/traj_vten_{t:.0f}min.pkl", 'wb')
        ten = {xvarname:xvort_term, yvarname:yvort_term, zvarname:zvort_term,
               hvarname:hvort_term, svarname:svort_term, cvarname:cvort_term}
        pickle.dump(ten, dbfile)
        dbfile.close()
    else:
        dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/traj_vten_{t:.0f}min.pkl", 'rb')
        ten = pickle.load(dbfile)
        dbfile.close()
        
        ten.update({xvarname:xvort_term, yvarname:yvort_term, zvarname:zvort_term,
                    hvarname:hvort_term, svarname:svort_term, cvarname:cvort_term})
        dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/traj_vten_{t:.0f}min.pkl", 'wb')
        pickle.dump(ten, dbfile)
        dbfile.close()
    #end if ~exists
#end for i in range(len(times))


#%%

times = np.arange(195,226)

fig,ax = plt.subplots(1, 1, figsize=(8,5), layout='constrained')

ax.fill_between(times[15:], np.percentile(svort_term, 25, axis=1),
                np.percentile(svort_term, 75, axis=1), color='deepskyblue', alpha=0.2)
ax.fill_between(times[15:], np.percentile(cvort_term, 25, axis=1),
                np.percentile(cvort_term, 75, axis=1), color='gold', alpha=0.3)
ax.fill_between(times[15:], np.percentile(zvort_term, 25, axis=1),
                np.percentile(zvort_term, 75, axis=1), color='red', alpha=0.2)

s1, = ax.plot(times[15:], np.median(svort_term, axis=1), 'deepskyblue', linewidth=2)
s2, = ax.plot(times[15:], np.median(cvort_term, axis=1), 'gold', linewidth=2)
s3, = ax.plot(times[15:], np.median(zvort_term, axis=1), 'red', linewidth=2)

ax.set_xlabel('Time (min)')
ax.set_ylabel('Friction term')
ax.set_xlim([210,225])
ax.legend(handles=[s1,s2,s3], labels=['sw','cw','z'], loc=1)

plt.show()


#%% Tendency plan views

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



calc_stretching = True
calc_tilting = True
calc_baroclinic = True
calc_friction = True

ip = '/Users/morgan.schneider/Documents/merger/merger-125m/cross_sections/MV1_vten/'

for fn in np.arange(28,59):
    print(f"cm1out_{fn:06d}")
    # Read output file
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












































