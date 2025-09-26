#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 15:24:56 2025

@author: morgan.schneider

- PAPER FIGURES -
composite_zvort (hvort, xvort, yvort, swvort, cwvort)

- Other figures -
Tendency plan views at single times
Parcel time series of tendency terms
"""

####################
### Load modules ###
####################

from CM1utils import *
from scipy.interpolate import RegularGridInterpolator

#%% Calculate tendency along trajectories

ds = nc.Dataset("/Volumes/Promise_Pegasus_70TB/merger/merger-125m/base/cm1out_000001.nc")
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
zh = ds.variables['z'][:].data
iz = np.where(zh >= 4.1)[0][1]
ds.close()

xlims = [-35,5]
ylims = [-135,-45]
zlims = [0,4.1]

ix1 = np.where(xh >= xlims[0])[0][0]
ix2 = np.where(xh >= xlims[1])[0][1]
iy1 = np.where(yh >= ylims[0])[0][0]
iy2 = np.where(yh >= ylims[1])[0][1]
ix = slice(ix1,ix2)
iy = slice(iy1,iy2)
# xx,yy = np.meshgrid(xh[ix], yh[iy], indexing='xy')


dbfile = open(f"/Users/morgan.schneider/Documents/merger/traj_MV1.pkl", 'rb')
traj = pickle.load(dbfile)
dbfile.close()

ds = nc.Dataset('/Volumes/Promise_Pegasus_70TB/merger/merger-125m/cm1out_pdata.nc')
ptime = ds.variables['time'][:].data
ds.close()

dbfile = open('/Users/morgan.schneider/Documents/merger/merger-125m/storm_motion.pkl', 'rb')
sm = pickle.load(dbfile)
u_storm = sm['u_storm']
v_storm = sm['v_storm']
dbfile.close()


times = [220]
# n = np.where(ptime/60 == 225)[0][0]
fnums = [53]

calc_stretching = True
calc_tilting = False
calc_baroclinic = False
calc_friction = False

for i in range(len(times)):
    t = times[i]
    # it1 = np.where(ptime/60 >= t-15)[0][0]
    # it2 = np.where(ptime/60 > t)[0][0]
    # its = slice(it1,it2)
    
    dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/traj_clusters/traj_clusters_{t:.0f}min_v2.pkl", 'rb')
    tmp = pickle.load(dbfile)
    cc = tmp['mv1']
    dbfile.close()
    
    pids_ml = traj[f"{t}min"]['pids'][(cc == 1)]
    x_ml = traj[f"{t}min"]['x'][:,(cc == 1)]/1000
    y_ml = traj[f"{t}min"]['y'][:,(cc == 1)]/1000
    z_ml = traj[f"{t}min"]['z'][:,(cc == 1)]/1000
    u_ml = traj[f"{t}min"]['u'][:,(cc == 1)]
    v_ml = traj[f"{t}min"]['v'][:,(cc == 1)]
    w_ml = traj[f"{t}min"]['w'][:,(cc == 1)]
    zvort_ml = traj[f"{t}min"]['zvort'][:,(cc == 1)]
    
    dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/hvort_traj_{t:.0f}min.pkl", 'rb')
    vort_traj = pickle.load(dbfile)
    xvort_ml = vort_traj['xvort_ml']
    yvort_ml = vort_traj['yvort_ml']
    hvort_ml = vort_traj['hvort_ml']
    svort_ml = vort_traj['vort_sw_ml']
    cvort_ml = vort_traj['vort_cw_ml']
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
        m = k-13-(i+2)*15
        print(f"cm1out_{k:06d}")
        
        if k == 13:
            fp = '/Volumes/Promise_Pegasus_70TB/merger/merger-125m/base/'
        else:
            fp = '/Volumes/Promise_Pegasus_70TB/merger/merger-125m/'
        #end if k == 13
        
        u_sr = u_ml - u_storm[k-13]
        v_sr = v_ml - v_storm[k-13]
        ws_sr = np.sqrt(u_sr**2 + v_sr**2)
        
        dudt = np.gradient(u_sr, ptime, axis=0)
        dvdt = np.gradient(v_sr, ptime, axis=0)
        dwsdt = 1/ws_sr * (u_sr * dudt + v_sr * dvdt)
        
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
                    
                    svort_term[m,p] = (1/ws_sr[it,p] * (u_sr[it,p] * xvort_term[m,p] + v_sr[it,p] * yvort_term[m,p]
                                                        + xvort_ml[it,p] * dudt[it,p] + yvort_ml[it,p] * dvdt[it,p]) 
                                       - ws_sr[it,p]**-2 * dwsdt[it,p] * (u_sr[it,p] * xvort_ml[it,p] + v_sr[it,p] * yvort_ml[it,p]))
                    cvort_term[m,p] = (1/ws_sr[it,p] * (-v_sr[it,p] * xvort_term[m,p] + u_sr[it,p] * yvort_term[m,p]
                                                        - xvort_ml[it,p] * dvdt[it,p] + yvort_ml[it,p] * dudt[it,p])
                                       - ws_sr[it,p]**-2 * dwsdt[it,p] * (-v_sr[it,p] * xvort_ml[it,p] + u_sr[it,p] * yvort_ml[it,p]))
                    
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
                    svort_term[m,p] = (1/ws_sr[it,p] * (u_sr[it,p] * xvort_term[m,p] + v_sr[it,p] * yvort_term[m,p]
                                                        + xvort_ml[it,p] * dudt[it,p] + yvort_ml[it,p] * dvdt[it,p]) 
                                       - ws_sr[it,p]**-2 * dwsdt[it,p] * (u_sr[it,p] * xvort_ml[it,p] + v_sr[it,p] * yvort_ml[it,p]))
                    cvort_term[m,p] = (1/ws_sr[it,p] * (-v_sr[it,p] * xvort_term[m,p] + u_sr[it,p] * yvort_term[m,p]
                                                        - xvort_ml[it,p] * dvdt[it,p] + yvort_ml[it,p] * dudt[it,p])
                                       - ws_sr[it,p]**-2 * dwsdt[it,p] * (-v_sr[it,p] * xvort_ml[it,p] + u_sr[it,p] * yvort_ml[it,p]))
                    
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
                svort_term[m,p] = (1/ws_sr[it,p] * (u_sr[it,p] * xvort_term[m,p] + v_sr[it,p] * yvort_term[m,p]
                                                    + xvort_ml[it,p] * dudt[it,p] + yvort_ml[it,p] * dvdt[it,p]) 
                                   - ws_sr[it,p]**-2 * dwsdt[it,p] * (u_sr[it,p] * xvort_ml[it,p] + v_sr[it,p] * yvort_ml[it,p]))
                cvort_term[m,p] = (1/ws_sr[it,p] * (-v_sr[it,p] * xvort_term[m,p] + u_sr[it,p] * yvort_term[m,p]
                                                    - xvort_ml[it,p] * dvdt[it,p] + yvort_ml[it,p] * dudt[it,p])
                                   - ws_sr[it,p]**-2 * dwsdt[it,p] * (-v_sr[it,p] * xvort_ml[it,p] + u_sr[it,p] * yvort_ml[it,p]))
                
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
                svort_term[m,p] = (1/ws_sr[it,p] * (u_sr[it,p] * xvort_term[m,p] + v_sr[it,p] * yvort_term[m,p]
                                                    + xvort_ml[it,p] * dudt[it,p] + yvort_ml[it,p] * dvdt[it,p]) 
                                   - ws_sr[it,p]**-2 * dwsdt[it,p] * (u_sr[it,p] * xvort_ml[it,p] + v_sr[it,p] * yvort_ml[it,p]))
                cvort_term[m,p] = (1/ws_sr[it,p] * (-v_sr[it,p] * xvort_term[m,p] + u_sr[it,p] * yvort_term[m,p]
                                                    - xvort_ml[it,p] * dvdt[it,p] + yvort_ml[it,p] * dudt[it,p])
                                   - ws_sr[it,p]**-2 * dwsdt[it,p] * (-v_sr[it,p] * xvort_ml[it,p] + u_sr[it,p] * yvort_ml[it,p]))
                
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


#%% Time series - single time, panels z/h/sw/cw, all terms

mvtime = 225
# times = np.arange(195,226)
t = np.linspace(mvtime-15, mvtime, 16)

dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/traj_vten_{mvtime}min.pkl", 'rb')
vten = pickle.load(dbfile)
dbfile.close()


tilt_z = vten['tilt_z']; stretch_z = vten['stretch_z']; bcl_z = vten['bcl_z']; fric_z = vten['fric_z']; total_z = vten['total_z']
tilt_h = vten['tilt_h']; stretch_h = vten['stretch_h']; bcl_h = vten['bcl_h']; fric_h = vten['fric_h']; total_h = vten['total_h']
tilt_sw = vten['tilt_sw']; stretch_sw = vten['stretch_sw']; bcl_sw = vten['bcl_sw']; fric_sw = vten['fric_sw']; total_sw = vten['total_sw']
tilt_cw = vten['tilt_cw']; stretch_cw = vten['stretch_cw']; bcl_cw = vten['bcl_cw']; fric_cw = vten['fric_cw']; total_cw = vten['total_cw']

if False:
    tp = t
    fric_z = total_z - tilt_z - stretch_z - bcl_z
    fric_h = total_h - tilt_h - stretch_h - bcl_h
    fric_sw = total_sw - tilt_sw - stretch_sw - bcl_sw
    fric_cw = total_cw - tilt_cw - stretch_cw - bcl_cw


if True:
    ds = nc.Dataset("/Volumes/Promise_Pegasus_70TB/merger/merger-125m/cm1out_pdata.nc")
    ptime = ds.variables['time'][:].data
    ds.close()
    
    dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/traj_clusters/traj_clusters_{mvtime:.0f}min_v2.pkl", 'rb')
    cc = pickle.load(dbfile)
    cc_mv1 = cc['mv1']
    dbfile.close()
    
    dbfile = open(f"/Users/morgan.schneider/Documents/merger/traj_MV1.pkl", 'rb')
    traj = pickle.load(dbfile)
    dbfile.close()
    pids_ml = traj[f"{mvtime}min"]['pids'][(cc_mv1 == 1)]
    zvort_ml = traj[f"{mvtime}min"]['zvort'][:,(cc_mv1 == 1)]
    
    dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/hvort_traj_{mvtime:.0f}min.pkl", 'rb')
    vort_traj = pickle.load(dbfile)
    xvort_ml = vort_traj['xvort_ml']
    yvort_ml = vort_traj['yvort_ml']
    hvort_ml = vort_traj['hvort_ml']
    svort_ml = vort_traj['vort_sw_ml']
    cvort_ml = vort_traj['vort_cw_ml_signed']
    dbfile.close()
    
    it1 = np.where(ptime/60 == mvtime-15)[0][0]
    it2 = np.where(ptime/60 == mvtime)[0][0]+1
    tp = ptime[it1:it2]/60
    
    total_z = np.gradient(zvort_ml[it1:it2,:], ptime[it1:it2], axis=0)
    total_h = np.gradient(hvort_ml[it1:it2,:], ptime[it1:it2], axis=0)
    total_sw = np.gradient(svort_ml[it1:it2,:], ptime[it1:it2], axis=0)
    total_cw = np.gradient(cvort_ml[it1:it2,:], ptime[it1:it2], axis=0)
    
    
    fric_z = total_z[::4,:] - tilt_z - stretch_z - bcl_z
    fric_h = total_h[::4,:] - tilt_h - stretch_h - bcl_h
    fric_sw = total_sw[::4,:] - tilt_sw - stretch_sw - bcl_sw
    fric_cw = total_cw[::4,:] - tilt_cw - stretch_cw - bcl_cw
    
    tp = t
    total_z = total_z[::4,:]
    total_h = total_h[::4,:]
    total_sw = total_sw[::4,:]
    total_cw = total_cw[::4,:]



fig,ax = plt.subplots(4, 1, figsize=(8,11), sharex=True, layout='constrained')

ax[0].axhline(0, color='gray', linewidth=1, linestyle='--')
# ax.fill_between(t, np.percentile(svort_term, 25, axis=1),
#                 np.percentile(svort_term, 75, axis=1), color='deepskyblue', alpha=0.2)
s1, = ax[0].plot(t, np.median(tilt_z, axis=1), 'mediumblue', linewidth=2)
s2, = ax[0].plot(t, np.median(stretch_z, axis=1), 'deepskyblue', linewidth=2)
s3, = ax[0].plot(t, np.median(bcl_z, axis=1), 'red', linewidth=2)
s4, = ax[0].plot(t, np.median(fric_z, axis=1), 'goldenrod', linewidth=2)
s5, = ax[0].plot(tp, np.median(total_z, axis=1), 'k', linewidth=2)
ax[0].set_title("Parcel \u03B6 tendency")
# ax[0].set_ylim([-0.0001, 0.0003])
ax[0].set_ylim([-0.0008, 0.0008])
ax[0].yaxis.set_major_locator(MultipleLocator(0.0002))
ax[0].legend(handles=[s1,s2,s3,s4,s5], labels=['Tilting','Stretching','Baroclinic','Friction','Total'], loc=2, ncol=2)
ax[0].grid(visible=True, which='major', color='darkgray', linestyle='-')
ax[0].grid(visible=True, which='minor', color='lightgray', linestyle='--')

ax[1].axhline(0, color='gray', linewidth=1, linestyle='--')
ax[1].plot(t, np.median(tilt_h, axis=1), 'mediumblue', linewidth=2)
ax[1].plot(t, np.median(stretch_h, axis=1), 'deepskyblue', linewidth=2)
ax[1].plot(t, np.median(bcl_h, axis=1), 'red', linewidth=2)
ax[1].plot(t, np.median(fric_h, axis=1), 'goldenrod', linewidth=2)
ax[1].plot(tp, np.median(total_h, axis=1), 'k', linewidth=2)
ax[1].set_title("Parcel \u03c9$_H$ tendency")
# ax[1].set_ylim([-0.0002, 0.0002])
ax[1].set_ylim([-0.0008, 0.0008])
ax[1].yaxis.set_major_locator(MultipleLocator(0.0002))
ax[1].grid(visible=True, which='major', color='darkgray', linestyle='-')
ax[1].grid(visible=True, which='minor', color='lightgray', linestyle='--')

ax[2].axhline(0, color='gray', linewidth=1, linestyle='--')
ax[2].plot(t, np.median(tilt_sw, axis=1), 'mediumblue', linewidth=2)
ax[2].plot(t, np.median(stretch_sw, axis=1), 'deepskyblue', linewidth=2)
ax[2].plot(t, np.median(bcl_sw, axis=1), 'red', linewidth=2)
ax[2].plot(t, np.median(fric_sw, axis=1), 'goldenrod', linewidth=2)
ax[2].plot(tp, np.median(total_sw, axis=1), 'k', linewidth=2)
ax[2].set_title("Parcel \u03c9$_{SW}$ tendency")
# ax[2].set_ylim([-0.0004, 0.0003])
ax[2].set_ylim([-0.0008, 0.0008])
ax[2].yaxis.set_major_locator(MultipleLocator(0.0002))
ax[2].grid(visible=True, which='major', color='darkgray', linestyle='-')
ax[2].grid(visible=True, which='minor', color='lightgray', linestyle='--')

ax[3].axhline(0, color='gray', linewidth=1, linestyle='--')
ax[3].plot(t, np.median(tilt_cw, axis=1), 'mediumblue', linewidth=2)
ax[3].plot(t, np.median(stretch_cw, axis=1), 'deepskyblue', linewidth=2)
ax[3].plot(t, np.median(bcl_cw, axis=1), 'red', linewidth=2)
ax[3].plot(t, np.median(fric_cw, axis=1), 'goldenrod', linewidth=2)
ax[3].plot(tp, np.median(total_cw, axis=1), 'k', linewidth=2)
ax[3].set_title("Parcel \u03c9$_{CW}$ tendency")
ax[3].set_xlabel('Time (min)', fontsize=12)
ax[3].set_xlim([t[0], t[-1]])
# ax[3].set_ylim([-0.0002, 0.0002])
ax[3].set_ylim([-0.0008, 0.0008])
ax[3].xaxis.set_major_locator(MultipleLocator(5))
ax[3].xaxis.set_minor_locator(MultipleLocator(1))
ax[3].yaxis.set_major_locator(MultipleLocator(0.0002))
ax[3].grid(visible=True, which='major', color='darkgray', linestyle='-')
ax[3].grid(visible=True, which='minor', color='lightgray', linestyle='--')

plt.suptitle(f"Parcel vorticity tendency at {mvtime} min", fontsize=12)

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






#%% Time series - both times, panels z/h/sw/cw, all terms

t = np.linspace(195, 225, 31)

# load parcel times
ds = nc.Dataset('/Volumes/Promise_Pegasus_70TB/merger/merger-125m/cm1out_pdata.nc')
ptime = ds.variables['time'][:].data
ds.close()

it1 = np.where(ptime/60 == 195)[0][0]
it2 = np.where(ptime/60 == 210)[0][0]
it3 = np.where(ptime/60 == 225)[0][0]

# load 210 min tendencies
dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/traj_vten_210min.pkl", 'rb')
vten1 = pickle.load(dbfile)
dbfile.close()
tilt_z1 = vten1['tilt_z']; stretch_z1 = vten1['stretch_z']; bcl_z1 = vten1['bcl_z']; fric_z1 = vten1['fric_z']
tilt_h1 = vten1['tilt_h']; stretch_h1 = vten1['stretch_h']; bcl_h1 = vten1['bcl_h']; fric_h1 = vten1['fric_h']
tilt_sw1 = vten1['tilt_sw']; stretch_sw1 = vten1['stretch_sw']; bcl_sw1 = vten1['bcl_sw']; fric_sw1 = vten1['fric_sw']
tilt_cw1 = vten1['tilt_cw']; stretch_cw1 = vten1['stretch_cw']; bcl_cw1 = vten1['bcl_cw']; fric_cw1 = vten1['fric_cw']
# load 225 min tendencies
dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/traj_vten_225min.pkl", 'rb')
vten2 = pickle.load(dbfile)
dbfile.close()
tilt_z2 = vten2['tilt_z']; stretch_z2 = vten2['stretch_z']; bcl_z2 = vten2['bcl_z']; fric_z2 = vten2['fric_z']
tilt_h2 = vten2['tilt_h']; stretch_h2 = vten2['stretch_h']; bcl_h2 = vten2['bcl_h']; fric_h2 = vten2['fric_h']
tilt_sw2 = vten2['tilt_sw']; stretch_sw2 = vten2['stretch_sw']; bcl_sw2 = vten2['bcl_sw']; fric_sw2 = vten2['fric_sw']
tilt_cw2 = vten2['tilt_cw']; stretch_cw2 = vten2['stretch_cw']; bcl_cw2 = vten2['bcl_cw']; fric_cw2 = vten2['fric_cw']

# load 210 min source regions
dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/traj_clusters/traj_clusters_210min_v2.pkl", 'rb')
cc = pickle.load(dbfile)
cc1 = cc['mv1']
dbfile.close()
# load 225 min source regions
dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/traj_clusters/traj_clusters_225min_v2.pkl", 'rb')
cc = pickle.load(dbfile)
cc2 = cc['mv1']
dbfile.close()

# load parcel vertical vorticity
dbfile = open(f"/Users/morgan.schneider/Documents/merger/traj_MV1.pkl", 'rb')
traj = pickle.load(dbfile)
dbfile.close()
zvort1 = traj[f"210min"]['zvort'][it1:it2+1,(cc1 == 1)]
zvort2 = traj[f"225min"]['zvort'][it2:it3+1,(cc2 == 1)]
# load 210 min parcel horizontal vorticity
dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/hvort_traj_210min.pkl", 'rb')
vort_traj = pickle.load(dbfile)
hvort1 = vort_traj['hvort_ml'][it1:it2+1,:]
svort1 = vort_traj['vort_sw_ml'][it1:it2+1,:]
cvort1 = vort_traj['vort_cw_ml_signed'][it1:it2+1,:]
dbfile.close()
# load 225 min parcel horizontal vorticity
dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/hvort_traj_225min.pkl", 'rb')
vort_traj = pickle.load(dbfile)
hvort2 = vort_traj['hvort_ml'][it2:it3+1,:]
svort2 = vort_traj['vort_sw_ml'][it2:it3+1,:]
cvort2 = vort_traj['vort_cw_ml_signed'][it2:it3+1,:]
dbfile.close()


if True:
    tsb_z1 = tilt_z1 + stretch_z1 + bcl_z1 + fric_z1
    tsb_h1 = tilt_h1 + stretch_h1 + bcl_h1 + fric_h1
    tsb_sw1 = tilt_sw1 + stretch_sw1 + bcl_sw1 + fric_sw1
    tsb_cw1 = tilt_cw1 + stretch_cw1 + bcl_cw1 + fric_cw1
    tsb_z2 = tilt_z2 + stretch_z2 + bcl_z2 + fric_z2
    tsb_h2 = tilt_h2 + stretch_h2 + bcl_h2 + fric_h2
    tsb_sw2 = tilt_sw2 + stretch_sw2 + bcl_sw2 + fric_sw2
    tsb_cw2 = tilt_cw2 + stretch_cw2 + bcl_cw2 + fric_cw2

if False:
    # it1 = np.where(ptime/60 == 195)[0][0]
    # it2 = np.where(ptime/60 == 210)[0][0]
    # it3 = np.where(ptime/60 == 225)[0][0]
    
    total_z1 = np.gradient(zvort1, ptime[it1:it2+1], axis=0)
    total_h1 = np.gradient(hvort1, ptime[it1:it2+1], axis=0)
    total_sw1 = np.gradient(svort1, ptime[it1:it2+1], axis=0)
    total_cw1 = np.gradient(cvort1, ptime[it1:it2+1], axis=0)
    fric_z1 = total_z1[::4,:] - tilt_z1 - stretch_z1 - bcl_z1
    fric_h1 = total_h1[::4,:] - tilt_h1 - stretch_h1 - bcl_h1
    fric_sw1 = total_sw1[::4,:] - tilt_sw1 - stretch_sw1 - bcl_sw1
    fric_cw1 = total_cw1[::4,:] - tilt_cw1 - stretch_cw1 - bcl_cw1
    
    total_z2 = np.gradient(zvort2, ptime[it2:it3+1], axis=0)
    total_h2 = np.gradient(hvort2, ptime[it2:it3+1], axis=0)
    total_sw2 = np.gradient(svort2, ptime[it2:it3+1], axis=0)
    total_cw2 = np.gradient(cvort2, ptime[it2:it3+1], axis=0)
    fric_z2 = total_z2[::4,:] - tilt_z2 - stretch_z2 - bcl_z2
    fric_h2 = total_h2[::4,:] - tilt_h2 - stretch_h2 - bcl_h2
    fric_sw2 = total_sw2[::4,:] - tilt_sw2 - stretch_sw2 - bcl_sw2
    fric_cw2 = total_cw2[::4,:] - tilt_cw2 - stretch_cw2 - bcl_cw2
    
    total_z1 = total_z1[::4,:]
    total_h1 = total_h1[::4,:]
    total_sw1 = total_sw1[::4,:]
    total_cw1 = total_cw1[::4,:]
    total_z2 = total_z2[::4,:]
    total_h2 = total_h2[::4,:]
    total_sw2 = total_sw2[::4,:]
    total_cw2 = total_cw2[::4,:]
    
if True:
    # it1 = np.where(ptime/60 == 195)[0][0]
    # it2 = np.where(ptime/60 == 210)[0][0]
    # it3 = np.where(ptime/60 == 225)[0][0]
    tp = ptime[it1:it3+1]/60
    
    total_z1 = np.gradient(zvort1[::4,:], ptime[it1:it2+1][::4], axis=0)
    total_h1 = np.gradient(hvort1[::4,:], ptime[it1:it2+1][::4], axis=0)
    total_sw1 = np.gradient(svort1[::4,:], ptime[it1:it2+1][::4], axis=0)
    total_cw1 = np.gradient(cvort1[::4,:], ptime[it1:it2+1][::4], axis=0)
    fric_z1 = total_z1 - tilt_z1 - stretch_z1 - bcl_z1
    fric_h1 = total_h1 - tilt_h1 - stretch_h1 - bcl_h1
    fric_sw1 = total_sw1 - tilt_sw1 - stretch_sw1 - bcl_sw1
    fric_cw1 = total_cw1 - tilt_cw1 - stretch_cw1 - bcl_cw1
    
    total_z2 = np.gradient(zvort2[::4,:], ptime[it2:it3+1][::4], axis=0)
    total_h2 = np.gradient(hvort2[::4,:], ptime[it2:it3+1][::4], axis=0)
    total_sw2 = np.gradient(svort2[::4,:], ptime[it2:it3+1][::4], axis=0)
    total_cw2 = np.gradient(cvort2[::4,:], ptime[it2:it3+1][::4], axis=0)
    
    fric_z2 = total_z2 - tilt_z2 - stretch_z2 - bcl_z2
    fric_h2 = total_h2 - tilt_h2 - stretch_h2 - bcl_h2
    fric_sw2 = total_sw2 - tilt_sw2 - stretch_sw2 - bcl_sw2
    fric_cw2 = total_cw2 - tilt_cw2 - stretch_cw2 - bcl_cw2


figsave = False

if True:
    fig,ax = plt.subplots(4, 1, figsize=(9,11), sharex=True, layout='constrained')
    
    ax[0].axhline(0, color='gray', linewidth=1.5, linestyle='--')
    s1, = ax[0].plot(t[:16], np.median(tilt_z1, axis=1), 'mediumblue', linewidth=2)
    s2, = ax[0].plot(t[:16], np.median(stretch_z1, axis=1), 'deepskyblue', linewidth=2)
    s3, = ax[0].plot(t[:16], np.median(bcl_z1, axis=1), 'red', linewidth=2)
    s4, = ax[0].plot(t[:16], np.median(fric_z1, axis=1), 'goldenrod', linewidth=2)
    s5, = ax[0].plot(t[:16], np.median(total_z1, axis=1), 'k', linewidth=3)
    s6, = ax[0].plot(t[:16], np.median(tsb_z1, axis=1), 'gray', linewidth=2)
    ax[0].plot(t[15:], np.median(tilt_z2, axis=1), 'mediumblue', linewidth=2)
    ax[0].plot(t[15:], np.median(stretch_z2, axis=1), 'deepskyblue', linewidth=2)
    ax[0].plot(t[15:], np.median(bcl_z2, axis=1), 'red', linewidth=2)
    ax[0].plot(t[15:], np.median(fric_z2, axis=1), 'goldenrod', linewidth=2)
    ax[0].plot(t[15:], np.median(total_z2, axis=1), 'k', linewidth=3)
    ax[0].plot(t[15:], np.median(tsb_z2, axis=1), 'gray', linewidth=2)
    ax[0].axvline(210, color='k', linewidth=1.5, linestyle='--')
    ax[0].set_title("Parcel \u03B6 tendency")
    ax[0].set_ylim([-0.0004, 0.0004])
    ax[0].legend(handles=[s1,s2,s3,s4,s5,s6], labels=['Tilting','Stretching','Baroclinic','Residual','Total','T+S+B'], loc=2, ncol=2)
    ax[0].grid(visible=True, which='major', color='darkgray', linestyle='-')
    ax[0].grid(visible=True, which='minor', color='lightgray', linestyle='--')
    ax[0].yaxis.set_minor_locator(MultipleLocator(0.0001))
    
    ax[1].axhline(0, color='gray', linewidth=1.5, linestyle='--')
    s1, = ax[1].plot(t[:16], np.median(tilt_h1, axis=1), 'mediumblue', linewidth=2)
    s2, = ax[1].plot(t[:16], np.median(stretch_h1, axis=1), 'deepskyblue', linewidth=2)
    s3, = ax[1].plot(t[:16], np.median(bcl_h1, axis=1), 'red', linewidth=2)
    s4, = ax[1].plot(t[:16], np.median(fric_h1, axis=1), 'goldenrod', linewidth=2)
    s5, = ax[1].plot(t[:16], np.median(total_h1, axis=1), 'k', linewidth=3)
    s6, = ax[1].plot(t[:16], np.median(tsb_h1, axis=1), 'gray', linewidth=2)
    ax[1].plot(t[15:], np.median(tilt_h2, axis=1), 'mediumblue', linewidth=2)
    ax[1].plot(t[15:], np.median(stretch_h2, axis=1), 'deepskyblue', linewidth=2)
    ax[1].plot(t[15:], np.median(bcl_h2, axis=1), 'red', linewidth=2)
    ax[1].plot(t[15:], np.median(fric_h2, axis=1), 'goldenrod', linewidth=2)
    ax[1].plot(t[15:], np.median(total_h2, axis=1), 'k', linewidth=3)
    ax[1].plot(t[15:], np.median(tsb_h2, axis=1), 'gray', linewidth=2)
    ax[1].axvline(210, color='k', linewidth=1.5, linestyle='--')
    ax[1].set_title("Parcel \u03c9$_H$ tendency")
    ax[1].set_ylim([-0.0004, 0.0008])
    ax[1].legend(handles=[s1,s2,s3,s4,s5,s6], labels=['Tilting','Stretching','Baroclinic','Residual','Total','T+S+B'], loc=2, ncol=2)
    ax[1].grid(visible=True, which='major', color='darkgray', linestyle='-')
    ax[1].grid(visible=True, which='minor', color='lightgray', linestyle='--')
    ax[1].yaxis.set_minor_locator(MultipleLocator(0.0001))
    
    ax[2].axhline(0, color='gray', linewidth=1.5, linestyle='--')
    s1, = ax[2].plot(t[:16], np.median(tilt_sw1, axis=1), 'mediumblue', linewidth=2)
    s2, = ax[2].plot(t[:16], np.median(stretch_sw1, axis=1), 'deepskyblue', linewidth=2)
    s3, = ax[2].plot(t[:16], np.median(bcl_sw1, axis=1), 'red', linewidth=2)
    s4, = ax[2].plot(t[:16], np.median(fric_sw1, axis=1), 'goldenrod', linewidth=2)
    s5, = ax[2].plot(t[:16], np.median(total_sw1, axis=1), 'k', linewidth=3)
    s6, = ax[2].plot(t[:16], np.median(tsb_sw1, axis=1), 'gray', linewidth=2)
    ax[2].plot(t[15:], np.median(tilt_sw2, axis=1), 'mediumblue', linewidth=2)
    ax[2].plot(t[15:], np.median(stretch_sw2, axis=1), 'deepskyblue', linewidth=2)
    ax[2].plot(t[15:], np.median(bcl_sw2, axis=1), 'red', linewidth=2)
    ax[2].plot(t[15:], np.median(fric_sw2, axis=1), 'goldenrod', linewidth=2)
    ax[2].plot(t[15:], np.median(total_sw2, axis=1), 'k', linewidth=3)
    ax[2].plot(t[15:], np.median(tsb_sw2, axis=1), 'gray', linewidth=2)
    ax[2].axvline(210, color='k', linewidth=1.5, linestyle='--')
    ax[2].set_title("Parcel \u03c9$_{SW}$ tendency")
    ax[2].set_ylim([-0.0008, 0.0006])
    ax[2].legend(handles=[s1,s2,s3,s4,s5,s6], labels=['Tilting','Stretching','Baroclinic','Residual','Total','T+S+B'], loc=2, ncol=2)
    ax[2].grid(visible=True, which='major', color='darkgray', linestyle='-')
    ax[2].grid(visible=True, which='minor', color='lightgray', linestyle='--')
    ax[2].yaxis.set_minor_locator(MultipleLocator(0.0001))
    
    ax[3].axhline(0, color='gray', linewidth=1.5, linestyle='--')
    s1, = ax[3].plot(t[:16], np.median(tilt_cw1, axis=1), 'mediumblue', linewidth=2)
    s2, = ax[3].plot(t[:16], np.median(stretch_cw1, axis=1), 'deepskyblue', linewidth=2)
    s3, = ax[3].plot(t[:16], np.median(bcl_cw1, axis=1), 'red', linewidth=2)
    s4, = ax[3].plot(t[:16], np.median(fric_cw1, axis=1), 'goldenrod', linewidth=2)
    s5, = ax[3].plot(t[:16], np.median(total_cw1, axis=1), 'k', linewidth=3)
    s6, = ax[3].plot(t[:16], np.median(tsb_cw1, axis=1), 'gray', linewidth=2)
    ax[3].plot(t[15:], np.median(tilt_cw2, axis=1), 'mediumblue', linewidth=2)
    ax[3].plot(t[15:], np.median(stretch_cw2, axis=1), 'deepskyblue', linewidth=2)
    ax[3].plot(t[15:], np.median(bcl_cw2, axis=1), 'red', linewidth=2)
    ax[3].plot(t[15:], np.median(fric_cw2, axis=1), 'goldenrod', linewidth=2)
    ax[3].plot(t[15:], np.median(total_cw2, axis=1), 'k', linewidth=3)
    ax[3].plot(t[15:], np.median(tsb_cw2, axis=1), 'gray', linewidth=2)
    ax[3].axvline(210, color='k', linewidth=1.5, linestyle='--')
    ax[3].set_title("Parcel \u03c9$_{CW}$ tendency")
    ax[3].set_xlabel('Time (min)', fontsize=12)
    ax[3].set_xlim([195, 225])
    ax[3].set_ylim([-0.0006, 0.0004])
    ax[3].legend(handles=[s1,s2,s3,s4,s5,s6], labels=['Tilting','Stretching','Baroclinic','Residual','Total','T+S+B'], loc='lower left', ncol=2)
    ax[3].grid(visible=True, which='major', color='darkgray', linestyle='-')
    ax[3].grid(visible=True, which='minor', color='lightgray', linestyle='--')
    ax[3].xaxis.set_major_locator(MultipleLocator(5))
    ax[3].xaxis.set_minor_locator(MultipleLocator(1))
    ax[3].yaxis.set_minor_locator(MultipleLocator(0.0001))
    
    if figsave:
        plt.savefig(f"/Users/morgan.schneider/Documents/merger/merger-125m/vort_tendency/traj_vten_timeseries1.png", dpi=300)




if True:
    fig,ax = plt.subplots(4, 1, figsize=(9,11), sharex=True, layout='constrained')
    
    ax[0].axhline(0, color='gray', linewidth=1.5, linestyle='--')
    ax[0].fill_between(t[:16], np.percentile(total_z1, 25, axis=1),
                    np.percentile(total_z1, 75, axis=1), color='k', alpha=0.2)
    ax[0].fill_between(t[15:], np.percentile(total_z2, 25, axis=1),
                    np.percentile(total_z2, 75, axis=1), color='k', alpha=0.2)
    s1, = ax[0].plot(t[:16], np.median(tilt_z1, axis=1), 'mediumblue', linewidth=2)
    s2, = ax[0].plot(t[:16], np.median(stretch_z1, axis=1), 'deepskyblue', linewidth=2)
    s3, = ax[0].plot(t[:16], np.median(bcl_z1, axis=1), 'red', linewidth=2)
    s4, = ax[0].plot(t[:16], np.median(fric_z1, axis=1), 'goldenrod', linewidth=2)
    s5, = ax[0].plot(t[:16], np.median(total_z1, axis=1), 'k', linewidth=3)
    # s6, = ax[0].plot(t[:16], np.median(tsb_z1, axis=1), 'dimgray', linewidth=2)
    ax[0].plot(t[15:], np.median(tilt_z2, axis=1), 'mediumblue', linewidth=2)
    ax[0].plot(t[15:], np.median(stretch_z2, axis=1), 'deepskyblue', linewidth=2)
    ax[0].plot(t[15:], np.median(bcl_z2, axis=1), 'red', linewidth=2)
    ax[0].plot(t[15:], np.median(fric_z2, axis=1), 'goldenrod', linewidth=2)
    ax[0].plot(t[15:], np.median(total_z2, axis=1), 'k', linewidth=3)
    # ax[0].plot(t[15:], np.median(tsb_z2, axis=1), 'dimgray', linewidth=2)
    ax[0].axvline(210, color='k', linewidth=1.5, linestyle='--')
    ax[0].set_title("Parcel \u03B6 tendency")
    ax[0].set_ylim([-0.0004, 0.0004])
    ax[0].legend(handles=[s1,s2,s3,s4,s5], labels=['Tilting','Stretching','Baroclinic','Residual','Total'], loc=2, ncol=2)
    ax[0].grid(visible=True, which='major', color='darkgray', linestyle='-')
    ax[0].grid(visible=True, which='minor', color='lightgray', linestyle='--')
    ax[0].yaxis.set_minor_locator(MultipleLocator(0.0001))
    
    ax[1].axhline(0, color='gray', linewidth=1.5, linestyle='--')
    ax[1].fill_between(t[:16], np.percentile(total_h1, 25, axis=1),
                    np.percentile(total_h1, 75, axis=1), color='k', alpha=0.2)
    ax[1].fill_between(t[15:], np.percentile(total_h2, 25, axis=1),
                    np.percentile(total_h2, 75, axis=1), color='k', alpha=0.2)
    s1, = ax[1].plot(t[:16], np.median(tilt_h1, axis=1), 'mediumblue', linewidth=2)
    s2, = ax[1].plot(t[:16], np.median(stretch_h1, axis=1), 'deepskyblue', linewidth=2)
    s3, = ax[1].plot(t[:16], np.median(bcl_h1, axis=1), 'red', linewidth=2)
    s4, = ax[1].plot(t[:16], np.median(fric_h1, axis=1), 'goldenrod', linewidth=2)
    s5, = ax[1].plot(t[:16], np.median(total_h1, axis=1), 'k', linewidth=3)
    # s6, = ax[1].plot(t[:16], np.median(tsb_h1, axis=1), 'dimgray', linewidth=2)
    ax[1].plot(t[15:], np.median(tilt_h2, axis=1), 'mediumblue', linewidth=2)
    ax[1].plot(t[15:], np.median(stretch_h2, axis=1), 'deepskyblue', linewidth=2)
    ax[1].plot(t[15:], np.median(bcl_h2, axis=1), 'red', linewidth=2)
    ax[1].plot(t[15:], np.median(fric_h2, axis=1), 'goldenrod', linewidth=2)
    ax[1].plot(t[15:], np.median(total_h2, axis=1), 'k', linewidth=3)
    # ax[1].plot(t[15:], np.median(tsb_h2, axis=1), 'dimgray', linewidth=2)
    ax[1].axvline(210, color='k', linewidth=1.5, linestyle='--')
    ax[1].set_title("Parcel \u03c9$_H$ tendency")
    ax[1].set_ylim([-0.0004, 0.0008])
    ax[1].legend(handles=[s1,s2,s3,s4,s5], labels=['Tilting','Stretching','Baroclinic','Residual','Total'], loc=2, ncol=2)
    ax[1].grid(visible=True, which='major', color='darkgray', linestyle='-')
    ax[1].grid(visible=True, which='minor', color='lightgray', linestyle='--')
    ax[1].yaxis.set_minor_locator(MultipleLocator(0.0001))
    
    ax[2].axhline(0, color='gray', linewidth=1.5, linestyle='--')
    ax[2].fill_between(t[:16], np.percentile(total_sw1, 25, axis=1),
                    np.percentile(total_sw1, 75, axis=1), color='k', alpha=0.2)
    ax[2].fill_between(t[15:], np.percentile(total_sw2, 25, axis=1),
                    np.percentile(total_sw2, 75, axis=1), color='k', alpha=0.2)
    s1, = ax[2].plot(t[:16], np.median(tilt_sw1, axis=1), 'mediumblue', linewidth=2)
    s2, = ax[2].plot(t[:16], np.median(stretch_sw1, axis=1), 'deepskyblue', linewidth=2)
    s3, = ax[2].plot(t[:16], np.median(bcl_sw1, axis=1), 'red', linewidth=2)
    s4, = ax[2].plot(t[:16], np.median(fric_sw1, axis=1), 'goldenrod', linewidth=2)
    s5, = ax[2].plot(t[:16], np.median(total_sw1, axis=1), 'k', linewidth=3)
    # s6, = ax[2].plot(t[:16], np.median(tsb_sw1, axis=1), 'dimgray', linewidth=2)
    ax[2].plot(t[15:], np.median(tilt_sw2, axis=1), 'mediumblue', linewidth=2)
    ax[2].plot(t[15:], np.median(stretch_sw2, axis=1), 'deepskyblue', linewidth=2)
    ax[2].plot(t[15:], np.median(bcl_sw2, axis=1), 'red', linewidth=2)
    ax[2].plot(t[15:], np.median(fric_sw2, axis=1), 'goldenrod', linewidth=2)
    ax[2].plot(t[15:], np.median(total_sw2, axis=1), 'k', linewidth=3)
    # ax[2].plot(t[15:], np.median(tsb_sw2, axis=1), 'dimgray', linewidth=2)
    ax[2].axvline(210, color='k', linewidth=1.5, linestyle='--')
    ax[2].set_title("Parcel \u03c9$_{SW}$ tendency")
    ax[2].set_ylim([-0.0008, 0.0006])
    ax[2].legend(handles=[s1,s2,s3,s4,s5], labels=['Tilting','Stretching','Baroclinic','Residual','Total'], loc=2, ncol=2)
    ax[2].grid(visible=True, which='major', color='darkgray', linestyle='-')
    ax[2].grid(visible=True, which='minor', color='lightgray', linestyle='--')
    ax[2].yaxis.set_minor_locator(MultipleLocator(0.0001))
    
    ax[3].axhline(0, color='gray', linewidth=1.5, linestyle='--')
    ax[3].fill_between(t[:16], np.percentile(total_cw1, 25, axis=1),
                    np.percentile(total_cw1, 75, axis=1), color='k', alpha=0.2)
    ax[3].fill_between(t[15:], np.percentile(total_cw2, 25, axis=1),
                    np.percentile(total_cw2, 75, axis=1), color='k', alpha=0.2)
    s1, = ax[3].plot(t[:16], np.median(tilt_cw1, axis=1), 'mediumblue', linewidth=2)
    s2, = ax[3].plot(t[:16], np.median(stretch_cw1, axis=1), 'deepskyblue', linewidth=2)
    s3, = ax[3].plot(t[:16], np.median(bcl_cw1, axis=1), 'red', linewidth=2)
    s4, = ax[3].plot(t[:16], np.median(fric_cw1, axis=1), 'goldenrod', linewidth=2)
    s5, = ax[3].plot(t[:16], np.median(total_cw1, axis=1), 'k', linewidth=3)
    # s6, = ax[3].plot(t[:16], np.median(tsb_cw1, axis=1), 'dimgray', linewidth=2)
    ax[3].plot(t[15:], np.median(tilt_cw2, axis=1), 'mediumblue', linewidth=2)
    ax[3].plot(t[15:], np.median(stretch_cw2, axis=1), 'deepskyblue', linewidth=2)
    ax[3].plot(t[15:], np.median(bcl_cw2, axis=1), 'red', linewidth=2)
    ax[3].plot(t[15:], np.median(fric_cw2, axis=1), 'goldenrod', linewidth=2)
    ax[3].plot(t[15:], np.median(total_cw2, axis=1), 'k', linewidth=3)
    # ax[3].plot(t[15:], np.median(tsb_cw2, axis=1), 'dimgray', linewidth=2)
    ax[3].axvline(210, color='k', linewidth=1.5, linestyle='--')
    ax[3].set_title("Parcel \u03c9$_{CW}$ tendency")
    ax[3].set_xlabel('Time (min)', fontsize=12)
    ax[3].set_xlim([195, 225])
    ax[3].set_ylim([-0.0006, 0.0004])
    ax[3].legend(handles=[s1,s2,s3,s4,s5], labels=['Tilting','Stretching','Baroclinic','Residual','Total'], loc='lower left', ncol=2)
    ax[3].grid(visible=True, which='major', color='darkgray', linestyle='-')
    ax[3].grid(visible=True, which='minor', color='lightgray', linestyle='--')
    ax[3].xaxis.set_major_locator(MultipleLocator(5))
    ax[3].xaxis.set_minor_locator(MultipleLocator(1))
    ax[3].yaxis.set_minor_locator(MultipleLocator(0.0001))
    
    if figsave:
        plt.savefig(f"/Users/morgan.schneider/Documents/merger/merger-125m/vort_tendency/traj_vten_timeseries2.png", dpi=300)


#%% Time series - both times, individual plots per vorticity component

# actual vorticity
fig,ax = plt.subplots(1, 1, figsize=(8,4), layout='constrained')

ax.axhline(0, color='gray', linewidth=1.5, linestyle='--')
s1, = ax.plot(np.linspace(195, 210, 16), zvort1[::4], 'goldenrod', linewidth=2)
s3, = ax.plot(np.linspace(195, 210, 16), svort1[::4], 'deepskyblue', linewidth=2)
s4, = ax.plot(np.linspace(195, 210, 16), cvort1[::4], 'mediumblue', linewidth=2)
s2, = ax.plot(np.linspace(195, 210, 16), hvort1[::4], 'red', linewidth=2)
ax.plot(np.linspace(210, 225, 16), zvort2[::4], 'goldenrod', linewidth=2)
ax.plot(np.linspace(210, 225, 16), svort2[::4], 'deepskyblue', linewidth=2)
ax.plot(np.linspace(210, 225, 16), cvort2[::4], 'mediumblue', linewidth=2)
ax.plot(np.linspace(210, 225, 16), hvort2[::4], 'red', linewidth=2)
ax.axvline(210, color='k', linewidth=1.5, linestyle='--')
ax.set_title("Parcel vorticity")
ax.set_xlim([195, 225])
ax.set_ylim([-0.04, 0.08])
ax.legend(handles=[s1,s2,s3,s4], labels=['Vertical','Horizontal','Streamwise','Crosswise'], loc=2)
ax.grid(visible=True, which='major', color='darkgray', linestyle='-')
ax.grid(visible=True, which='minor', color='lightgray', linestyle='--')
ax.xaxis.set_major_locator(MultipleLocator(5))
ax.xaxis.set_minor_locator(MultipleLocator(1))


# fig,ax = plt.subplots(1, 1, figsize=(8,4), layout='constrained')

# ax.axhline(0, color='gray', linewidth=1.5, linestyle='--')
# s1, = ax.plot(np.linspace(195, 210, 61), zvort1, 'goldenrod', linewidth=2)
# s3, = ax.plot(np.linspace(195, 210, 61), svort1, 'deepskyblue', linewidth=2)
# s4, = ax.plot(np.linspace(195, 210, 61), cvort1, 'mediumblue', linewidth=2)
# s2, = ax.plot(np.linspace(195, 210, 61), hvort1, 'red', linewidth=2)
# ax.plot(np.linspace(210, 225, 61), zvort2, 'goldenrod', linewidth=2)
# ax.plot(np.linspace(210, 225, 61), svort2, 'deepskyblue', linewidth=2)
# ax.plot(np.linspace(210, 225, 61), cvort2, 'mediumblue', linewidth=2)
# ax.plot(np.linspace(210, 225, 61), hvort2, 'red', linewidth=2)
# ax.axvline(210, color='k', linewidth=1.5, linestyle='--')
# ax.set_title("Parcel vorticity")
# ax.set_xlim([195, 225])
# ax.set_ylim([-0.06, 0.08])
# ax.legend(handles=[s1,s2,s3,s4], labels=['Vertical','Horizontal','Streamwise','Crosswise'], loc=2)



# vertical tendency
fig,ax = plt.subplots(1, 1, figsize=(8,4), layout='constrained')

ax.axhline(0, color='gray', linewidth=1.5, linestyle='--')
s1, = ax.plot(t[:16], np.median(tilt_z1, axis=1), 'mediumblue', linewidth=2)
s2, = ax.plot(t[:16], np.median(stretch_z1, axis=1), 'deepskyblue', linewidth=2)
s3, = ax.plot(t[:16], np.median(bcl_z1, axis=1), 'red', linewidth=2)
s4, = ax.plot(t[:16], np.median(fric_z1, axis=1), 'goldenrod', linewidth=2)
ax.plot(t[15:], np.median(tilt_z2, axis=1), 'mediumblue', linewidth=2)
ax.plot(t[15:], np.median(stretch_z2, axis=1), 'deepskyblue', linewidth=2)
ax.plot(t[15:], np.median(bcl_z2, axis=1), 'red', linewidth=2)
ax.plot(t[15:], np.median(fric_z2, axis=1), 'goldenrod', linewidth=2)
ax.axvline(210, color='k', linewidth=1.5, linestyle='--')
ax.set_title("Parcel vertical vorticity tendency")
ax.set_xlim([195, 225])
ax.set_ylim([-0.0004, 0.0004])
ax.legend(handles=[s1,s2,s3,s4], labels=['Tilting','Stretching','Baroclinic','Friction'], loc=2)
ax.grid(visible=True, which='major', color='darkgray', linestyle='-')
ax.grid(visible=True, which='minor', color='lightgray', linestyle='--')
ax.xaxis.set_major_locator(MultipleLocator(5))
ax.xaxis.set_minor_locator(MultipleLocator(1))

# horizontal tendency
fig,ax = plt.subplots(1, 1, figsize=(8,4), layout='constrained')

ax.axhline(0, color='gray', linewidth=1.5, linestyle='--')
s1, = ax.plot(t[:16], np.median(tilt_h1, axis=1), 'mediumblue', linewidth=2)
s2, = ax.plot(t[:16], np.median(stretch_h1, axis=1), 'deepskyblue', linewidth=2)
s3, = ax.plot(t[:16], np.median(bcl_h1, axis=1), 'red', linewidth=2)
s4, = ax.plot(t[:16], np.median(fric_h1, axis=1), 'goldenrod', linewidth=2)
ax.plot(t[15:], np.median(tilt_h2, axis=1), 'mediumblue', linewidth=2)
ax.plot(t[15:], np.median(stretch_h2, axis=1), 'deepskyblue', linewidth=2)
ax.plot(t[15:], np.median(bcl_h2, axis=1), 'red', linewidth=2)
ax.plot(t[15:], np.median(fric_h2, axis=1), 'goldenrod', linewidth=2)
ax.axvline(210, color='k', linewidth=1.5, linestyle='--')
ax.set_title("Parcel horizontal vorticity tendency")
ax.set_xlim([195, 225])
ax.set_ylim([-0.0004, 0.0008])
ax.legend(handles=[s1,s2,s3,s4], labels=['Tilting','Stretching','Baroclinic','Friction'], loc=2)
ax.grid(visible=True, which='major', color='darkgray', linestyle='-')
ax.grid(visible=True, which='minor', color='lightgray', linestyle='--')
ax.xaxis.set_major_locator(MultipleLocator(5))
ax.xaxis.set_minor_locator(MultipleLocator(1))


# streamwise tendency
fig,ax = plt.subplots(1, 1, figsize=(8,4), layout='constrained')

ax.axhline(0, color='gray', linewidth=1.5, linestyle='--')
s1, = ax.plot(t[:16], np.median(tilt_sw1, axis=1), 'mediumblue', linewidth=2)
s2, = ax.plot(t[:16], np.median(stretch_sw1, axis=1), 'deepskyblue', linewidth=2)
s3, = ax.plot(t[:16], np.median(bcl_sw1, axis=1), 'red', linewidth=2)
s4, = ax.plot(t[:16], np.median(fric_sw1, axis=1), 'goldenrod', linewidth=2)
ax.plot(t[15:], np.median(tilt_sw2, axis=1), 'mediumblue', linewidth=2)
ax.plot(t[15:], np.median(stretch_sw2, axis=1), 'deepskyblue', linewidth=2)
ax.plot(t[15:], np.median(bcl_sw2, axis=1), 'red', linewidth=2)
ax.plot(t[15:], np.median(fric_sw2, axis=1), 'goldenrod', linewidth=2)
ax.axvline(210, color='k', linewidth=1.5, linestyle='--')
ax.set_title("Parcel streamwise vorticity tendency")
ax.set_xlim([195, 225])
ax.set_ylim([-0.0008, 0.0006])
ax.legend(handles=[s1,s2,s3,s4], labels=['Tilting','Stretching','Baroclinic','Friction'], loc=2)
ax.grid(visible=True, which='major', color='darkgray', linestyle='-')
ax.grid(visible=True, which='minor', color='lightgray', linestyle='--')
ax.xaxis.set_major_locator(MultipleLocator(5))
ax.xaxis.set_minor_locator(MultipleLocator(1))


# crosswise tendency
fig,ax = plt.subplots(1, 1, figsize=(8,4), layout='constrained')

ax.axhline(0, color='gray', linewidth=1.5, linestyle='--')
s1, = ax.plot(t[:16], np.median(tilt_cw1, axis=1), 'mediumblue', linewidth=2)
s2, = ax.plot(t[:16], np.median(stretch_cw1, axis=1), 'deepskyblue', linewidth=2)
s3, = ax.plot(t[:16], np.median(bcl_cw1, axis=1), 'red', linewidth=2)
s4, = ax.plot(t[:16], np.median(fric_cw1, axis=1), 'goldenrod', linewidth=2)
ax.plot(t[15:], np.median(tilt_cw2, axis=1), 'mediumblue', linewidth=2)
ax.plot(t[15:], np.median(stretch_cw2, axis=1), 'deepskyblue', linewidth=2)
ax.plot(t[15:], np.median(bcl_cw2, axis=1), 'red', linewidth=2)
ax.plot(t[15:], np.median(fric_cw2, axis=1), 'goldenrod', linewidth=2)
ax.axvline(210, color='k', linewidth=1.5, linestyle='--')
ax.set_title("Parcel crosswise vorticity tendency")
ax.set_xlim([195, 225])
ax.set_ylim([-0.0006, 0.0004])
ax.legend(handles=[s1,s2,s3,s4], labels=['Tilting','Stretching','Baroclinic','Friction'], loc=2)
ax.grid(visible=True, which='major', color='darkgray', linestyle='-')
ax.grid(visible=True, which='minor', color='lightgray', linestyle='--')
ax.xaxis.set_major_locator(MultipleLocator(5))
ax.xaxis.set_minor_locator(MultipleLocator(1))

plt.show()

#%% Time series - both times, individual plots per term for all components

# actual vorticity
fig,ax = plt.subplots(1, 1, figsize=(8,4), layout='constrained')

ax.axhline(0, color='gray', linewidth=1.5, linestyle='--')
s1, = ax.plot(np.linspace(195, 210, 16), zvort1[::4], 'goldenrod', linewidth=2)
s3, = ax.plot(np.linspace(195, 210, 16), svort1[::4], 'deepskyblue', linewidth=2)
s4, = ax.plot(np.linspace(195, 210, 16), cvort1[::4], 'mediumblue', linewidth=2)
s2, = ax.plot(np.linspace(195, 210, 16), hvort1[::4], 'red', linewidth=2)
ax.plot(np.linspace(210, 225, 16), zvort2[::4], 'goldenrod', linewidth=2)
ax.plot(np.linspace(210, 225, 16), svort2[::4], 'deepskyblue', linewidth=2)
ax.plot(np.linspace(210, 225, 16), cvort2[::4], 'mediumblue', linewidth=2)
ax.plot(np.linspace(210, 225, 16), hvort2[::4], 'red', linewidth=2)
ax.axvline(210, color='k', linewidth=1.5, linestyle='--')
ax.set_title("Parcel vorticity")
ax.set_xlim([195, 225])
ax.set_ylim([-0.04, 0.08])
ax.legend(handles=[s1,s2,s3,s4], labels=['Vertical','Horizontal','Streamwise','Crosswise'], loc=2)
ax.grid(visible=True, which='major', color='darkgray', linestyle='-')
ax.grid(visible=True, which='minor', color='lightgray', linestyle='--')
ax.xaxis.set_major_locator(MultipleLocator(5))
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(0.02))
ax.yaxis.set_minor_locator(MultipleLocator(0.01))


# tilting term
fig,ax = plt.subplots(1, 1, figsize=(8,4), layout='constrained')

ax.axhline(0, color='gray', linewidth=1.5, linestyle='--')
s1, = ax.plot(t[:16], np.median(tilt_z1, axis=1), 'goldenrod', linewidth=2)
s3, = ax.plot(t[:16], np.median(tilt_sw1, axis=1), 'deepskyblue', linewidth=2)
s4, = ax.plot(t[:16], np.median(tilt_cw1, axis=1), 'mediumblue', linewidth=2)
s2, = ax.plot(t[:16], np.median(tilt_h1, axis=1), 'red', linewidth=2)
ax.plot(t[15:], np.median(tilt_z2, axis=1), 'goldenrod', linewidth=2)
ax.plot(t[15:], np.median(tilt_sw2, axis=1), 'deepskyblue', linewidth=2)
ax.plot(t[15:], np.median(tilt_cw2, axis=1), 'mediumblue', linewidth=2)
ax.plot(t[15:], np.median(tilt_h2, axis=1), 'red', linewidth=2)
ax.axvline(210, color='k', linewidth=1.5, linestyle='--')
ax.set_title("Parcel vorticity tendency - Tilting term")
ax.set_xlim([195, 225])
ax.set_ylim([-0.0008, 0.0004])
# ax.set_ylim([-0.0004, 0.0003])
ax.legend(handles=[s1,s2,s3,s4], labels=['Vertical','Horizontal','Streamwise','Crosswise'], loc=2)
ax.grid(visible=True, which='major', color='darkgray', linestyle='-')
ax.grid(visible=True, which='minor', color='lightgray', linestyle='--')
ax.xaxis.set_major_locator(MultipleLocator(5))
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(MultipleLocator(0.0001))


# stretching term
fig,ax = plt.subplots(1, 1, figsize=(8,4), layout='constrained')

ax.axhline(0, color='gray', linewidth=1.5, linestyle='--')
s1, = ax.plot(t[:16], np.median(stretch_z1, axis=1), 'goldenrod', linewidth=2)
s3, = ax.plot(t[:16], np.median(stretch_sw1, axis=1), 'deepskyblue', linewidth=2)
s4, = ax.plot(t[:16], np.median(stretch_cw1, axis=1), 'mediumblue', linewidth=2)
s2, = ax.plot(t[:16], np.median(stretch_h1, axis=1), 'red', linewidth=2)
ax.plot(t[15:], np.median(stretch_z2, axis=1), 'goldenrod', linewidth=2)
ax.plot(t[15:], np.median(stretch_sw2, axis=1), 'deepskyblue', linewidth=2)
ax.plot(t[15:], np.median(stretch_cw2, axis=1), 'mediumblue', linewidth=2)
ax.plot(t[15:], np.median(stretch_h2, axis=1), 'red', linewidth=2)
ax.axvline(210, color='k', linewidth=1.5, linestyle='--')
ax.set_title("Parcel vorticity tendency - Stretching term")
ax.set_xlim([195, 225])
ax.set_ylim([-0.0006, 0.0008])
# ax.set_ylim([-0.0002, 0.0002])
ax.legend(handles=[s1,s2,s3,s4], labels=['Vertical','Horizontal','Streamwise','Crosswise'], loc=2)
ax.grid(visible=True, which='major', color='darkgray', linestyle='-')
ax.grid(visible=True, which='minor', color='lightgray', linestyle='--')
ax.xaxis.set_major_locator(MultipleLocator(5))
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(MultipleLocator(0.0001))


# baroclinic term
fig,ax = plt.subplots(1, 1, figsize=(8,4), layout='constrained')

ax.axhline(0, color='gray', linewidth=1.5, linestyle='--')
s1, = ax.plot(t[:16], np.median(bcl_z1, axis=1), 'goldenrod', linewidth=2)
s3, = ax.plot(t[:16], np.median(bcl_sw1, axis=1), 'deepskyblue', linewidth=2)
s4, = ax.plot(t[:16], np.median(bcl_cw1, axis=1), 'mediumblue', linewidth=2)
s2, = ax.plot(t[:16], np.median(bcl_h1, axis=1), 'red', linewidth=2)
ax.plot(t[15:], np.median(bcl_z2, axis=1), 'goldenrod', linewidth=2)
ax.plot(t[15:], np.median(bcl_sw2, axis=1), 'deepskyblue', linewidth=2)
ax.plot(t[15:], np.median(bcl_cw2, axis=1), 'mediumblue', linewidth=2)
ax.plot(t[15:], np.median(bcl_h2, axis=1), 'red', linewidth=2)
ax.axvline(210, color='k', linewidth=1.5, linestyle='--')
ax.set_title("Parcel vorticity tendency - Baroclinic term")
ax.set_xlim([195, 225])
ax.set_ylim([-0.0002, 0.0003])
# ax.set_ylim([-0.00005, 0.00005])
ax.legend(handles=[s1,s2,s3,s4], labels=['Vertical','Horizontal','Streamwise','Crosswise'], loc=2)
ax.grid(visible=True, which='major', color='darkgray', linestyle='-')
ax.grid(visible=True, which='minor', color='lightgray', linestyle='--')
ax.xaxis.set_major_locator(MultipleLocator(5))
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(MultipleLocator(0.000025))


# friction term
fig,ax = plt.subplots(1, 1, figsize=(8,4), layout='constrained')

ax.axhline(0, color='gray', linewidth=1.5, linestyle='--')
s1, = ax.plot(t[:16], np.median(fric_z1, axis=1), 'goldenrod', linewidth=2)
s3, = ax.plot(t[:16], np.median(fric_sw1, axis=1), 'deepskyblue', linewidth=2)
s4, = ax.plot(t[:16], np.median(fric_cw1, axis=1), 'mediumblue', linewidth=2)
s2, = ax.plot(t[:16], np.median(fric_h1, axis=1), 'red', linewidth=2)
ax.plot(t[15:], np.median(fric_z2, axis=1), 'goldenrod', linewidth=2)
ax.plot(t[15:], np.median(fric_sw2, axis=1), 'deepskyblue', linewidth=2)
ax.plot(t[15:], np.median(fric_cw2, axis=1), 'mediumblue', linewidth=2)
ax.plot(t[15:], np.median(fric_h2, axis=1), 'red', linewidth=2)
ax.axvline(210, color='k', linewidth=1.5, linestyle='--')
ax.set_title("Parcel vorticity tendency - Friction term")
ax.set_xlim([195, 225])
ax.set_ylim([-1.5e-10, 2e-10])
# ax.set_ylim([-0.5e-10, 0.5e-10])
ax.legend(handles=[s1,s2,s3,s4], labels=['Vertical','Horizontal','Streamwise','Crosswise'], loc=2)
ax.grid(visible=True, which='major', color='darkgray', linestyle='-')
ax.grid(visible=True, which='minor', color='lightgray', linestyle='--')
ax.xaxis.set_major_locator(MultipleLocator(5))
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(MultipleLocator(0.25e-10))

plt.show()







#%% Calculate tendency plan views - median parcel

dbfile = open('/Users/morgan.schneider/Documents/merger/merger-125m/storm_motion.pkl', 'rb')
sm = pickle.load(dbfile)
u_storm = sm['u_storm']
v_storm = sm['v_storm']
dbfile.close()


calc_stretching = False
calc_tilting = False
calc_baroclinic = False
calc_friction = False

ip = '/Users/morgan.schneider/Documents/merger/merger-125m/cross_sections/MV1_vten/'

# Use median parcel
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
        stretch_cw = (-v_sr/ws_sr) * stretch_x + (u_sr/ws_sr) * stretch_y
        
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
        tilt_cw = (-v_sr/ws_sr) * tilt_x + (u_sr/ws_sr) * tilt_y
        
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
        bcl_cw = (-v_sr/ws_sr) * bcl_x + (u_sr/ws_sr) * bcl_y
        
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
        fric_cw = (-v_sr/ws_sr) * fric_x + (u_sr/ws_sr) * fric_y
        
        dat = {'x':xh[ix], 'y':yh[iy], 'z':zh[iz],
               'fric_x':fric_x, 'fric_y':fric_y, 'fric_z':fric_z,
               'fric_h':fric_h, 'fric_sw':fric_sw, 'fric_cw':fric_cw}
        save_to_pickle(dat, ip+f"plan{stime/60:.0f}_friction.pkl")
        del dat,fric_x,fric_y,fric_z,fric_h,fric_sw,fric_cw
    
    del u,v,w,xvort,yvort,zvort,u_sr,v_sr,ws_sr,hvort







#%% Calculate tendency plan views - parcel-centered composites ***FOR PAPER FIGS***

#         210 min   220 min
#
# min     204 m     184 m
# max     998 m     988 m
# mean    628 m     587 m
# median  698 m     599 m
# 25%     341 m     394 m
# 75%     867 m     740 m

# Total pids ml: 65 pids at 210 min, 51 pids at 220 min
# within +/- 100 m: 13 pids at 210 min (20%) - 44th %ile to 62nd %ile ->  -6 to +12
#                   18 pids at 220 min (35%) - 38th %ile to 72nd %ile -> -12 to +22
# within +/- 150 m: 21 pids at 210 min (32%) - 42nd %ile to 72nd %ile ->  -8 to +22
#                   22 pids at 220 min (43%) - 34th %ile to 76th %ile -> -16 to +26
# within +/- 200 m: 30 pids at 210 min (46%) - 36th %ile to 80th %ile -> -14 to +30
#                   27 pids at 220 min (53%) - 26th %ile to 78th %ile -> -24 to +28



mvtime = 220

if mvtime == 210:
    fnum = 43
elif mvtime == 220:
    fnum = 53

# ds = nc.Dataset('/Volumes/Promise_Pegasus_70TB/merger/merger-125m/cm1out_000013.nc')
ds = nc.Dataset('/Users/morgan.schneider/Documents/merger/merger-125m/temporary/cm1out_000043.nc')
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
zh = ds.variables['z'][:].data
ds.close()

xh = np.linspace(-179.9375, 179.9375, 2880)
yh = np.linspace(-179.9375, 179.9375, 2880)

dbfile = open('/Users/morgan.schneider/Documents/merger/merger-125m/storm_motion.pkl', 'rb')
sm = pickle.load(dbfile)
u_storm = sm['u_storm']
v_storm = sm['v_storm']
dbfile.close()

dbfile = open(f"/Users/morgan.schneider/Documents/merger/traj_MV1.pkl", 'rb')
traj = pickle.load(dbfile)
dbfile.close()

# ds = nc.Dataset('/Volumes/Promise_Pegasus_70TB/merger/merger-125m/cm1out_pdata.nc')
ds = nc.Dataset('/Users/morgan.schneider/Documents/merger/merger-125m/temporary/cm1out_pdata.nc')
ptime = ds.variables['time'][:].data
ds.close()

dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/traj_clusters_{mvtime}min_v2.pkl", 'rb')
ccs = pickle.load(dbfile)
cc = ccs['mv1']
dbfile.close()

pids_ml = traj[f"{mvtime}min"]['pids'][(cc == 1)]
x_ml = traj[f"{mvtime}min"]['x'][:,(cc == 1)]/1000
y_ml = traj[f"{mvtime}min"]['y'][:,(cc == 1)]/1000
z_ml = traj[f"{mvtime}min"]['z'][:,(cc == 1)]/1000

# # for using a specific layer
# ti = np.where(ptime == mvtime*60)[0][0]
# zm = np.round(np.median(z_ml[ti,:]), decimals=-1)
# dz = 150

# pids_ml = pids_ml[(z_ml[ti,:] >= zm-dz) & (z_ml[ti,:] <= zm+dz)]
# x_ml = x_ml[:, (z_ml[ti,:] >= zm-dz) & (z_ml[ti,:] <= zm+dz)]
# y_ml = y_ml[:, (z_ml[ti,:] >= zm-dz) & (z_ml[ti,:] <= zm+dz)]
# z_ml = z_ml[:, (z_ml[ti,:] >= zm-dz) & (z_ml[ti,:] <= zm+dz)]


ip = f"/Users/morgan.schneider/Documents/merger/merger-125m/"


times = np.zeros(shape=(11,), dtype=float)

sx = np.zeros(shape=(len(pids_ml), 11, 33, 33), dtype=float)
sy = np.zeros(shape=(len(pids_ml), 11, 33, 33), dtype=float)
sz = np.zeros(shape=(len(pids_ml), 11, 33, 33), dtype=float)
sh = np.zeros(shape=(len(pids_ml), 11, 33, 33), dtype=float)
ssw = np.zeros(shape=(len(pids_ml), 11, 33, 33), dtype=float)
scw = np.zeros(shape=(len(pids_ml), 11, 33, 33), dtype=float)
tx = np.zeros(shape=(len(pids_ml), 11, 33, 33), dtype=float)
ty = np.zeros(shape=(len(pids_ml), 11, 33, 33), dtype=float)
tz = np.zeros(shape=(len(pids_ml), 11, 33, 33), dtype=float)
th = np.zeros(shape=(len(pids_ml), 11, 33, 33), dtype=float)
tsw = np.zeros(shape=(len(pids_ml), 11, 33, 33), dtype=float)
tcw = np.zeros(shape=(len(pids_ml), 11, 33, 33), dtype=float)
bx = np.zeros(shape=(len(pids_ml), 11, 33, 33), dtype=float)
by = np.zeros(shape=(len(pids_ml), 11, 33, 33), dtype=float)
bh = np.zeros(shape=(len(pids_ml), 11, 33, 33), dtype=float)
bsw = np.zeros(shape=(len(pids_ml), 11, 33, 33), dtype=float)
bcw = np.zeros(shape=(len(pids_ml), 11, 33, 33), dtype=float)

tilt_components = False

if tilt_components: # individual tilting directions
    t_xy = np.zeros(shape=(len(pids_ml), 11, 33, 33), dtype=float)
    t_xz = np.zeros(shape=(len(pids_ml), 11, 33, 33), dtype=float)
    t_yx = np.zeros(shape=(len(pids_ml), 11, 33, 33), dtype=float)
    t_yz = np.zeros(shape=(len(pids_ml), 11, 33, 33), dtype=float)
    t_zx = np.zeros(shape=(len(pids_ml), 11, 33, 33), dtype=float)
    t_zy = np.zeros(shape=(len(pids_ml), 11, 33, 33), dtype=float)
    t_zsw = np.zeros(shape=(len(pids_ml), 11, 33, 33), dtype=float)
    t_zcw = np.zeros(shape=(len(pids_ml), 11, 33, 33), dtype=float)


m = 0
for fn in np.arange(fnum-10,fnum+1):
    print(f"cm1out_{fn:06d}")
    
    # ds = nc.Dataset(f"/Volumes/Promise_Pegasus_70TB/merger/merger-125m/cm1out_{fn:06d}.nc")
    ds = nc.Dataset(f"/Users/morgan.schneider/Documents/merger/merger-125m/temporary/cm1out_{fn:06d}.nc")
    stime = ds.variables['time'][:].data[0]
    it = np.where(ptime == stime)[0][0]
    
    xmin = np.min(x_ml[it,:]); ix1 = np.abs(xh - xmin).argmin()
    xmax = np.max(x_ml[it,:]); ix2 = np.abs(xh - xmax).argmin()
    ymin = np.min(y_ml[it,:]); iy1 = np.abs(yh - ymin).argmin()
    ymax = np.max(y_ml[it,:]); iy2 = np.abs(yh - ymax).argmin()
    zmax = np.max(z_ml[it,:]); iz1 = np.abs(zh - zmax).argmin()
    ix = slice(ix1-20, ix2+21)
    iy = slice(iy1-20, iy2+21)
    iz = slice(0, iz1+8)
    
    u = ds.variables['uinterp'][:].data[0,iz,iy,ix]
    v = ds.variables['vinterp'][:].data[0,iz,iy,ix]
    w = ds.variables['winterp'][:].data[0,iz,iy,ix]
    xvort = ds.variables['xvort'][:].data[0,iz,iy,ix]
    yvort = ds.variables['yvort'][:].data[0,iz,iy,ix]
    zvort = ds.variables['zvort'][:].data[0,iz,iy,ix]
    rho = ds.variables['rho'][:].data[0,iz,iy,ix]
    prs = ds.variables['prs'][:].data[0,iz,iy,ix]
    
    u_sr = u - u_storm[fn-13]
    v_sr = v - v_storm[fn-13]
    ws_sr = np.sqrt(u_sr**2 + v_sr**2)
    hvort = np.sqrt(xvort**2 + yvort**2)
    svort = (u_sr/ws_sr) * xvort + (v_sr/ws_sr) * yvort
    cvort = (-v_sr/ws_sr) * xvort + (u_sr/ws_sr) * yvort
    ds.close()
    
    
    dudx = np.gradient(u, xh[ix]*1000, axis=2)
    dudy = np.gradient(u, yh[iy]*1000, axis=1)
    dudz = np.gradient(u, zh[iz]*1000, axis=0)
    dvdx = np.gradient(v, xh[ix]*1000, axis=2)
    dvdy = np.gradient(v, yh[iy]*1000, axis=1)
    dvdz = np.gradient(v, zh[iz]*1000, axis=0)
    dwdx = np.gradient(w, xh[ix]*1000, axis=2)
    dwdy = np.gradient(w, yh[iy]*1000, axis=1)
    dwdz = np.gradient(w, zh[iz]*1000, axis=0)
    drdx = np.gradient(rho, xh[ix]*1000, axis=2)
    drdy = np.gradient(rho, yh[iy]*1000, axis=1)
    drdz = np.gradient(rho, zh[iz]*1000, axis=0)
    dpdx = np.gradient(prs, xh[ix]*1000, axis=2)
    dpdy = np.gradient(prs, yh[iy]*1000, axis=1)
    dpdz = np.gradient(prs, zh[iz]*1000, axis=0)
    del u,v,w,rho,prs
    
    # zz,yy,xx = np.meshgrid(zh[iz]*1000, yh[iy]*1000, xh[ix]*1000, indexing='ij')
    # du = mc.gradient(u, coordinates=(zz,yy,xx))
    # dv = mc.gradient(v, coordinates=(zz,yy,xx))
    # dw = mc.gradient(w, coordinates=(zz,yy,xx))
    # drho = mc.gradient(rho, coordinates=(zz,yy,xx))
    # dprs = mc.gradient(prs, coordinates=(zz,yy,xx))
    # dudx = du[2]; dudy = du[1]; dudz = du[0]
    # dvdx = dv[2]; dvdy = dv[1]; dvdz = dv[0]
    # dwdx = dw[2]; dwdy = dw[1]; dwdz = dw[0]
    # drdx = drho[2]; drdy = drho[1]; drdz = drho[0]
    # dpdx = dprs[2]; dpdy = dprs[1]; dpdz = dprs[0]
    # del u,v,w,rho,prs,du,dv,dw,drho,dprs
    
    
    for p in range(len(pids_ml)):
        xp = x_ml[it,p]
        yp = y_ml[it,p]
        zp = z_ml[it,p]
        
        ixp = np.abs(xh[ix]-xp).argmin()
        iyp = np.abs(yh[iy]-yp).argmin()
        k = np.abs(zh[iz]-zp).argmin()
        i = slice(ixp-16,ixp+17)
        j = slice(iyp-16,iyp+17)
        
        
        # Stretching
        sx[p,m,:,:] = xvort[k,j,i] * dudx[k,j,i]
        sy[p,m,:,:] = yvort[k,j,i] * dvdy[k,j,i]
        sz[p,m,:,:] = zvort[k,j,i] * dwdz[k,j,i]
        sh[p,m,:,:] = (xvort[k,j,i]/hvort[k,j,i]) * sx[p,m,:,:] + (yvort[k,j,i]/hvort[k,j,i]) * sy[p,m,:,:]
        ssw[p,m,:,:] = (u_sr[k,j,i]/ws_sr[k,j,i]) * sx[p,m,:,:] + (v_sr[k,j,i]/ws_sr[k,j,i]) * sy[p,m,:,:]
        scw[p,m,:,:] = (-v_sr[k,j,i]/ws_sr[k,j,i]) * sx[p,m,:,:] + (u_sr[k,j,i]/ws_sr[k,j,i]) * sy[p,m,:,:]
        
        # Tilting
        tx[p,m,:,:] = yvort[k,j,i] * dudy[k,j,i] + zvort[k,j,i] * dudz[k,j,i]
        ty[p,m,:,:] = xvort[k,j,i] * dvdx[k,j,i] + zvort[k,j,i] * dvdz[k,j,i]
        tz[p,m,:,:] = xvort[k,j,i] * dwdx[k,j,i] + yvort[k,j,i] * dwdy[k,j,i]
        th[p,m,:,:] = (xvort[k,j,i]/hvort[k,j,i]) * tx[p,m,:,:] + (yvort[k,j,i]/hvort[k,j,i]) * ty[p,m,:,:]
        tsw[p,m,:,:] = (u_sr[k,j,i]/ws_sr[k,j,i]) * tx[p,m,:,:] + (v_sr[k,j,i]/ws_sr[k,j,i]) * ty[p,m,:,:]
        tcw[p,m,:,:] = (-v_sr[k,j,i]/ws_sr[k,j,i]) * tx[p,m,:,:] + (u_sr[k,j,i]/ws_sr[k,j,i]) * ty[p,m,:,:]
        
        if tilt_components:
            # x vorticity components
            t_xy[p,m,:,:] = yvort[k,j,i] * dudy[k,j,i]
            t_xz[p,m,:,:] = zvort[k,j,i] * dudz[k,j,i]
            # y vorticity components
            t_yx[p,m,:,:] = xvort[k,j,i] * dvdx[k,j,i]
            t_yz[p,m,:,:] = zvort[k,j,i] * dvdz[k,j,i]
            # z vorticity components
            t_zx[p,m,:,:] = xvort[k,j,i] * dwdx[k,j,i]
            t_zy[p,m,:,:] = yvort[k,j,i] * dwdy[k,j,i]
            t_zsw[p,m,:,:] = svort[k,j,i] * ((u_sr[k,j,i]/ws_sr[k,j,i]) * dwdx[k,j,i] + 
                                           (v_sr[k,j,i]/ws_sr[k,j,i]) * dwdy[k,j,i])
            t_zcw[p,m,:,:] = cvort[k,j,i] * ((-v_sr[k,j,i]/ws_sr[k,j,i]) * dwdx[k,j,i] +
                                           (u_sr[k,j,i]/ws_sr[k,j,i]) * dwdy[k,j,i])
            
        # Baroclinic
        bx[p,m,:,:] = (1/1.1)**2 * (drdy[k,j,i] * dpdz[k,j,i] - drdz[k,j,i] * dpdy[k,j,i])
        by[p,m,:,:] = (1/1.1)**2 * (drdz[k,j,i] * dpdx[k,j,i] - drdx[k,j,i] * dpdz[k,j,i])
        bh[p,m,:,:] = (xvort[k,j,i]/hvort[k,j,i]) * bx[p,m,:,:] + (yvort[k,j,i]/hvort[k,j,i]) * by[p,m,:,:]
        bsw[p,m,:,:] = (u_sr[k,j,i]/ws_sr[k,j,i]) * bx[p,m,:,:] + (v_sr[k,j,i]/ws_sr[k,j,i]) * by[p,m,:,:]
        bcw[p,m,:,:] = (-v_sr[k,j,i]/ws_sr[k,j,i]) * bx[p,m,:,:] + (u_sr[k,j,i]/ws_sr[k,j,i]) * by[p,m,:,:]
        
    times[m] = stime/60
    m = m + 1
    
    del u_sr,v_sr,ws_sr,xvort,yvort,zvort,hvort,svort,cvort
    del dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,drdx,drdy,drdz,dpdx,dpdy,dpdz


# Save with individual parcels
if False:
    data = {'time':times, 'stretch_x':sx, 'stretch_y':sy, 'stretch_z':sz,
            'stretch_h':sh, 'stretch_sw':ssw, 'stretch_cw':scw,
            'tilt_x':tx, 'tilt_y':ty, 'tilt_z':tz, 'tilt_h':th,
            'tilt_sw':tsw, 'tilt_cw':tcw, 'bcl_x':bx, 'bcl_y':by,
            'bcl_h':bh, 'bcl_sw':bsw, 'bcl_cw':bcw}
    dbfile = open(ip+f"vten_traj_{mvtime}min_parcels.pkl", 'wb')
    pickle.dump(data, dbfile)
    dbfile.close()


stretch_x = np.mean(sx, axis=0)
stretch_y = np.mean(sy, axis=0)
stretch_z = np.mean(sz, axis=0)
stretch_h = np.mean(sh, axis=0)
stretch_sw = np.mean(ssw, axis=0)
stretch_cw = np.mean(scw, axis=0)
    
tilt_x = np.mean(tx, axis=0)
tilt_y = np.mean(ty, axis=0)
tilt_z = np.mean(tz, axis=0)
tilt_h = np.mean(th, axis=0)
tilt_sw = np.mean(tsw, axis=0)
tilt_cw = np.mean(tcw, axis=0)
    
bcl_x = np.mean(bx, axis=0)
bcl_y = np.mean(by, axis=0)
bcl_h = np.mean(bh, axis=0)
bcl_sw = np.mean(bsw, axis=0)
bcl_cw = np.mean(bcw, axis=0)

if tilt_components:
    tilt_xy = np.mean(t_xy, axis=0)
    tilt_xz = np.mean(t_xz, axis=0)
    tilt_yx = np.mean(t_yx, axis=0)
    tilt_yz = np.mean(t_yz, axis=0)
    tilt_zx = np.mean(t_zx, axis=0)
    tilt_zy = np.mean(t_zy, axis=0)
    tilt_zsw = np.mean(t_zsw, axis=0)
    tilt_zcw = np.mean(t_zcw, axis=0)
    #del t_xy,t_xz,t_yx,t_yz,t_zx,t_zy,t_zsw,t_zcw
    
# del sx,sy,sz,sh,ssw,scw,tx,ty,tz,th,tsw,tcw,bx,by,bh,bsw,bcw

# Save composited terms
if False:
    data = {'time':times, 'stretch_x':stretch_x, 'stretch_y':stretch_y, 'stretch_z':stretch_z,
            'stretch_h':stretch_h, 'stretch_sw':stretch_sw, 'stretch_cw':stretch_cw,
            'tilt_x':tilt_x, 'tilt_y':tilt_y, 'tilt_z':tilt_z, 'tilt_h':tilt_h,
            'tilt_sw':tilt_sw, 'tilt_cw':tilt_cw, 'bcl_x':bcl_x, 'bcl_y':bcl_y,
            'bcl_h':bcl_h, 'bcl_sw':bcl_sw, 'bcl_cw':bcl_cw}
    save_to_pickle(data, ip+f"vten_traj_{mvtime}min.pkl", new_pkl=True)




#%% Load data for plan views of tendency - median parcel

ip = '/Users/morgan.schneider/Documents/merger/merger-125m/'

t = 203
mvtime = 210

fnum = t - 180 + 13

dbfile = open(ip+f"plan{t}_tilt.pkl", 'rb')
tmp = pickle.load(dbfile)
xh = tmp['x']
yh = tmp['y']
zh = tmp['z']
tilt_x = tmp['tilt_x']
tilt_y = tmp['tilt_y']
tilt_z = tmp['tilt_z']
tilt_h = tmp['tilt_h']
tilt_sw = tmp['tilt_sw']
tilt_cw = tmp['tilt_cw']
dbfile.close()

dbfile = open(ip+f"plan{t}_stretch.pkl", 'rb')
tmp = pickle.load(dbfile)
stretch_x = tmp['stretch_x']
stretch_y = tmp['stretch_y']
stretch_z = tmp['stretch_z']
stretch_h = tmp['stretch_h']
stretch_sw = tmp['stretch_sw']
stretch_cw = tmp['stretch_cw']
dbfile.close()

dbfile = open(ip+f"plan{t}_baroclinic.pkl", 'rb')
tmp = pickle.load(dbfile)
bcl_x = tmp['bcl_x']
bcl_y = tmp['bcl_y']
bcl_z = tmp['bcl_z']
bcl_h = tmp['bcl_h']
bcl_sw = tmp['bcl_sw']
bcl_cw = tmp['bcl_cw']
dbfile.close()

dbfile = open(ip+f"plan{t}_friction.pkl", 'rb')
tmp = pickle.load(dbfile)
fric_x = tmp['fric_x']
fric_y = tmp['fric_y']
fric_z = tmp['fric_z']
fric_h = tmp['fric_h']
fric_sw = tmp['fric_sw']
fric_cw = tmp['fric_cw']
dbfile.close()


# Read parcel time
ds = nc.Dataset('/Volumes/Promise_Pegasus_70TB/merger/merger-125m/cm1out_pdata.nc')
ptime = ds.variables['time'][:].data
# xp = ds.variables['x'][:].data
# yp = ds.variables['y'][:].data
# zp = ds.variables['z'][:].data
ds.close()

# Load filtered parcel data
dbfile = open('/Users/morgan.schneider/Documents/merger/traj_MV1.pkl', 'rb')
traj = pickle.load(dbfile)
x_mv = traj[f"{mvtime}min"]['x']/1000
y_mv = traj[f"{mvtime}min"]['y']/1000
z_mv = traj[f"{mvtime}min"]['z']/1000
# w_mv = traj[f"{mvtime}min"]['w']
# zvort_mv = traj[f"{mvtime}min"]['zvort']
dbfile.close()

# Load source regions
dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/traj_clusters/traj_clusters_{mvtime}min_v2.pkl", 'rb')
tmp = pickle.load(dbfile)
cc = tmp['mv1']
dbfile.close()

# Mid-level source
x_ml = x_mv[:,(cc==1)]
y_ml = y_mv[:,(cc==1)]
z_ml = z_mv[:,(cc==1)]

# match model time and height to parcel time and median height
it = np.where(ptime/60 == t)[0][0]
iz = np.abs(zh - np.median(z_ml[it,:])).argmin()


# if np.abs(zh[iz] - np.median(z_ml[it,:])) > 0.05:
#     from scipy.interpolate import interp1d
    
#     tz_interp = interp1d(zh, tilt_z, kind='cubic', axis=0)
#     sz_interp = interp1d(zh, stretch_z, kind='cubic', axis=0)
#     bz_interp = interp1d(zh, bcl_z, kind='cubic', axis=0)
#     fz_interp = interp1d(zh, fric_z, kind='cubic', axis=0)
    
#     tilt_z_new = tz_interp(np.median(z_ml[it,:]))
#     stretch_z_new = sz_interp(np.median(z_ml[it,:]))
#     bcl_z_new = bz_interp(np.median(z_ml[it,:]))
#     fric_z_new = fz_interp(np.median(z_ml[it,:]))
    
#     th_interp = interp1d(zh, tilt_h, kind='cubic', axis=0)
#     sh_interp = interp1d(zh, stretch_h, kind='cubic', axis=0)
#     bh_interp = interp1d(zh, bcl_h, kind='cubic', axis=0)
#     fh_interp = interp1d(zh, fric_h, kind='cubic', axis=0)
    
#     tilt_h_new = th_interp(np.median(z_ml[it,:]))
#     stretch_h_new = sh_interp(np.median(z_ml[it,:]))
#     bcl_h_new = bh_interp(np.median(z_ml[it,:]))
#     fric_h_new = fh_interp(np.median(z_ml[it,:]))
    

ixp = np.abs(xh - np.median(x_ml[it,:])).argmin()
iyp = np.abs(yh - np.median(y_ml[it,:])).argmin()



# print(f"...zvort limits, {t} min, {zh[iz]*1000:.0f} m...")
# print(f"Tilting: {np.min(tilt_z[iz,:,:]):.6f}, {np.max(tilt_z[iz,:,:]):.6f}")
# print(f"Stretching: {np.min(stretch_z[iz,:,:]):.6f}, {np.max(stretch_z[iz,:,:]):.6f}")
# print(f"Baroclinic: {np.min(bcl_z[iz,:,:])}, {np.max(bcl_z[iz,:,:])}")
# print(f"Friction: {np.min(fric_z[iz,:,:])}, {np.max(fric_z[iz,:,:])}\n")

# print(f"...hvort limits, {t} min, {zh[iz]*1000:.0f} m...")
# print(f"Tilting: {np.min(tilt_h[iz,:,:]):.6f}, {np.max(tilt_h[iz,:,:]):.6f}")
# print(f"Stretching: {np.min(stretch_h[iz,:,:]):.6f}, {np.max(stretch_h[iz,:,:]):.6f}")
# print(f"Baroclinic: {np.min(bcl_h[iz,:,:])}, {np.max(bcl_h[iz,:,:])}")
# print(f"Friction: {np.min(fric_h[iz,:,:])}, {np.max(fric_h[iz,:,:])}")


# # Load model data from parcel time in a 5 km box roughly around median parcel location
# ds = nc.Dataset(f"/Volumes/Promise_Pegasus_70TB/merger/merger-125m/cm1out_{fnum:06d}.nc")
# ix1 = np.abs(ds.variables['xh'][:].data - xh[0]).argmin()
# ix2 = np.abs(ds.variables['xh'][:].data - xh[-1]).argmin() + 1
# iy1 = np.abs(ds.variables['yh'][:].data - yh[0]).argmin()
# iy2 = np.abs(ds.variables['yh'][:].data - yh[-1]).argmin() + 1
# dbz = ds.variables['dbz'][:].data[0,0,slice(iy1,iy2),slice(ix1,ix2)]
# w = np.max(ds.variables['winterp'][:].data[0,slice(0,iz+1),slice(iy1,iy2),slice(ix1,ix2)], axis=0)
# zvort = ds.variables['zvort'][:].data[0,iz,slice(iy1,iy2),slice(ix1,ix2)]
# hvort = np.sqrt(ds.variables['xvort'][:].data[0,iz,slice(iy1,iy2),slice(ix1,ix2)]**2 + 
#                 ds.variables['yvort'][:].data[0,iz,slice(iy1,iy2),slice(ix1,ix2)]**2)
# ds.close()





#%% Plot plan views of tendency terms - median parcel, single time

figsave = False


xlims = [-19, -14] # [-19,-14] at 205-210 / [-11,-5] at 220-225
ylims = [-73.5, -68.5] # [-72,-67]+0.75 at 205-210 / ?[-55,-49]+1 at 220-225

sz = 30

m = 4

tlims = [-0.002, 0.002]; tlims = [m*-1e-4, m*1e-4]; tlevs = np.linspace(tlims[0], tlims[1], 41) # 0.0008
slims = [-0.002, 0.002]; slims = [m*-1e-4, m*1e-4]; slevs = np.linspace(slims[0], slims[1], 41) # 0.0008
blims = [-1e-3, 1e-3]; blims = [m*-1e-4, m*1e-4]; blevs = np.linspace(blims[0], blims[1], 41) # 2e-4
flims = [-1e-9, 1e-9]; flims = [m*-1e-4, m*1e-4]; flevs = np.linspace(flims[0], flims[1], 41) # 2e-10

zvort_levs = [-0.04, -0.02, 0.02, 0.04]

tlevs = np.append(np.append([-0.01], tlevs), [0.01])
slevs = np.append(np.append([-0.01], slevs), [0.01])
blevs = np.append(np.append([-0.01], blevs), [0.01])


fig,ax = plt.subplots(2,2, figsize=(11,8), sharex=True, sharey=True, layout='constrained', subplot_kw=dict(box_aspect=1))

plot_contourf(xh, yh, tilt_z[iz,:,:], 'zvort', ax[0,0], levels=tlevs, datalims=tlims, xlims=xlims, ylims=ylims)
# ax[0,0].contour(xh, yh, w, levels=[5,10,15], colors=['gray','k','k'], linewidths=[1.5,1.5,3], linestyles='-')
# c = ax[0,0].contour(xh, yh, zvort, levels=zvort_levs, colors='k', linewidths=1)
# ax[0,0].clabel(c, fontsize=8, inline=True)
ax[0,0].scatter(x_ml[it,:], y_ml[it,:], s=sz, color='k', marker='.')
ax[0,0].set_title(f"Tilting")

plot_contourf(xh, yh, stretch_z[iz,:,:], 'zvort', ax[0,1], levels=slevs, datalims=slims, xlims=xlims, ylims=ylims)
# ax[0,1].contour(xh, yh, w, levels=[5,10,15], colors=['gray','k','k'], linewidths=[1.5,1.5,3], linestyles='-')
# ax[0,1].contour(xh, yh, zvort, levels=zvort_levs, colors='k', linewidths=1)
ax[0,1].scatter(x_ml[it,:], y_ml[it,:], s=sz, color='k', marker='.')
ax[0,1].set_title(f"Stretching")

plot_contourf(xh, yh, bcl_z[iz,:,:], 'zvort', ax[1,0], levels=blevs, datalims=blims, xlims=xlims, ylims=ylims)
# ax[1,0].contour(xh, yh, w, levels=[5,10,15], colors=['gray','k','k'], linewidths=[1.5,1.5,3], linestyles='-')
# ax[1,0].contour(xh, yh, zvort, levels=zvort_levs, colors='k', linewidths=1)
ax[1,0].scatter(x_ml[it,:], y_ml[it,:], s=sz, color='k', marker='.')
ax[1,0].set_title(f"Baroclinic")

plot_contourf(xh, yh, fric_z[iz,:,:], 'zvort', ax[1,1], levels=flevs, datalims=flims, xlims=xlims, ylims=ylims)
# ax[1,1].contour(xh, yh, w, levels=[5,10,15], colors=['gray','k','k'], linewidths=[1.5,1.5,3], linestyles='-')
# ax[1,1].contour(xh, yh, zvort, levels=zvort_levs, colors='k', linewidths=1)
ax[1,1].scatter(x_ml[it,:], y_ml[it,:], s=sz, color='k', marker='.')
ax[1,1].set_title(f"Friction")

plt.suptitle(f"\u03B6 tendency - {t} min, {zh[iz]*1000:.0f} m")

if figsave:
    plt.savefig(f"/Users/morgan.schneider/Documents/merger/merger-125m/vort_tendency/zvort_{t}min.png", dpi=300)


# tlims = [-0.001, 0.001]; tlevs = np.linspace(tlims[0], tlims[1], 21)
# slims = [-0.001, 0.001]; slevs = np.linspace(slims[0], slims[1], 21)
# blims = [-1e-3, 1e-3]; blevs = np.linspace(blims[0], blims[1], 21)
# flims = [-2e-9, 2e-9]; flevs = np.linspace(flims[0], flims[1], 21)


fig,ax = plt.subplots(2,2, figsize=(11,8), sharex=True, sharey=True, layout='constrained', subplot_kw=dict(box_aspect=1))

cf = plot_contourf(xh, yh, tilt_x[iz,:,:], 'xvort', ax[0,0], levels=tlevs, datalims=tlims, xlims=xlims, ylims=ylims, cmap='balance')
# ax[0,0].contour(xh, yh, w, levels=[5,10,15], colors=['gray','k','k'], linewidths=[1.5,1.5,3], linestyles='-')
# c = ax[0,0].contour(xh, yh, hvort, levels=hvort_levs, colors='k', linewidths=1)
ax[0,0].scatter(x_ml[it,:], y_ml[it,:], s=sz, color='k', marker='.')
ax[0,0].set_title(f"Tilting")

plot_contourf(xh, yh, stretch_x[iz,:,:], 'xvort', ax[0,1], levels=slevs, datalims=slims, xlims=xlims, ylims=ylims, cmap='balance')
# ax[0,1].contour(xh, yh, w, levels=[5,10,15], colors=['gray','k','k'], linewidths=[1.5,1.5,3], linestyles='-')
# ax[0,1].contour(xh, yh, hvort, levels=hvort_levs, colors='k', linewidths=1)
ax[0,1].scatter(x_ml[it,:], y_ml[it,:], s=sz, color='k', marker='.')
ax[0,1].set_title(f"Stretching")

plot_contourf(xh, yh, bcl_x[iz,:,:], 'xvort', ax[1,0], levels=blevs, datalims=blims, xlims=xlims, ylims=ylims, cmap='balance')
# ax[1,0].contour(xh, yh, w, levels=[5,10,15], colors=['gray','k','k'], linewidths=[1.5,1.5,3], linestyles='-')
# ax[1,0].contour(xh, yh, hvort, levels=hvort_levs, colors='k', linewidths=1)
ax[1,0].scatter(x_ml[it,:], y_ml[it,:], s=sz, color='k', marker='.')
ax[1,0].set_title(f"Baroclinic")

plot_contourf(xh, yh, fric_x[iz,:,:], 'xvort', ax[1,1], levels=flevs, datalims=flims, xlims=xlims, ylims=ylims, cmap='balance')
# ax[1,1].contour(xh, yh, w, levels=[5,10,15], colors=['gray','k','k'], linewidths=[1.5,1.5,3], linestyles='-')
# ax[1,1].contour(xh, yh, hvort, levels=hvort_levs, colors='k', linewidths=1)
ax[1,1].scatter(x_ml[it,:], y_ml[it,:], s=sz, color='k', marker='.')
ax[1,1].set_title(f"Friction")

plt.suptitle(f"\u03c9$_x$ tendency - {t} min, {zh[iz]*1000:.0f} m")

if figsave:
    plt.savefig(f"/Users/morgan.schneider/Documents/merger/merger-125m/vort_tendency/xvort_{t}min.png", dpi=300)


fig,ax = plt.subplots(2,2, figsize=(11,8), sharex=True, sharey=True, layout='constrained', subplot_kw=dict(box_aspect=1))

cf = plot_contourf(xh, yh, tilt_y[iz,:,:], 'yvort', ax[0,0], levels=tlevs, datalims=tlims, xlims=xlims, ylims=ylims, cmap='balance')
# ax[0,0].contour(xh, yh, w, levels=[5,10,15], colors=['gray','k','k'], linewidths=[1.5,1.5,3], linestyles='-')
# c = ax[0,0].contour(xh, yh, hvort, levels=hvort_levs, colors='k', linewidths=1)
ax[0,0].scatter(x_ml[it,:], y_ml[it,:], s=sz, color='k', marker='.')
ax[0,0].set_title(f"Tilting")

plot_contourf(xh, yh, stretch_y[iz,:,:], 'yvort', ax[0,1], levels=slevs, datalims=slims, xlims=xlims, ylims=ylims, cmap='balance')
# ax[0,1].contour(xh, yh, w, levels=[5,10,15], colors=['gray','k','k'], linewidths=[1.5,1.5,3], linestyles='-')
# ax[0,1].contour(xh, yh, hvort, levels=hvort_levs, colors='k', linewidths=1)
ax[0,1].scatter(x_ml[it,:], y_ml[it,:], s=sz, color='k', marker='.')
ax[0,1].set_title(f"Stretching")

plot_contourf(xh, yh, bcl_y[iz,:,:], 'yvort', ax[1,0], levels=blevs, datalims=blims, xlims=xlims, ylims=ylims, cmap='balance')
# ax[1,0].contour(xh, yh, w, levels=[5,10,15], colors=['gray','k','k'], linewidths=[1.5,1.5,3], linestyles='-')
# ax[1,0].contour(xh, yh, hvort, levels=hvort_levs, colors='k', linewidths=1)
ax[1,0].scatter(x_ml[it,:], y_ml[it,:], s=sz, color='k', marker='.')
ax[1,0].set_title(f"Baroclinic")

plot_contourf(xh, yh, fric_y[iz,:,:], 'yvort', ax[1,1], levels=flevs, datalims=flims, xlims=xlims, ylims=ylims, cmap='balance')
# ax[1,1].contour(xh, yh, w, levels=[5,10,15], colors=['gray','k','k'], linewidths=[1.5,1.5,3], linestyles='-')
# ax[1,1].contour(xh, yh, hvort, levels=hvort_levs, colors='k', linewidths=1)
ax[1,1].scatter(x_ml[it,:], y_ml[it,:], s=sz, color='k', marker='.')
ax[1,1].set_title(f"Friction")

plt.suptitle(f"\u03c9$_y$ tendency - {t} min, {zh[iz]*1000:.0f} m")

if figsave:
    plt.savefig(f"/Users/morgan.schneider/Documents/merger/merger-125m/vort_tendency/yvort_{t}min.png", dpi=300)





hvort_levs = [0.05, 0.1]

fig,ax = plt.subplots(2,2, figsize=(11,8), sharex=True, sharey=True, layout='constrained', subplot_kw=dict(box_aspect=1))

cf = plot_contourf(xh, yh, tilt_h[iz,:,:], 'hvort', ax[0,0], levels=tlevs, datalims=tlims, xlims=xlims, ylims=ylims, cmap='balance')
# ax[0,0].contour(xh, yh, w, levels=[5,10,15], colors=['gray','k','k'], linewidths=[1.5,1.5,3], linestyles='-')
# c = ax[0,0].contour(xh, yh, hvort, levels=hvort_levs, colors='k', linewidths=1)
ax[0,0].scatter(x_ml[it,:], y_ml[it,:], s=sz, color='k', marker='.')
ax[0,0].set_title(f"Tilting")

plot_contourf(xh, yh, stretch_h[iz,:,:], 'hvort', ax[0,1], levels=slevs, datalims=slims, xlims=xlims, ylims=ylims, cmap='balance')
# ax[0,1].contour(xh, yh, w, levels=[5,10,15], colors=['gray','k','k'], linewidths=[1.5,1.5,3], linestyles='-')
# ax[0,1].contour(xh, yh, hvort, levels=hvort_levs, colors='k', linewidths=1)
ax[0,1].scatter(x_ml[it,:], y_ml[it,:], s=sz, color='k', marker='.')
ax[0,1].set_title(f"Stretching")

plot_contourf(xh, yh, bcl_h[iz,:,:], 'hvort', ax[1,0], levels=blevs, datalims=blims, xlims=xlims, ylims=ylims, cmap='balance')
# ax[1,0].contour(xh, yh, w, levels=[5,10,15], colors=['gray','k','k'], linewidths=[1.5,1.5,3], linestyles='-')
# ax[1,0].contour(xh, yh, hvort, levels=hvort_levs, colors='k', linewidths=1)
ax[1,0].scatter(x_ml[it,:], y_ml[it,:], s=sz, color='k', marker='.')
ax[1,0].set_title(f"Baroclinic")

plot_contourf(xh, yh, fric_h[iz,:,:], 'hvort', ax[1,1], levels=flevs, datalims=flims, xlims=xlims, ylims=ylims, cmap='balance')
# ax[1,1].contour(xh, yh, w, levels=[5,10,15], colors=['gray','k','k'], linewidths=[1.5,1.5,3], linestyles='-')
# ax[1,1].contour(xh, yh, hvort, levels=hvort_levs, colors='k', linewidths=1)
ax[1,1].scatter(x_ml[it,:], y_ml[it,:], s=sz, color='k', marker='.')
ax[1,1].set_title(f"Friction")

plt.suptitle(f"\u03c9$_H$ tendency - {t} min, {zh[iz]*1000:.0f} m")

if figsave:
    plt.savefig(f"/Users/morgan.schneider/Documents/merger/merger-125m/vort_tendency/hvort_{t}min.png", dpi=300)







fig,ax = plt.subplots(2,2, figsize=(11,8), sharex=True, sharey=True, layout='constrained', subplot_kw=dict(box_aspect=1))

plot_contourf(xh, yh, tilt_sw[iz,:,:], 'hvort', ax[0,0], levels=tlevs, datalims=tlims, xlims=xlims, ylims=ylims, cmap='balance')
# ax[0,0].contour(xh, yh, w, levels=[5,10,15], colors=['gray','k','k'], linewidths=[1.5,1.5,3], linestyles='-')
# c = ax[0,0].contour(xh, yh, hvort, levels=hvort_levs, colors='k', linewidths=1)
# ax[0,0].clabel(c, fontsize=8, inline=True)
ax[0,0].scatter(x_ml[it,:], y_ml[it,:], s=sz, color='k', marker='.')
ax[0,0].set_title(f"Tilting")

plot_contourf(xh, yh, stretch_sw[iz,:,:], 'hvort', ax[0,1], levels=slevs, datalims=slims, xlims=xlims, ylims=ylims, cmap='balance')
# ax[0,1].contour(xh, yh, w, levels=[5,10,15], colors=['gray','k','k'], linewidths=[1.5,1.5,3], linestyles='-')
# ax[0,1].contour(xh, yh, hvort, levels=hvort_levs, colors='k', linewidths=1)
ax[0,1].scatter(x_ml[it,:], y_ml[it,:], s=sz, color='k', marker='.')
ax[0,1].set_title(f"Stretching")

plot_contourf(xh, yh, bcl_sw[iz,:,:], 'hvort', ax[1,0], levels=blevs, datalims=blims, xlims=xlims, ylims=ylims, cmap='balance')
# ax[1,0].contour(xh, yh, w, levels=[5,10,15], colors=['gray','k','k'], linewidths=[1.5,1.5,3], linestyles='-')
# ax[1,0].contour(xh, yh, hvort, levels=hvort_levs, colors='k', linewidths=1)
ax[1,0].scatter(x_ml[it,:], y_ml[it,:], s=sz, color='k', marker='.')
ax[1,0].set_title(f"Baroclinic")

plot_contourf(xh, yh, fric_sw[iz,:,:], 'hvort', ax[1,1], levels=flevs, datalims=flims, xlims=xlims, ylims=ylims, cmap='balance')
# ax[1,1].contour(xh, yh, w, levels=[5,10,15], colors=['gray','k','k'], linewidths=[1.5,1.5,3], linestyles='-')
# ax[1,1].contour(xh, yh, hvort, levels=hvort_levs, colors='k', linewidths=1)
ax[1,1].scatter(x_ml[it,:], y_ml[it,:], s=sz, color='k', marker='.')
ax[1,1].set_title(f"Friction")

plt.suptitle(f"Streamwise \u03c9$_H$ tendency - {t} min, {zh[iz]*1000:.0f} m")

if figsave:
    plt.savefig(f"/Users/morgan.schneider/Documents/merger/merger-125m/vort_tendency/swvort_{t}min.png", dpi=300)




fig,ax = plt.subplots(2,2, figsize=(11,8), sharex=True, sharey=True, layout='constrained', subplot_kw=dict(box_aspect=1))

plot_contourf(xh, yh, tilt_cw[iz,:,:], 'hvort', ax[0,0], levels=tlevs, datalims=tlims, xlims=xlims, ylims=ylims, cmap='balance')
# ax[0,0].contour(xh, yh, w, levels=[5,10,15], colors=['gray','k','k'], linewidths=[1.5,1.5,3], linestyles='-')
# c = ax[0,0].contour(xh, yh, hvort, levels=hvort_levs, colors='k', linewidths=1)
# ax[0,0].clabel(c, fontsize=8, inline=True)
ax[0,0].scatter(x_ml[it,:], y_ml[it,:], s=sz, color='k', marker='.')
ax[0,0].set_title(f"Tilting")

plot_contourf(xh, yh, stretch_cw[iz,:,:], 'hvort', ax[0,1], levels=slevs, datalims=slims, xlims=xlims, ylims=ylims, cmap='balance')
# ax[0,1].contour(xh, yh, w, levels=[5,10,15], colors=['gray','k','k'], linewidths=[1.5,1.5,3], linestyles='-')
# ax[0,1].contour(xh, yh, hvort, levels=hvort_levs, colors='k', linewidths=1)
ax[0,1].scatter(x_ml[it,:], y_ml[it,:], s=sz, color='k', marker='.')
ax[0,1].set_title(f"Stretching")

plot_contourf(xh, yh, bcl_cw[iz,:,:], 'hvort', ax[1,0], levels=blevs, datalims=blims, xlims=xlims, ylims=ylims, cmap='balance')
# ax[1,0].contour(xh, yh, w, levels=[5,10,15], colors=['gray','k','k'], linewidths=[1.5,1.5,3], linestyles='-')
# ax[1,0].contour(xh, yh, hvort, levels=hvort_levs, colors='k', linewidths=1)
ax[1,0].scatter(x_ml[it,:], y_ml[it,:], s=sz, color='k', marker='.')
ax[1,0].set_title(f"Baroclinic")

plot_contourf(xh, yh, fric_cw[iz,:,:], 'hvort', ax[1,1], levels=flevs, datalims=flims, xlims=xlims, ylims=ylims, cmap='balance')
# ax[1,1].contour(xh, yh, w, levels=[5,10,15], colors=['gray','k','k'], linewidths=[1.5,1.5,3], linestyles='-')
# ax[1,1].contour(xh, yh, hvort, levels=hvort_levs, colors='k', linewidths=1)
ax[1,1].scatter(x_ml[it,:], y_ml[it,:], s=sz, color='k', marker='.')
ax[1,1].set_title(f"Friction")

plt.suptitle(f"Crosswise \u03c9$_H$ tendency - {t} min, {zh[iz]*1000:.0f} m")

if figsave:
    plt.savefig(f"/Users/morgan.schneider/Documents/merger/merger-125m/vort_tendency/cwvort_{t}min.png", dpi=300)




#%% Load data for time composites of tendencies, single time - median parcel


ip = '/Users/morgan.schneider/Documents/merger/merger-125m/cross_sections/MV1_vten/'

times = np.arange(218, 220) # 203-208, 217-223/220-225, 214-219/218-221
mvtime = 220

# Read parcel time
ds = nc.Dataset('/Volumes/Promise_Pegasus_70TB/merger/merger-125m/cm1out_pdata.nc')
ptime = ds.variables['time'][:].data
ds.close()

# Load filtered parcel data
dbfile = open('/Users/morgan.schneider/Documents/merger/traj_MV1.pkl', 'rb')
traj = pickle.load(dbfile)
x_mv = traj[f"{mvtime}min"]['x']/1000
y_mv = traj[f"{mvtime}min"]['y']/1000
z_mv = traj[f"{mvtime}min"]['z']/1000
dbfile.close()

# Load source regions
dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/traj_clusters/traj_clusters_{mvtime}min_v2.pkl", 'rb')
tmp = pickle.load(dbfile)
cc = tmp['mv1']
dbfile.close()

# Mid-level source
x_ml = x_mv[:,(cc==1)]
y_ml = y_mv[:,(cc==1)]
z_ml = z_mv[:,(cc==1)]

x_median = np.median(x_ml, axis=1); x_max = np.max(x_ml, axis=1); x_min = np.min(x_ml, axis=1)
y_median = np.median(y_ml, axis=1); y_max = np.max(y_ml, axis=1); y_min = np.min(y_ml, axis=1)
z_median = np.median(z_ml, axis=1); z_max = np.max(z_ml, axis=1); z_min = np.min(z_ml, axis=1)


dbfile = open(ip+f"plan{mvtime}_tilt.pkl", 'rb')
tmp = pickle.load(dbfile)
xh = tmp['x']
yh = tmp['y']
zh = tmp['z']
dbfile.close()



tilt_x = np.zeros(shape=(len(times),17,17), dtype=float)
tilt_y = np.zeros(shape=(len(times),17,17), dtype=float)
tilt_z = np.zeros(shape=(len(times),17,17), dtype=float)
tilt_h = np.zeros(shape=(len(times),17,17), dtype=float)
tilt_sw = np.zeros(shape=(len(times),17,17), dtype=float)
tilt_cw = np.zeros(shape=(len(times),17,17), dtype=float)

stretch_x = np.zeros(shape=(len(times),17,17), dtype=float)
stretch_y = np.zeros(shape=(len(times),17,17), dtype=float)
stretch_z = np.zeros(shape=(len(times),17,17), dtype=float)
stretch_h = np.zeros(shape=(len(times),17,17), dtype=float)
stretch_sw = np.zeros(shape=(len(times),17,17), dtype=float)
stretch_cw = np.zeros(shape=(len(times),17,17), dtype=float)

bcl_x = np.zeros(shape=(len(times),17,17), dtype=float)
bcl_y = np.zeros(shape=(len(times),17,17), dtype=float)
bcl_h = np.zeros(shape=(len(times),17,17), dtype=float)
bcl_sw = np.zeros(shape=(len(times),17,17), dtype=float)
bcl_cw = np.zeros(shape=(len(times),17,17), dtype=float)


for i in range(len(times)):
    t = times[i]
    it = np.where(ptime/60 == t)[0][0]
    iz = np.abs(zh - z_median[it]).argmin()
    i1 = np.abs(xh - x_median[it]).argmin()
    i2 = np.abs(yh - y_median[it]).argmin()
    ix = slice(i1-8,i1+9)
    iy = slice(i2-8,i2+9)
    
    
    dbfile = open(ip+f"plan{t}_tilt.pkl", 'rb')
    tmp = pickle.load(dbfile)
    tilt_x[i,:,:] = tmp['tilt_x'][iz,iy,ix]
    tilt_y[i,:,:] = tmp['tilt_y'][iz,iy,ix]
    tilt_z[i,:,:] = tmp['tilt_z'][iz,iy,ix]
    tilt_h[i,:,:] = tmp['tilt_h'][iz,iy,ix]
    tilt_sw[i,:,:] = tmp['tilt_sw'][iz,iy,ix]
    tilt_cw[i,:,:] = tmp['tilt_cw'][iz,iy,ix]
    dbfile.close()
    
    dbfile = open(ip+f"plan{t}_stretch.pkl", 'rb')
    tmp = pickle.load(dbfile)
    stretch_x[i,:,:] = tmp['stretch_x'][iz,iy,ix]
    stretch_y[i,:,:] = tmp['stretch_y'][iz,iy,ix]
    stretch_z[i,:,:] = tmp['stretch_z'][iz,iy,ix]
    stretch_h[i,:,:] = tmp['stretch_h'][iz,iy,ix]
    stretch_sw[i,:,:] = tmp['stretch_sw'][iz,iy,ix]
    stretch_cw[i,:,:] = tmp['stretch_cw'][iz,iy,ix]
    dbfile.close()
    
    dbfile = open(ip+f"plan{t}_baroclinic.pkl", 'rb')
    tmp = pickle.load(dbfile)
    bcl_x[i,:,:] = tmp['bcl_x'][iz,iy,ix]
    bcl_y[i,:,:] = tmp['bcl_y'][iz,iy,ix]
    bcl_h[i,:,:] = tmp['bcl_h'][iz,iy,ix]
    bcl_sw[i,:,:] = tmp['bcl_sw'][iz,iy,ix]
    bcl_cw[i,:,:] = tmp['bcl_cw'][iz,iy,ix]
    dbfile.close()
    


tilt_x_comp = np.mean(tilt_x, axis=0)
tilt_y_comp = np.mean(tilt_y, axis=0)
tilt_z_comp = np.mean(tilt_z, axis=0)
tilt_h_comp = np.mean(tilt_h, axis=0)
tilt_sw_comp = np.mean(tilt_sw, axis=0)
tilt_cw_comp = np.mean(tilt_cw, axis=0)

stretch_x_comp = np.mean(stretch_x, axis=0)
stretch_y_comp = np.mean(stretch_y, axis=0)
stretch_z_comp = np.mean(stretch_z, axis=0)
stretch_h_comp = np.mean(stretch_h, axis=0)
stretch_sw_comp = np.mean(stretch_sw, axis=0)
stretch_cw_comp = np.mean(stretch_cw, axis=0)

bcl_x_comp = np.mean(bcl_x, axis=0)
bcl_y_comp = np.mean(bcl_y, axis=0)
bcl_h_comp = np.mean(bcl_h, axis=0)
bcl_sw_comp = np.mean(bcl_sw, axis=0)
bcl_cw_comp = np.mean(bcl_cw, axis=0)


itp = np.where(ptime/60 == times[-1])[0][0]

xp = x_ml[itp,:]
yp = y_ml[itp,:]
zp = z_ml[itp,:]
xmp = x_median[itp]
ymp = y_median[itp]
zmp = z_median[itp]

iz = np.abs(zh - zmp).argmin()
i1 = np.abs(xh - xmp).argmin()
i2 = np.abs(yh - ymp).argmin()
ix = slice(i1-8,i1+9)
iy = slice(i2-8,i2+9)


#%% Load data for time composites of tendencies, single time - parcel-centered composites

ip = '/Users/morgan.schneider/Documents/merger/merger-125m/cross_sections/MV1_vten/'

mvtime = 220

times = np.arange(214,219) # 203-208, 217-223/220-225, 214-219/218-221

# Read parcel time
ds = nc.Dataset('/Volumes/Promise_Pegasus_70TB/merger/merger-125m/cm1out_pdata.nc')
ptime = ds.variables['time'][:].data
ds.close()

# Load filtered parcel data
dbfile = open('/Users/morgan.schneider/Documents/merger/traj_MV1.pkl', 'rb')
traj = pickle.load(dbfile)
dbfile.close()

# Load source regions
dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/traj_clusters/traj_clusters_{mvtime}min_v2.pkl", 'rb')
tmp = pickle.load(dbfile)
cc = tmp['mv1']
dbfile.close()

# Mid-level source
x_ml = traj[f"{mvtime}min"]['x'][:,(cc==1)]/1000
y_ml = traj[f"{mvtime}min"]['y'][:,(cc==1)]/1000
z_ml = traj[f"{mvtime}min"]['z'][:,(cc==1)]/1000

x_median = np.median(x_ml, axis=1)
y_median = np.median(y_ml, axis=1)
z_median = np.median(z_ml, axis=1)


ds = nc.Dataset('/Volumes/Promise_Pegasus_70TB/merger/merger-125m/cm1out_000014.nc')
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
zh = ds.variables['z'][:].data
ds.close()

itp = np.where(ptime/60 == times[-1])[0][0]
ixp = np.abs(xh - x_median[itp]).argmin()
iyp = np.abs(yh - y_median[itp]).argmin()
ix = slice(ixp-8,ixp+9)
iy = slice(iyp-8,iyp+9)

xp = x_ml[itp,:]
yp = y_ml[itp,:]
zp = z_ml[itp,:]


dbfile = open(ip+f"vten_traj_{mvtime}min.pkl", 'rb')
vten = pickle.load(dbfile)
stime = vten['time']
dbfile.close()

# it1 = np.where(stime == times[0])[0][0]
# it2 = np.where(stime == times[-1])[0][0]
# it = slice(it1,it2)

# stretch_x_comp = np.mean(vten['stretch_x'][it,:,:], axis=0)
# stretch_y_comp = np.mean(vten['stretch_y'][it,:,:], axis=0)
# stretch_z_comp = np.mean(vten['stretch_z'][it,:,:], axis=0)
# stretch_h_comp = np.mean(vten['stretch_h'][it,:,:], axis=0)
# stretch_sw_comp = np.mean(vten['stretch_sw'][it,:,:], axis=0)
# stretch_cw_comp = np.mean(vten['stretch_cw'][it,:,:], axis=0)

# tilt_x_comp = np.mean(vten['tilt_x'][it,:,:], axis=0)
# tilt_y_comp = np.mean(vten['tilt_y'][it,:,:], axis=0)
# tilt_z_comp = np.mean(vten['tilt_z'][it,:,:], axis=0)
# tilt_h_comp = np.mean(vten['tilt_h'][it,:,:], axis=0)
# tilt_sw_comp = np.mean(vten['tilt_sw'][it,:,:], axis=0)
# tilt_cw_comp = np.mean(vten['tilt_cw'][it,:,:], axis=0)

# bcl_x_comp = np.mean(vten['bcl_x'][it,:,:], axis=0)
# bcl_y_comp = np.mean(vten['bcl_y'][it,:,:], axis=0)
# bcl_h_comp = np.mean(vten['bcl_h'][it,:,:], axis=0)
# bcl_sw_comp = np.mean(vten['bcl_sw'][it,:,:], axis=0)
# bcl_cw_comp = np.mean(vten['bcl_cw'][it,:,:], axis=0)

time = 219

it = np.where(stime == time)[0][0]

itp = np.where(ptime/60 == time)[0][0]
ixp = np.abs(xh - x_median[itp]).argmin()
iyp = np.abs(yh - y_median[itp]).argmin()
ix = slice(ixp-8,ixp+9)
iy = slice(iyp-8,iyp+9)

xp = x_ml[itp,:]
yp = y_ml[itp,:]
zp = z_ml[itp,:]

tilt_z_comp = vten['tilt_z'][it,:,:]
stretch_z_comp = vten['stretch_z'][it,:,:]
tilt_h_comp = vten['tilt_h'][it,:,:]
stretch_h_comp = vten['stretch_h'][it,:,:]
bcl_h_comp = vten['bcl_h'][it,:,:]


#%% Make plan view plots of tendency terms - parcel-centered composites, single time

figsave = False


xlims = [xh[ix][0], xh[ix][-1]] # [-19,-14] at 205-210 / [-11,-5] at 220-225
ylims = [yh[iy][0], yh[iy][-1]] # [-72,-67]+0.75 at 205-210 / ?[-55,-49]+1 at 220-225

sz = 80

m = 2

tlims = [m*-1e-4, m*1e-4]; tlevs = np.linspace(tlims[0], tlims[1], 41)
slims = [m*-1e-4, m*1e-4]; slevs = np.linspace(slims[0], slims[1], 41)
blims = [m*-1e-4, m*1e-4]; blevs = np.linspace(blims[0], blims[1], 41)

zvort_levs = [-0.04, -0.02, 0.02, 0.04]

tlevs = np.append(np.append([-0.01], tlevs), [0.01])
slevs = np.append(np.append([-0.01], slevs), [0.01])
blevs = np.append(np.append([-0.01], blevs), [0.01])




fig,ax = plt.subplots(1,2, figsize=(8.5,4), sharex=True, sharey=True, layout='constrained', subplot_kw=dict(box_aspect=1))

c = plot_contourf(xh[ix], yh[iy], tilt_z_comp, 'zvort', ax[0], levels=tlevs, datalims=tlims, xlims=xlims, ylims=ylims, cbar=False)
# cb = plt.colorbar(c, ax=ax[0], extend='both')
# cb.set_label("d\u03B6/dt (s$^{-2}$)", fontsize=10)
# cb.formatter.set_powerlimits((0,0))
ax[0].scatter(np.median(xp), np.median(yp), s=sz, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
# ax[0].set_title(f"Tilting")
ax[0].set_xlabel('x (km)', fontsize=12)
ax[0].set_ylabel('y (km)', fontsize=12)

c = plot_contourf(xh[ix], yh[iy], stretch_z_comp, 'zvort', ax[1], levels=slevs, datalims=slims, xlims=xlims, ylims=ylims, cbar=False)
cb = plt.colorbar(c, ax=ax[1], extend='both')
cb.set_label("d\u03B6/dt (s$^{-2}$)", fontsize=10)
cb.formatter.set_powerlimits((0,0))
ax[1].scatter(np.median(xp), np.median(yp), s=sz, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
# ax[1].set_title(f"Stretching")
ax[1].set_xlabel('x (km)', fontsize=12)

# plt.suptitle(f"Composite \u03B6 tendency ({times[0]}-{times[-1]} min)")
plt.suptitle(f"Composite \u03B6 tendency ({time} min)", fontsize=14)

if figsave:
    plt.savefig(f"/Users/morgan.schneider/Documents/merger/merger-125m/vort_tendency/comp_zvort_{time}min.png", dpi=300)


fig,ax = plt.subplots(1, 3, figsize=(11.5,4), sharex=True, sharey=True, layout='constrained', subplot_kw=dict(box_aspect=1))

c = plot_contourf(xh[ix], yh[iy], tilt_h_comp, 'zvort', ax[0], levels=tlevs, datalims=tlims, xlims=xlims, ylims=ylims, cbar=False)
# cb = plt.colorbar(c, ax=ax[0], extend='both')
# cb.set_label("d\u03B6/dt (s$^{-2}$)", fontsize=10)
# cb.formatter.set_powerlimits((0,0))
ax[0].scatter(np.median(xp), np.median(yp), s=sz, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
# ax[0].set_title(f"Tilting")
ax[0].set_xlabel('x (km)', fontsize=12)
ax[0].set_ylabel('y (km)', fontsize=12)

c = plot_contourf(xh[ix], yh[iy], stretch_h_comp, 'zvort', ax[1], levels=slevs, datalims=slims, xlims=xlims, ylims=ylims, cbar=False)
# cb = plt.colorbar(c, ax=ax[1], extend='both')
# cb.set_label("d\u03B7/dt (s$^{-2}$)", fontsize=10)
# cb.formatter.set_powerlimits((0,0))
ax[1].scatter(np.median(xp), np.median(yp), s=sz, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
# ax[1].set_title(f"Stretching")
ax[1].set_xlabel('x (km)', fontsize=12)

c = plot_contourf(xh[ix], yh[iy], bcl_h_comp, 'zvort', ax[2], levels=slevs, datalims=slims, xlims=xlims, ylims=ylims, cbar=False)
cb = plt.colorbar(c, ax=ax[2], extend='both')
cb.set_label("d|\u03c9$_H$|/dt (s$^{-2}$)", fontsize=12)
cb.formatter.set_powerlimits((0,0))
ax[2].scatter(np.median(xp), np.median(yp), s=sz, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
# ax[2].set_title(f"Baroclinic")
ax[2].set_xlabel('x (km)', fontsize=12)

plt.suptitle(f"Composite |\u03c9$_H$| tendency ({time} min)", fontsize=16)

if figsave:
    plt.savefig(f"/Users/morgan.schneider/Documents/merger/merger-125m/vort_tendency/comp_hvort_{time}min.png", dpi=300)



#%

fig,ax = plt.subplots(1, 3, figsize=(11.5,4), sharex=True, sharey=True, layout='constrained', subplot_kw=dict(box_aspect=1))

c = plot_contourf(xh[ix], yh[iy], tilt_x_comp, 'xvort', ax[0], levels=tlevs, datalims=tlims, xlims=xlims, ylims=ylims, cbar=False)
# cb = plt.colorbar(c, ax=ax[0], extend='both')
# cb.set_label("d\u03B6/dt (s$^{-2}$)", fontsize=10)
# cb.formatter.set_powerlimits((0,0))
ax[0].scatter(xp, yp, s=sz, color='k', marker='.')
ax[0].set_title(f"Tilting")
ax[0].set_xlabel('x (km)', fontsize=12)
ax[0].set_ylabel('y (km)', fontsize=12)

c = plot_contourf(xh[ix], yh[iy], stretch_x_comp, 'xvort', ax[1], levels=slevs, datalims=slims, xlims=xlims, ylims=ylims, cbar=False)
# cb = plt.colorbar(c, ax=ax[1], extend='both')
# cb.set_label("d\u03BE/dt (s$^{-2}$)", fontsize=10)
# cb.formatter.set_powerlimits((0,0))
ax[1].scatter(xp, yp, s=sz, color='k', marker='.')
ax[1].set_title(f"Stretching")
ax[1].set_xlabel('x (km)', fontsize=12)

c = plot_contourf(xh[ix], yh[iy], bcl_x_comp, 'xvort', ax[2], levels=slevs, datalims=slims, xlims=xlims, ylims=ylims, cbar=False)
cb = plt.colorbar(c, ax=ax[2], extend='both')
cb.set_label("d\u03BE/dt (s$^{-2}$)", fontsize=12)
cb.formatter.set_powerlimits((0,0))
ax[2].scatter(xp, yp, s=sz, color='k', marker='.')
ax[2].set_title(f"Baroclinic")
ax[2].set_xlabel('x (km)', fontsize=12)

plt.suptitle(f"Composite \u03BE tendency ({times[0]}-{times[-1]} min)")

if figsave:
    plt.savefig(f"/Users/morgan.schneider/Documents/merger/merger-125m/vort_tendency/xvort_comp_{mvtime}min.png", dpi=300)




fig,ax = plt.subplots(1, 3, figsize=(11.5,4), sharex=True, sharey=True, layout='constrained', subplot_kw=dict(box_aspect=1))

c = plot_contourf(xh[ix], yh[iy], tilt_y_comp, 'yvort', ax[0], levels=tlevs, datalims=tlims, xlims=xlims, ylims=ylims, cbar=False)
# cb = plt.colorbar(c, ax=ax[0], extend='both')
# cb.set_label("d\u03B6/dt (s$^{-2}$)", fontsize=10)
# cb.formatter.set_powerlimits((0,0))
ax[0].scatter(xp, yp, s=sz, color='k', marker='.')
ax[0].set_title(f"Tilting")
ax[0].set_xlabel('x (km)', fontsize=12)
ax[0].set_ylabel('y (km)', fontsize=12)

c = plot_contourf(xh[ix], yh[iy], stretch_y_comp, 'yvort', ax[1], levels=slevs, datalims=slims, xlims=xlims, ylims=ylims, cbar=False)
# cb = plt.colorbar(c, ax=ax[1], extend='both')
# cb.set_label("d\u03B7/dt (s$^{-2}$)", fontsize=10)
# cb.formatter.set_powerlimits((0,0))
ax[1].scatter(xp, yp, s=sz, color='k', marker='.')
ax[1].set_title(f"Stretching")
ax[1].set_xlabel('x (km)', fontsize=12)

c = plot_contourf(xh[ix], yh[iy], bcl_y_comp, 'yvort', ax[2], levels=slevs, datalims=slims, xlims=xlims, ylims=ylims, cbar=False)
cb = plt.colorbar(c, ax=ax[2], extend='both')
cb.set_label("d\u03B7/dt (s$^{-2}$)", fontsize=12)
cb.formatter.set_powerlimits((0,0))
ax[2].scatter(xp, yp, s=sz, color='k', marker='.')
ax[2].set_title(f"Baroclinic")
ax[2].set_xlabel('x (km)', fontsize=12)

plt.suptitle(f"Composite \u03B7 tendency ({times[0]}-{times[-1]} min)")

if figsave:
    plt.savefig(f"/Users/morgan.schneider/Documents/merger/merger-125m/vort_tendency/yvort_comp_{mvtime}min.png", dpi=300)



fig,ax = plt.subplots(1, 3, figsize=(11.5,4), sharex=True, sharey=True, layout='constrained', subplot_kw=dict(box_aspect=1))

c = plot_contourf(xh[ix], yh[iy], tilt_sw_comp, 'xvort', ax[0], levels=tlevs, datalims=tlims, xlims=xlims, ylims=ylims, cbar=False)
# cb = plt.colorbar(c, ax=ax[0], extend='both')
# cb.set_label("d\u03B6/dt (s$^{-2}$)", fontsize=10)
# cb.formatter.set_powerlimits((0,0))
ax[0].scatter(xp, yp, s=sz, color='k', marker='.')
ax[0].set_title(f"Tilting")
ax[0].set_xlabel('x (km)', fontsize=12)
ax[0].set_ylabel('y (km)', fontsize=12)

c = plot_contourf(xh[ix], yh[iy], stretch_sw_comp, 'xvort', ax[1], levels=slevs, datalims=slims, xlims=xlims, ylims=ylims, cbar=False)
# cb = plt.colorbar(c, ax=ax[1], extend='both')
# cb.set_label("d\u03B7/dt (s$^{-2}$)", fontsize=10)
# cb.formatter.set_powerlimits((0,0))
ax[1].scatter(xp, yp, s=sz, color='k', marker='.')
ax[1].set_title(f"Stretching")
ax[1].set_xlabel('x (km)', fontsize=12)

c = plot_contourf(xh[ix], yh[iy], bcl_sw_comp, 'xvort', ax[2], levels=slevs, datalims=slims, xlims=xlims, ylims=ylims, cbar=False)
cb = plt.colorbar(c, ax=ax[2], extend='both')
cb.set_label("d\u03c9$_{SW}$/dt (s$^{-2}$)", fontsize=12)
cb.formatter.set_powerlimits((0,0))
ax[2].scatter(xp, yp, s=sz, color='k', marker='.')
ax[2].set_title(f"Baroclinic")
ax[2].set_xlabel('x (km)', fontsize=12)

plt.suptitle(f"Composite streamwise \u03c9 tendency ({times[0]}-{times[-1]} min)")

if figsave:
    plt.savefig(f"/Users/morgan.schneider/Documents/merger/merger-125m/vort_tendency/swvort_comp_{mvtime}min.png", dpi=300)




fig,ax = plt.subplots(1, 3, figsize=(11.5,4), sharex=True, sharey=True, layout='constrained', subplot_kw=dict(box_aspect=1))

c = plot_contourf(xh[ix], yh[iy], tilt_cw_comp, 'yvort', ax[0], levels=tlevs, datalims=tlims, xlims=xlims, ylims=ylims, cbar=False)
# cb = plt.colorbar(c, ax=ax[0], extend='both')
# cb.set_label("d\u03B6/dt (s$^{-2}$)", fontsize=10)
# cb.formatter.set_powerlimits((0,0))
ax[0].scatter(xp, yp, s=sz, color='k', marker='.')
ax[0].set_title(f"Tilting")
ax[0].set_xlabel('x (km)', fontsize=12)
ax[0].set_ylabel('y (km)', fontsize=12)

c = plot_contourf(xh[ix], yh[iy], stretch_cw_comp, 'yvort', ax[1], levels=slevs, datalims=slims, xlims=xlims, ylims=ylims, cbar=False)
# cb = plt.colorbar(c, ax=ax[1], extend='both')
# cb.set_label("d\u03B7/dt (s$^{-2}$)", fontsize=10)
# cb.formatter.set_powerlimits((0,0))
ax[1].scatter(xp, yp, s=sz, color='k', marker='.')
ax[1].set_title(f"Stretching")
ax[1].set_xlabel('x (km)', fontsize=12)

c = plot_contourf(xh[ix], yh[iy], bcl_cw_comp, 'yvort', ax[2], levels=slevs, datalims=slims, xlims=xlims, ylims=ylims, cbar=False)
cb = plt.colorbar(c, ax=ax[2], extend='both')
cb.set_label("d\u03c9$_{CW}$/dt (s$^{-2}$)", fontsize=12)
cb.formatter.set_powerlimits((0,0))
ax[2].scatter(xp, yp, s=sz, color='k', marker='.')
ax[2].set_title(f"Baroclinic")
ax[2].set_xlabel('x (km)', fontsize=12)

plt.suptitle(f"Composite crosswise \u03c9 tendency ({times[0]}-{times[-1]} min)")

if figsave:
    plt.savefig(f"/Users/morgan.schneider/Documents/merger/merger-125m/vort_tendency/cwvort_comp_{mvtime}min.png", dpi=300)


#%% Load data for composite vorticity tendency plots, all times - median parcel

ip = '/Users/morgan.schneider/Documents/merger/merger-125m/cross_sections/MV1_vten/'


# Read parcel time
ds = nc.Dataset('/Volumes/Promise_Pegasus_70TB/merger/merger-125m/cm1out_pdata.nc')
ptime = ds.variables['time'][:].data
ds.close()

# Load filtered parcel data
dbfile = open('/Users/morgan.schneider/Documents/merger/traj_MV1.pkl', 'rb')
traj = pickle.load(dbfile)
x1_mv = traj[f"210min"]['x']/1000
y1_mv = traj[f"210min"]['y']/1000
z1_mv = traj[f"210min"]['z']/1000
x2_mv = traj[f"220min"]['x']/1000
y2_mv = traj[f"220min"]['y']/1000
z2_mv = traj[f"220min"]['z']/1000
dbfile.close()

# Load source regions and pull mid-level source
dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/traj_clusters/traj_clusters_210min_v2.pkl", 'rb')
tmp = pickle.load(dbfile)
cc = tmp['mv1']
dbfile.close()

x1_ml = x1_mv[:,(cc==1)]
y1_ml = y1_mv[:,(cc==1)]
z1_ml = z1_mv[:,(cc==1)]

dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/traj_clusters/traj_clusters_220min_v2.pkl", 'rb')
tmp = pickle.load(dbfile)
cc = tmp['mv1']
dbfile.close()

x2_ml = x2_mv[:,(cc==1)]
y2_ml = y2_mv[:,(cc==1)]
z2_ml = z2_mv[:,(cc==1)]


x1_median = np.median(x1_ml, axis=1); x1_max = np.max(x1_ml, axis=1); x1_min = np.min(x1_ml, axis=1)
y1_median = np.median(y1_ml, axis=1); y1_max = np.max(y1_ml, axis=1); y1_min = np.min(y1_ml, axis=1)
z1_median = np.median(z1_ml, axis=1); z1_max = np.max(z1_ml, axis=1); z1_min = np.min(z1_ml, axis=1)
x2_median = np.median(x2_ml, axis=1); x2_max = np.max(x2_ml, axis=1); x2_min = np.min(x2_ml, axis=1)
y2_median = np.median(y2_ml, axis=1); y2_max = np.max(y2_ml, axis=1); y2_min = np.min(y2_ml, axis=1)
z2_median = np.median(z2_ml, axis=1); z2_max = np.max(z2_ml, axis=1); z2_min = np.min(z2_ml, axis=1)


dbfile = open(ip+f"plan210_tilt.pkl", 'rb')
tmp = pickle.load(dbfile)
xh1 = tmp['x']
yh1 = tmp['y']
zh1 = tmp['z']
dbfile.close()

dbfile = open(ip+f"plan220_tilt.pkl", 'rb')
tmp = pickle.load(dbfile)
xh2 = tmp['x']
yh2 = tmp['y']
zh2 = tmp['z']
dbfile.close()




### Time 1 ###
times1 = np.arange(203, 208) # 210 min, downdraft
tilt_x = np.zeros(shape=(len(times1),17,17), dtype=float)
tilt_y = np.zeros(shape=(len(times1),17,17), dtype=float)
tilt_z = np.zeros(shape=(len(times1),17,17), dtype=float)
tilt_h = np.zeros(shape=(len(times1),17,17), dtype=float)
tilt_sw = np.zeros(shape=(len(times1),17,17), dtype=float)
tilt_cw = np.zeros(shape=(len(times1),17,17), dtype=float)
stretch_x = np.zeros(shape=(len(times1),17,17), dtype=float)
stretch_y = np.zeros(shape=(len(times1),17,17), dtype=float)
stretch_z = np.zeros(shape=(len(times1),17,17), dtype=float)
stretch_h = np.zeros(shape=(len(times1),17,17), dtype=float)
stretch_sw = np.zeros(shape=(len(times1),17,17), dtype=float)
stretch_cw = np.zeros(shape=(len(times1),17,17), dtype=float)
bcl_x = np.zeros(shape=(len(times1),17,17), dtype=float)
bcl_y = np.zeros(shape=(len(times1),17,17), dtype=float)
bcl_h = np.zeros(shape=(len(times1),17,17), dtype=float)
bcl_sw = np.zeros(shape=(len(times1),17,17), dtype=float)
bcl_cw = np.zeros(shape=(len(times1),17,17), dtype=float)

for i in range(len(times1)):
    t = times1[i]
    it = np.where(ptime/60 == t)[0][0]
    iz = np.abs(zh1 - z1_median[it]).argmin()
    i1 = np.abs(xh1 - x1_median[it]).argmin()
    i2 = np.abs(yh1 - y1_median[it]).argmin()
    ix = slice(i1-8,i1+9)
    iy = slice(i2-8,i2+9)
    
    dbfile = open(ip+f"plan{t}_tilt.pkl", 'rb')
    tmp = pickle.load(dbfile)
    tilt_x[i,:,:] = tmp['tilt_x'][iz,iy,ix]
    tilt_y[i,:,:] = tmp['tilt_y'][iz,iy,ix]
    tilt_z[i,:,:] = tmp['tilt_z'][iz,iy,ix]
    tilt_h[i,:,:] = tmp['tilt_h'][iz,iy,ix]
    tilt_sw[i,:,:] = tmp['tilt_sw'][iz,iy,ix]
    tilt_cw[i,:,:] = tmp['tilt_cw'][iz,iy,ix]
    dbfile.close()
    
    dbfile = open(ip+f"plan{t}_stretch.pkl", 'rb')
    tmp = pickle.load(dbfile)
    stretch_x[i,:,:] = tmp['stretch_x'][iz,iy,ix]
    stretch_y[i,:,:] = tmp['stretch_y'][iz,iy,ix]
    stretch_z[i,:,:] = tmp['stretch_z'][iz,iy,ix]
    stretch_h[i,:,:] = tmp['stretch_h'][iz,iy,ix]
    stretch_sw[i,:,:] = tmp['stretch_sw'][iz,iy,ix]
    stretch_cw[i,:,:] = tmp['stretch_cw'][iz,iy,ix]
    dbfile.close()
    
    dbfile = open(ip+f"plan{t}_baroclinic.pkl", 'rb')
    tmp = pickle.load(dbfile)
    bcl_x[i,:,:] = tmp['bcl_x'][iz,iy,ix]
    bcl_y[i,:,:] = tmp['bcl_y'][iz,iy,ix]
    bcl_h[i,:,:] = tmp['bcl_h'][iz,iy,ix]
    bcl_sw[i,:,:] = tmp['bcl_sw'][iz,iy,ix]
    bcl_cw[i,:,:] = tmp['bcl_cw'][iz,iy,ix]
    dbfile.close()

tilt_x1 = np.mean(tilt_x, axis=0)
tilt_y1 = np.mean(tilt_y, axis=0)
tilt_z1 = np.mean(tilt_z, axis=0)
tilt_h1 = np.mean(tilt_h, axis=0)
tilt_sw1 = np.mean(tilt_sw, axis=0)
tilt_cw1 = np.mean(tilt_cw, axis=0)
stretch_x1 = np.mean(stretch_x, axis=0)
stretch_y1 = np.mean(stretch_y, axis=0)
stretch_z1 = np.mean(stretch_z, axis=0)
stretch_h1 = np.mean(stretch_h, axis=0)
stretch_sw1 = np.mean(stretch_sw, axis=0)
stretch_cw1 = np.mean(stretch_cw, axis=0)
bcl_x1 = np.mean(bcl_x, axis=0)
bcl_y1 = np.mean(bcl_y, axis=0)
bcl_h1 = np.mean(bcl_h, axis=0)
bcl_sw1 = np.mean(bcl_sw, axis=0)
bcl_cw1 = np.mean(bcl_cw, axis=0)
del tilt_x,tilt_y,tilt_z,tilt_h,tilt_sw,tilt_cw
del stretch_x,stretch_y,stretch_z,stretch_h,stretch_sw,stretch_cw
del bcl_x,bcl_y,bcl_h,bcl_sw,bcl_cw



### Time 2 ###
times2 = np.arange(214, 220) # [217, 222] for 225 min downdraft
tilt_x = np.zeros(shape=(len(times2),17,17), dtype=float)
tilt_y = np.zeros(shape=(len(times2),17,17), dtype=float)
tilt_z = np.zeros(shape=(len(times2),17,17), dtype=float)
tilt_h = np.zeros(shape=(len(times2),17,17), dtype=float)
tilt_sw = np.zeros(shape=(len(times2),17,17), dtype=float)
tilt_cw = np.zeros(shape=(len(times2),17,17), dtype=float)
stretch_x = np.zeros(shape=(len(times2),17,17), dtype=float)
stretch_y = np.zeros(shape=(len(times2),17,17), dtype=float)
stretch_z = np.zeros(shape=(len(times2),17,17), dtype=float)
stretch_h = np.zeros(shape=(len(times2),17,17), dtype=float)
stretch_sw = np.zeros(shape=(len(times2),17,17), dtype=float)
stretch_cw = np.zeros(shape=(len(times2),17,17), dtype=float)
bcl_x = np.zeros(shape=(len(times2),17,17), dtype=float)
bcl_y = np.zeros(shape=(len(times2),17,17), dtype=float)
bcl_h = np.zeros(shape=(len(times2),17,17), dtype=float)
bcl_sw = np.zeros(shape=(len(times2),17,17), dtype=float)
bcl_cw = np.zeros(shape=(len(times2),17,17), dtype=float)

for i in range(len(times2)):
    t = times2[i]
    it = np.where(ptime/60 == t)[0][0]
    iz = np.abs(zh2 - z2_median[it]).argmin()
    i1 = np.abs(xh2 - x2_median[it]).argmin()
    i2 = np.abs(yh2 - y2_median[it]).argmin()
    ix = slice(i1-8,i1+9)
    iy = slice(i2-8,i2+9)
    
    dbfile = open(ip+f"plan{t}_tilt.pkl", 'rb')
    tmp = pickle.load(dbfile)
    tilt_x[i,:,:] = tmp['tilt_x'][iz,iy,ix]
    tilt_y[i,:,:] = tmp['tilt_y'][iz,iy,ix]
    tilt_z[i,:,:] = tmp['tilt_z'][iz,iy,ix]
    tilt_h[i,:,:] = tmp['tilt_h'][iz,iy,ix]
    tilt_sw[i,:,:] = tmp['tilt_sw'][iz,iy,ix]
    tilt_cw[i,:,:] = tmp['tilt_cw'][iz,iy,ix]
    dbfile.close()
    
    dbfile = open(ip+f"plan{t}_stretch.pkl", 'rb')
    tmp = pickle.load(dbfile)
    stretch_x[i,:,:] = tmp['stretch_x'][iz,iy,ix]
    stretch_y[i,:,:] = tmp['stretch_y'][iz,iy,ix]
    stretch_z[i,:,:] = tmp['stretch_z'][iz,iy,ix]
    stretch_h[i,:,:] = tmp['stretch_h'][iz,iy,ix]
    stretch_sw[i,:,:] = tmp['stretch_sw'][iz,iy,ix]
    stretch_cw[i,:,:] = tmp['stretch_cw'][iz,iy,ix]
    dbfile.close()
    
    dbfile = open(ip+f"plan{t}_baroclinic.pkl", 'rb')
    tmp = pickle.load(dbfile)
    bcl_x[i,:,:] = tmp['bcl_x'][iz,iy,ix]
    bcl_y[i,:,:] = tmp['bcl_y'][iz,iy,ix]
    bcl_h[i,:,:] = tmp['bcl_h'][iz,iy,ix]
    bcl_sw[i,:,:] = tmp['bcl_sw'][iz,iy,ix]
    bcl_cw[i,:,:] = tmp['bcl_cw'][iz,iy,ix]
    dbfile.close()

tilt_x2 = np.mean(tilt_x, axis=0)
tilt_y2 = np.mean(tilt_y, axis=0)
tilt_z2 = np.mean(tilt_z, axis=0)
tilt_h2 = np.mean(tilt_h, axis=0)
tilt_sw2 = np.mean(tilt_sw, axis=0)
tilt_cw2 = np.mean(tilt_cw, axis=0)
stretch_x2 = np.mean(stretch_x, axis=0)
stretch_y2 = np.mean(stretch_y, axis=0)
stretch_z2 = np.mean(stretch_z, axis=0)
stretch_h2 = np.mean(stretch_h, axis=0)
stretch_sw2 = np.mean(stretch_sw, axis=0)
stretch_cw2 = np.mean(stretch_cw, axis=0)
bcl_x2 = np.mean(bcl_x, axis=0)
bcl_y2 = np.mean(bcl_y, axis=0)
bcl_h2 = np.mean(bcl_h, axis=0)
bcl_sw2 = np.mean(bcl_sw, axis=0)
bcl_cw2 = np.mean(bcl_cw, axis=0)
del tilt_x,tilt_y,tilt_z,tilt_h,tilt_sw,tilt_cw
del stretch_x,stretch_y,stretch_z,stretch_h,stretch_sw,stretch_cw
del bcl_x,bcl_y,bcl_h,bcl_sw,bcl_cw



### Time 3 ###
times3 = np.arange(218, 220) # [220, 225] for 225 min rotor
tilt_x = np.zeros(shape=(len(times3),17,17), dtype=float)
tilt_y = np.zeros(shape=(len(times3),17,17), dtype=float)
tilt_z = np.zeros(shape=(len(times3),17,17), dtype=float)
tilt_h = np.zeros(shape=(len(times3),17,17), dtype=float)
tilt_sw = np.zeros(shape=(len(times3),17,17), dtype=float)
tilt_cw = np.zeros(shape=(len(times3),17,17), dtype=float)
stretch_x = np.zeros(shape=(len(times3),17,17), dtype=float)
stretch_y = np.zeros(shape=(len(times3),17,17), dtype=float)
stretch_z = np.zeros(shape=(len(times3),17,17), dtype=float)
stretch_h = np.zeros(shape=(len(times3),17,17), dtype=float)
stretch_sw = np.zeros(shape=(len(times3),17,17), dtype=float)
stretch_cw = np.zeros(shape=(len(times3),17,17), dtype=float)
bcl_x = np.zeros(shape=(len(times3),17,17), dtype=float)
bcl_y = np.zeros(shape=(len(times3),17,17), dtype=float)
bcl_h = np.zeros(shape=(len(times3),17,17), dtype=float)
bcl_sw = np.zeros(shape=(len(times3),17,17), dtype=float)
bcl_cw = np.zeros(shape=(len(times3),17,17), dtype=float)

for i in range(len(times3)):
    t = times3[i]
    it = np.where(ptime/60 == t)[0][0]
    iz = np.abs(zh2 - z2_median[it]).argmin()
    i1 = np.abs(xh2 - x2_median[it]).argmin()
    i2 = np.abs(yh2 - y2_median[it]).argmin()
    ix = slice(i1-8,i1+9)
    iy = slice(i2-8,i2+9)
    
    dbfile = open(ip+f"plan{t}_tilt.pkl", 'rb')
    tmp = pickle.load(dbfile)
    tilt_x[i,:,:] = tmp['tilt_x'][iz,iy,ix]
    tilt_y[i,:,:] = tmp['tilt_y'][iz,iy,ix]
    tilt_z[i,:,:] = tmp['tilt_z'][iz,iy,ix]
    tilt_h[i,:,:] = tmp['tilt_h'][iz,iy,ix]
    tilt_sw[i,:,:] = tmp['tilt_sw'][iz,iy,ix]
    tilt_cw[i,:,:] = tmp['tilt_cw'][iz,iy,ix]
    dbfile.close()
    
    dbfile = open(ip+f"plan{t}_stretch.pkl", 'rb')
    tmp = pickle.load(dbfile)
    stretch_x[i,:,:] = tmp['stretch_x'][iz,iy,ix]
    stretch_y[i,:,:] = tmp['stretch_y'][iz,iy,ix]
    stretch_z[i,:,:] = tmp['stretch_z'][iz,iy,ix]
    stretch_h[i,:,:] = tmp['stretch_h'][iz,iy,ix]
    stretch_sw[i,:,:] = tmp['stretch_sw'][iz,iy,ix]
    stretch_cw[i,:,:] = tmp['stretch_cw'][iz,iy,ix]
    dbfile.close()
    
    dbfile = open(ip+f"plan{t}_baroclinic.pkl", 'rb')
    tmp = pickle.load(dbfile)
    bcl_x[i,:,:] = tmp['bcl_x'][iz,iy,ix]
    bcl_y[i,:,:] = tmp['bcl_y'][iz,iy,ix]
    bcl_h[i,:,:] = tmp['bcl_h'][iz,iy,ix]
    bcl_sw[i,:,:] = tmp['bcl_sw'][iz,iy,ix]
    bcl_cw[i,:,:] = tmp['bcl_cw'][iz,iy,ix]
    dbfile.close()

tilt_x3 = np.mean(tilt_x, axis=0)
tilt_y3 = np.mean(tilt_y, axis=0)
tilt_z3 = np.mean(tilt_z, axis=0)
tilt_h3 = np.mean(tilt_h, axis=0)
tilt_sw3 = np.mean(tilt_sw, axis=0)
tilt_cw3 = np.mean(tilt_cw, axis=0)
stretch_x3 = np.mean(stretch_x, axis=0)
stretch_y3 = np.mean(stretch_y, axis=0)
stretch_z3 = np.mean(stretch_z, axis=0)
stretch_h3 = np.mean(stretch_h, axis=0)
stretch_sw3 = np.mean(stretch_sw, axis=0)
stretch_cw3 = np.mean(stretch_cw, axis=0)
bcl_x3 = np.mean(bcl_x, axis=0)
bcl_y3 = np.mean(bcl_y, axis=0)
bcl_h3 = np.mean(bcl_h, axis=0)
bcl_sw3 = np.mean(bcl_sw, axis=0)
bcl_cw3 = np.mean(bcl_cw, axis=0)
del tilt_x,tilt_y,tilt_z,tilt_h,tilt_sw,tilt_cw
del stretch_x,stretch_y,stretch_z,stretch_h,stretch_sw,stretch_cw
del bcl_x,bcl_y,bcl_h,bcl_sw,bcl_cw



itp1 = np.where(ptime/60 == times1[-1])[0][0]
itp2 = np.where(ptime/60 == times2[-1])[0][0]
itp3 = np.where(ptime/60 == times3[-1])[0][0]

xp1 = x1_ml[itp1,:]
yp1 = y1_ml[itp1,:]
zp1 = z1_ml[itp1,:]
xmp1 = x1_median[itp1]
ymp1 = y1_median[itp1]
zmp1 = z1_median[itp1]

xp2 = x2_ml[itp2,:]
yp2 = y2_ml[itp2,:]
zp2 = z2_ml[itp2,:]
xmp2 = x2_median[itp2]
ymp2 = y2_median[itp2]
zmp2 = z2_median[itp2]

xp3 = x2_ml[itp3,:]
yp3 = y2_ml[itp3,:]
zp3 = z2_ml[itp3,:]
xmp3 = x2_median[itp3]
ymp3 = y2_median[itp3]
zmp3 = z2_median[itp3]

iz1 = np.abs(zh1 - zmp1).argmin()
i1 = np.abs(xh1 - xmp1).argmin()
i2 = np.abs(yh1 - ymp1).argmin()
ix1 = slice(i1-8,i1+9)
iy1 = slice(i2-8,i2+9)

iz2 = np.abs(zh2 - zmp2).argmin()
i1 = np.abs(xh2 - xmp2).argmin()
i2 = np.abs(yh2 - ymp2).argmin()
ix2 = slice(i1-8,i1+9)
iy2 = slice(i2-8,i2+9)

iz3 = np.abs(zh2 - zmp3).argmin()
i1 = np.abs(xh2 - xmp3).argmin()
i2 = np.abs(yh2 - ymp3).argmin()
ix3 = slice(i1-8,i1+9)
iy3 = slice(i2-8,i2+9)


#%% Load data for composite vorticity tendency plots, all times - parcel-centered composites ***FOR PAPER FIGS***

ip = '/Users/morgan.schneider/Documents/merger/merger-125m/'

# fp = '/Volumes/Promise_Pegasus_70TB/merger/merger-125m/'
fp = '/Users/morgan.schneider/Documents/merger/merger-125m/temporary/'

# Load parcel time
# if 'ptime' not in locals():
#     ds = nc.Dataset(fp+'cm1out_pdata.nc')
#     ptime = ds.variables['time'][:].data
#     ds.close()

ptime = np.linspace(10800,14400,241)

# Load filtered parcel data
dbfile = open(ip+'traj_MV1.pkl', 'rb')
traj = pickle.load(dbfile)
dbfile.close()

# Load 210 min source regions and pull mid-level source
dbfile = open(ip+f"traj_clusters_210min_v2.pkl", 'rb')
tmp = pickle.load(dbfile)
cc = tmp['mv1']
dbfile.close()

x1_ml = traj[f"210min"]['x'][:,(cc==1)]
y1_ml = traj[f"210min"]['y'][:,(cc==1)]
z1_ml = traj[f"210min"]['z'][:,(cc==1)]

# Load 220 min source regions and pull mid-level source
dbfile = open(ip+f"traj_clusters_220min_v2.pkl", 'rb')
tmp = pickle.load(dbfile)
cc = tmp['mv1']
dbfile.close()

x2_ml = traj[f"220min"]['x'][:,(cc==1)]
y2_ml = traj[f"220min"]['y'][:,(cc==1)]
z2_ml = traj[f"220min"]['z'][:,(cc==1)]

x1_median = np.median(x1_ml, axis=1)
y1_median = np.median(y1_ml, axis=1)
z1_median = np.median(z1_ml, axis=1)
x2_median = np.median(x2_ml, axis=1)
y2_median = np.median(y2_ml, axis=1)
z2_median = np.median(z2_ml, axis=1)


# Load model grid
ds = nc.Dataset(fp+'cm1out_000043.nc')
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
zh = ds.variables['z'][:].data
ds.close()

# Load vorticity tendencies
dbfile = open(ip+f"vten_traj_210min_parcels.pkl", 'rb')
vten1 = pickle.load(dbfile)
stime1 = vten1['time']
dbfile.close()

dbfile = open(ip+f"vten_traj_220min_parcels.pkl", 'rb')
vten2 = pickle.load(dbfile)
stime2 = vten2['time']
dbfile.close()


# Choose parcel ending layer
dz = 150
it210 = np.where(ptime == 210*60)[0][0]
it220 = np.where(ptime == 220*60)[0][0]
iz1 = ((z1_ml[it210,:] >= z1_median[it210]-dz) & (z1_ml[it210,:] <= z1_median[it210]+dz))
iz2 = ((z2_ml[it220,:] >= z2_median[it220]-dz) & (z2_ml[it220,:] <= z2_median[it220]+dz))

### Choose averaging times ###
times1 = np.arange(203,209)
times2 = np.arange(214,219)
times3 = np.arange(218,220)


### Time 1 ###
it0 = np.where(stime1 == times1[0])[0][0]
itf = np.where(stime1 == times1[-1])[0][0]
it1 = slice(it0,itf)

stretch_x1 = np.mean(vten1['stretch_x'][iz1,it1,:,:], axis=(0,1))
stretch_y1 = np.mean(vten1['stretch_y'][iz1,it1,:,:], axis=(0,1))
stretch_z1 = np.mean(vten1['stretch_z'][iz1,it1,:,:], axis=(0,1))
stretch_h1 = np.mean(vten1['stretch_h'][iz1,it1,:,:], axis=(0,1))
stretch_sw1 = np.mean(vten1['stretch_sw'][iz1,it1,:,:], axis=(0,1))
stretch_cw1 = np.mean(vten1['stretch_cw'][iz1,it1,:,:], axis=(0,1))

tilt_x1 = np.mean(vten1['tilt_x'][iz1,it1,:,:], axis=(0,1))
tilt_y1 = np.mean(vten1['tilt_y'][iz1,it1,:,:], axis=(0,1))
tilt_z1 = np.mean(vten1['tilt_z'][iz1,it1,:,:], axis=(0,1))
tilt_h1 = np.mean(vten1['tilt_h'][iz1,it1,:,:], axis=(0,1))
tilt_sw1 = np.mean(vten1['tilt_sw'][iz1,it1,:,:], axis=(0,1))
tilt_cw1 = np.mean(vten1['tilt_cw'][iz1,it1,:,:], axis=(0,1))

bcl_x1 = np.mean(vten1['bcl_x'][iz1,it1,:,:], axis=(0,1))
bcl_y1 = np.mean(vten1['bcl_y'][iz1,it1,:,:], axis=(0,1))
bcl_h1 = np.mean(vten1['bcl_h'][iz1,it1,:,:], axis=(0,1))
bcl_sw1 = np.mean(vten1['bcl_sw'][iz1,it1,:,:], axis=(0,1))
bcl_cw1 = np.mean(vten1['bcl_cw'][iz1,it1,:,:], axis=(0,1))


### Time 2 ###
it0 = np.where(stime2 == times2[0])[0][0]
itf = np.where(stime2 == times2[-1])[0][0]
it2 = slice(it0,itf)

stretch_x2 = np.mean(vten2['stretch_x'][iz2,it2,:,:], axis=(0,1))
stretch_y2 = np.mean(vten2['stretch_y'][iz2,it2,:,:], axis=(0,1))
stretch_z2 = np.mean(vten2['stretch_z'][iz2,it2,:,:], axis=(0,1))
stretch_h2 = np.mean(vten2['stretch_h'][iz2,it2,:,:], axis=(0,1))
stretch_sw2 = np.mean(vten2['stretch_sw'][iz2,it2,:,:], axis=(0,1))
stretch_cw2 = np.mean(vten2['stretch_cw'][iz2,it2,:,:], axis=(0,1))

tilt_x2 = np.mean(vten2['tilt_x'][iz2,it2,:,:], axis=(0,1))
tilt_y2 = np.mean(vten2['tilt_y'][iz2,it2,:,:], axis=(0,1))
tilt_z2 = np.mean(vten2['tilt_z'][iz2,it2,:,:], axis=(0,1))
tilt_h2 = np.mean(vten2['tilt_h'][iz2,it2,:,:], axis=(0,1))
tilt_sw2 = np.mean(vten2['tilt_sw'][iz2,it2,:,:], axis=(0,1))
tilt_cw2 = np.mean(vten2['tilt_cw'][iz2,it2,:,:], axis=(0,1))

bcl_x2 = np.mean(vten2['bcl_x'][iz2,it2,:,:], axis=(0,1))
bcl_y2 = np.mean(vten2['bcl_y'][iz2,it2,:,:], axis=(0,1))
bcl_h2 = np.mean(vten2['bcl_h'][iz2,it2,:,:], axis=(0,1))
bcl_sw2 = np.mean(vten2['bcl_sw'][iz2,it2,:,:], axis=(0,1))
bcl_cw2 = np.mean(vten2['bcl_cw'][iz2,it2,:,:], axis=(0,1))


### Time 3 ###
it0 = np.where(stime2 == times3[0])[0][0]
itf = np.where(stime2 == times3[-1])[0][0]
it3 = slice(it0,itf)

stretch_x3 = np.mean(vten2['stretch_x'][iz2,it3,:,:], axis=(0,1))
stretch_y3 = np.mean(vten2['stretch_y'][iz2,it3,:,:], axis=(0,1))
stretch_z3 = np.mean(vten2['stretch_z'][iz2,it3,:,:], axis=(0,1))
stretch_h3 = np.mean(vten2['stretch_h'][iz2,it3,:,:], axis=(0,1))
stretch_sw3 = np.mean(vten2['stretch_sw'][iz2,it3,:,:], axis=(0,1))
stretch_cw3 = np.mean(vten2['stretch_cw'][iz2,it3,:,:], axis=(0,1))

tilt_x3 = np.mean(vten2['tilt_x'][iz2,it3,:,:], axis=(0,1))
tilt_y3 = np.mean(vten2['tilt_y'][iz2,it3,:,:], axis=(0,1))
tilt_z3 = np.mean(vten2['tilt_z'][iz2,it3,:,:], axis=(0,1))
tilt_h3 = np.mean(vten2['tilt_h'][iz2,it3,:,:], axis=(0,1))
tilt_sw3 = np.mean(vten2['tilt_sw'][iz2,it3,:,:], axis=(0,1))
tilt_cw3 = np.mean(vten2['tilt_cw'][iz2,it3,:,:], axis=(0,1))

bcl_x3 = np.mean(vten2['bcl_x'][iz2,it3,:,:], axis=(0,1))
bcl_y3 = np.mean(vten2['bcl_y'][iz2,it3,:,:], axis=(0,1))
bcl_h3 = np.mean(vten2['bcl_h'][iz2,it3,:,:], axis=(0,1))
bcl_sw3 = np.mean(vten2['bcl_sw'][iz2,it3,:,:], axis=(0,1))
bcl_cw3 = np.mean(vten2['bcl_cw'][iz2,it3,:,:], axis=(0,1))


itp1 = np.where(ptime/60 == times1[-1])[0][0]
xp1 = x1_ml[itp1,:]
yp1 = y1_ml[itp1,:]
ixp = np.abs(xh - x1_median[itp1]).argmin()
iyp = np.abs(yh - y1_median[itp1]).argmin()
ix1 = slice(ixp-8,ixp+9)
iy1 = slice(iyp-8,iyp+9)

itp2 = np.where(ptime/60 == times2[-1])[0][0]
xp2 = x2_ml[itp2,:]
yp2 = y2_ml[itp2,:]
ixp = np.abs(xh - x2_median[itp2]).argmin()
iyp = np.abs(yh - y2_median[itp2]).argmin()
ix2 = slice(ixp-8,ixp+9)
iy2 = slice(iyp-8,iyp+9)

itp3 = np.where(ptime/60 == times3[-1])[0][0]
xp3 = x2_ml[itp3,:]
yp3 = y2_ml[itp3,:]
ixp = np.abs(xh - x2_median[itp3]).argmin()
iyp = np.abs(yh - y2_median[itp3]).argmin()
ix3 = slice(ixp-8,ixp+9)
iy3 = slice(iyp-8,iyp+9)


#%% Make the tendency plan view plots ***PAPER FIG***

from matplotlib.ticker import MultipleLocator


figsave = False

# if using parcel-centered composites
# xh1 = xh; yh1 = yh; xh2 = xh; yh2 = yh

# xl1 = [xh[ix1][0], xh[ix1][-1]]
# yl1 = [yh[iy1][0], yh[iy1][-1]]
# xl2 = [xh[ix2][0], xh[ix2][-1]]
# yl2 = [yh[iy2][0], yh[iy2][-1]]
# xl3 = [xh[ix3][0], xh[ix3][-1]]
# yl3 = [yh[iy3][0], yh[iy3][-1]]
# x1 = xh[ix1]; y1 = yh[iy1]
# x2 = xh[ix2]; y2 = yh[iy2]
# x3 = xh[ix3]; y3 = yh[iy3]

xl = [-2, 2]; yl = [-2, 2]
xx = np.linspace(xl[0], xl[1], 33); yy = np.linspace(yl[0], yl[1], 33)
xp = 0; yp = 0

xl1 = xl; xl2 = xl; xl3 = xl; yl1 = yl; yl2 = yl; yl3 = yl
xp1 = xp; xp2 = xp; xp3 = xp; yp1 = yp; yp2 = yp; yp3 = yp
x1 = xx; x2 = xx; x3 = xx; y1 = yy; y2 = yy; y3 = yy



lims = [-3e-4, 3e-4]
levs = np.linspace(lims[0], lims[1], 41)
levs = np.append(np.append([-0.01], levs), [0.01])

lims1 = lims; lims2 = lims; lims3 = lims
levs1 = levs; levs2 = levs; levs3 = levs

# cb_ticks = np.linspace()

# lims1 = [-4e-4, 4e-4]; levs1 = np.linspace(lims1[0], lims1[1], 41)
# lims2 = [-4e-4, 4e-4]; levs2 = np.linspace(lims2[0], lims2[1], 41)
# lims3 = [-4e-4, 4e-4]; levs3 = np.linspace(lims3[0], lims3[1], 41)
# levs1 = np.append(np.append([-0.01], levs1), [0.01])
# levs2 = np.append(np.append([-0.01], levs2), [0.01])
# levs3 = np.append(np.append([-0.01], levs3), [0.01])

cm = 'balance'

xl1 = [-1, 1]; xl2 = [-1, 1]; xl3 = [-1, 1]


### Vertical ###
fig,ax = plt.subplots(3, 2, figsize=(7,9), sharex=True, sharey=True, layout='constrained', subplot_kw=dict(box_aspect=1))

plot_contourf(x1, y1, tilt_z1, 'zvort', ax[0,0], levels=levs1, datalims=lims1, xlims=xl1, ylims=yl1, cmap=cm, cbar=False)
# ax[0,0].scatter(xp1, yp1, s=30, color='k', marker='.')
ax[0,0].scatter(xp1, yp1, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[0,0].set_title(f"                                       {times1[0]}-{times1[-1]} min", fontsize=16)
ax[0,0].set_ylabel('y (km)', fontsize=14)

c1 = plot_contourf(x1, y1, stretch_z1, 'zvort', ax[0,1], levels=levs1, datalims=lims1, xlims=xl1, ylims=yl1, cmap=cm, cbar=False)
cb1 = plt.colorbar(c1, ax=ax[0,1], extend='both')
cb1.set_label("d\u03B6/dt (s$^{-2}$)", fontsize=12)
# cb1.set_ticks(np.linspace(-3e-4, 3e-4, 7))
cb1.formatter.set_powerlimits((0,0))
# ax[0,1].scatter(xp1, yp1, s=30, color='k', marker='.')
ax[0,1].scatter(xp1, yp1, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)

plot_contourf(x2, y2, tilt_z2, 'zvort', ax[1,0], levels=levs2, datalims=lims2, xlims=xl2, ylims=yl2, cmap=cm, cbar=False)
# ax[1,0].scatter(xp2, yp2, s=30, color='k', marker='.')
ax[1,0].scatter(xp2, yp2, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[1,0].set_title(f"                                       {times2[0]}-{times2[-1]} min", fontsize=16)
ax[1,0].set_ylabel('y (km)', fontsize=14)

c2 = plot_contourf(x2, y2, stretch_z2, 'zvort', ax[1,1], levels=levs2, datalims=lims2, xlims=xl2, ylims=yl2, cmap=cm, cbar=False)
cb2 = plt.colorbar(c2, ax=ax[1,1], extend='both')
cb2.set_label("d\u03B6/dt (s$^{-2}$)", fontsize=12)
# cb2.set_ticks(np.linspace(-3e-4, 3e-4, 7))
cb2.formatter.set_powerlimits((0,0))
# ax[1,1].scatter(xp2, yp2, s=30, color='k', marker='.')
ax[1,1].scatter(xp2, yp2, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)

plot_contourf(x3, y3, tilt_z3, 'zvort', ax[2,0], levels=levs3, datalims=lims3, xlims=xl3, ylims=yl3, cmap=cm, cbar=False)
# ax[2,0].scatter(xp3, yp3, s=30, color='k', marker='.')
ax[2,0].scatter(xp3, yp3, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[2,0].set_title(f"                                       {times3[0]}-{times3[-1]} min", fontsize=16)
ax[2,0].set_xlabel('x (km)', fontsize=14)
ax[2,0].set_ylabel('y (km)', fontsize=14)

c3 = plot_contourf(x3, y3, stretch_z3, 'zvort', ax[2,1], levels=levs3, datalims=lims3, xlims=xl3, ylims=yl3, cmap=cm, cbar=False)
cb3 = plt.colorbar(c3, ax=ax[2,1], extend='both')
cb3.set_label("d\u03B6/dt (s$^{-2}$)", fontsize=12)
# cb3.set_ticks(np.linspace(-3e-4, 3e-4, 7))
cb3.formatter.set_powerlimits((0,0))
# ax[2,1].scatter(xp3, yp3, s=30, color='k', marker='.')
ax[2,1].scatter(xp3, yp3, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[2,1].set_xlabel('x (km)', fontsize=14)

for i in range(3):
    for j in range(2):
        ax[i,j].xaxis.set_major_locator(MultipleLocator(0.5))
        ax[i,j].xaxis.set_minor_locator(MultipleLocator(0.25))
        ax[i,j].yaxis.set_major_locator(MultipleLocator(0.5))
        ax[i,j].yaxis.set_minor_locator(MultipleLocator(0.25))
        ax[i,j].tick_params(axis='both', labelsize=12)

plt.suptitle(f"Composite \u03B6 tendency (parcel-centered)", fontsize=16)

if figsave:
    plt.savefig(ip+f"vort_tendency/zvort_composite_v2.png", dpi=300)




### Horizontal ###
fig,ax = plt.subplots(3, 3, figsize=(9.5,9), sharex=True, sharey=True, layout='constrained', subplot_kw=dict(box_aspect=1))

plot_contourf(x1, y1, tilt_h1, 'zvort', ax[0,0], levels=levs1, datalims=lims1, xlims=xl1, ylims=yl1, cmap=cm, cbar=False)
# ax[0,0].scatter(xp1, yp1, s=30, color='k', marker='.')
ax[0,0].scatter(xp1, yp1, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[0,0].set_ylabel('y (km)', fontsize=14)

plot_contourf(x1, y1, stretch_h1, 'zvort', ax[0,1], levels=levs1, datalims=lims1, xlims=xl1, ylims=yl1, cmap=cm, cbar=False)
# ax[0,1].scatter(xp1, yp1, s=30, color='k', marker='.')
ax[0,1].scatter(xp1, yp1, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[0,1].set_title(f"{times1[0]}-{times1[-1]} min", fontsize=16)

c1 = plot_contourf(x1, y1, bcl_h1, 'zvort', ax[0,2], levels=levs1, datalims=lims1, xlims=xl1, ylims=yl1, cmap=cm, cbar=False)
cb1 = plt.colorbar(c1, ax=ax[0,2], extend='both')
cb1.set_label("d\u03c9$_H$/dt (s$^{-2}$)", fontsize=12)
# cb1.set_ticks(np.linspace(-2.5e-4, 2.5e-4, 5))
cb1.formatter.set_powerlimits((0,0))
# ax[0,2].scatter(xp1, yp1, s=30, color='k', marker='.')
ax[0,2].scatter(xp1, yp1, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)


plot_contourf(x2, y2, tilt_h2, 'zvort', ax[1,0], levels=levs2, datalims=lims2, xlims=xl2, ylims=yl2, cmap=cm, cbar=False)
# ax[1,0].scatter(xp2, yp2, s=30, color='k', marker='.')
ax[1,0].scatter(xp2, yp2, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[1,0].set_ylabel('y (km)', fontsize=14)

plot_contourf(x2, y2, stretch_h2, 'zvort', ax[1,1], levels=levs2, datalims=lims2, xlims=xl2, ylims=yl2, cmap=cm, cbar=False)
# ax[1,1].scatter(xp2, yp2, s=30, color='k', marker='.')
ax[1,1].scatter(xp2, yp2, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[1,1].set_title(f"{times2[0]}-{times2[-1]} min", fontsize=16)

c2 = plot_contourf(x2, y2, bcl_h2, 'zvort', ax[1,2], levels=levs2, datalims=lims2, xlims=xl2, ylims=yl2, cmap=cm, cbar=False)
cb2 = plt.colorbar(c2, ax=ax[1,2], extend='both')
cb2.set_label("d\u03c9$_H$/dt (s$^{-2}$)", fontsize=12)
# cb2.set_ticks(np.linspace(-2.5e-4, 2.5e-4, 5))
cb2.formatter.set_powerlimits((0,0))
# ax[1,2].scatter(xp2, yp2, s=30, color='k', marker='.')
ax[1,2].scatter(xp2, yp2, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)


plot_contourf(x3, y3, tilt_h3, 'zvort', ax[2,0], levels=levs3, datalims=lims3, xlims=xl3, ylims=yl3, cmap=cm, cbar=False)
# ax[2,0].scatter(xp3, yp3, s=30, color='k', marker='.')
ax[2,0].scatter(xp3, yp3, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[2,0].set_xlabel('x (km)', fontsize=14)
ax[2,0].set_ylabel('y (km)', fontsize=14)

plot_contourf(x3, y3, stretch_h3, 'zvort', ax[2,1], levels=levs3, datalims=lims3, xlims=xl3, ylims=yl3, cmap=cm, cbar=False)
# ax[2,1].scatter(xp3, yp3, s=30, color='k', marker='.')
ax[2,1].scatter(xp3, yp3, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[2,1].set_title(f"{times3[0]}-{times3[-1]} min", fontsize=16)
ax[2,1].set_xlabel('x (km)', fontsize=14)

c3 = plot_contourf(x3, y3, bcl_h3, 'zvort', ax[2,2], levels=levs3, datalims=lims3, xlims=xl3, ylims=yl3, cmap=cm, cbar=False)
cb3 = plt.colorbar(c3, ax=ax[2,2], extend='both')
cb3.set_label("d\u03c9$_H$/dt (s$^{-2}$)", fontsize=12)
# cb3.set_ticks(np.linspace(-2.5e-4, 2.5e-4, 5))
cb3.formatter.set_powerlimits((0,0))
# ax[2,2].scatter(xp3, yp3, s=30, color='k', marker='.')
ax[2,2].scatter(xp3, yp3, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[2,2].set_xlabel('x (km)', fontsize=14)

for i in range(3):
    for j in range(3):
        ax[i,j].xaxis.set_major_locator(MultipleLocator(0.5))
        ax[i,j].xaxis.set_minor_locator(MultipleLocator(0.25))
        ax[i,j].yaxis.set_major_locator(MultipleLocator(0.5))
        ax[i,j].yaxis.set_minor_locator(MultipleLocator(0.25))
        ax[i,j].tick_params(axis='both', labelsize=12)

plt.suptitle(f"Composite |\u03c9$_H$| tendency (parcel-centered)", fontsize=18)

if figsave:
    plt.savefig(ip+f"vort_tendency/hvort_composite_v2.png", dpi=300)





#%% Other parcel-centered plan views (x, y, sw, cw)

### X vorticity ###
fig,ax = plt.subplots(3, 3, figsize=(9.5,9), sharex=True, sharey=True, layout='constrained', subplot_kw=dict(box_aspect=1))

plot_contourf(x1, y1, tilt_x1, 'zvort', ax[0,0], levels=levs1, datalims=lims1, xlims=xl1, ylims=yl1, cmap=cm, cbar=False)
# ax[0,0].scatter(xp1, yp1, s=30, color='k', marker='.')
ax[0,0].scatter(xp1, yp1, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[0,0].set_ylabel('y (km)', fontsize=14)

plot_contourf(x1, y1, stretch_x1, 'zvort', ax[0,1], levels=levs1, datalims=lims1, xlims=xl1, ylims=yl1, cmap=cm, cbar=False)
# ax[0,1].scatter(xp1, yp1, s=30, color='k', marker='.')
ax[0,1].scatter(xp1, yp1, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[0,1].set_title(f"{times1[0]}-{times1[-1]} min", fontsize=16)

c1 = plot_contourf(x1, y1, bcl_x1, 'zvort', ax[0,2], levels=levs1, datalims=lims1, xlims=xl1, ylims=yl1, cmap=cm, cbar=False)
cb1 = plt.colorbar(c1, ax=ax[0,2], extend='both')
cb1.set_label("d\u03BE/dt (s$^{-2}$)", fontsize=12)
cb1.formatter.set_powerlimits((0,0))
# ax[0,2].scatter(xp1, yp1, s=30, color='k', marker='.')
ax[0,2].scatter(xp1, yp1, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)


plot_contourf(x2, y2, tilt_x2, 'zvort', ax[1,0], levels=levs2, datalims=lims2, xlims=xl2, ylims=yl2, cmap=cm, cbar=False)
# ax[1,0].scatter(xp2, yp2, s=30, color='k', marker='.')
ax[1,0].scatter(xp2, yp2, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[1,0].set_ylabel('y (km)', fontsize=14)

plot_contourf(x2, y2, stretch_x2, 'zvort', ax[1,1], levels=levs2, datalims=lims2, xlims=xl2, ylims=yl2, cmap=cm, cbar=False)
# ax[1,1].scatter(xp2, yp2, s=30, color='k', marker='.')
ax[1,1].scatter(xp2, yp2, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[1,1].set_title(f"{times2[0]}-{times2[-1]} min", fontsize=16)

c2 = plot_contourf(x2, y2, bcl_x2, 'zvort', ax[1,2], levels=levs2, datalims=lims2, xlims=xl2, ylims=yl2, cmap=cm, cbar=False)
cb2 = plt.colorbar(c2, ax=ax[1,2], extend='both')
cb2.set_label("d\u03BE/dt (s$^{-2}$)", fontsize=12)
cb2.formatter.set_powerlimits((0,0))
# ax[1,2].scatter(xp2, yp2, s=30, color='k', marker='.')
ax[1,2].scatter(xp2, yp2, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)


plot_contourf(x3, y3, tilt_x3, 'zvort', ax[2,0], levels=levs3, datalims=lims3, xlims=xl3, ylims=yl3, cmap=cm, cbar=False)
# ax[2,0].scatter(xp3, yp3, s=30, color='k', marker='.')
ax[2,0].scatter(xp3, yp3, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[2,0].set_xlabel('x (km)', fontsize=14)
ax[2,0].set_ylabel('y (km)', fontsize=14)

plot_contourf(x3, y3, stretch_x3, 'zvort', ax[2,1], levels=levs3, datalims=lims3, xlims=xl3, ylims=yl3, cmap=cm, cbar=False)
# ax[2,1].scatter(xp3, yp3, s=30, color='k', marker='.')
ax[2,1].scatter(xp3, yp3, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[2,1].set_title(f"{times3[0]}-{times3[-1]} min", fontsize=16)
ax[2,1].set_xlabel('x (km)', fontsize=14)

c3 = plot_contourf(x3, y3, bcl_x3, 'zvort', ax[2,2], levels=levs3, datalims=lims3, xlims=xl3, ylims=yl3, cmap=cm, cbar=False)
cb3 = plt.colorbar(c3, ax=ax[2,2], extend='both')
cb3.set_label("d\u03BE/dt (s$^{-2}$)", fontsize=12)
cb3.formatter.set_powerlimits((0,0))
# ax[2,2].scatter(xp3, yp3, s=30, color='k', marker='.')
ax[2,2].scatter(xp3, yp3, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[2,2].set_xlabel('x (km)', fontsize=14)

plt.suptitle(f"Composite \u03BE tendency (parcel-centered)", fontsize=18)

if figsave:
    plt.savefig(f"/Users/morgan.schneider/Documents/merger/merger-125m/vort_tendency/xvort_composite.png", dpi=300)




### Y vorticity ###
fig,ax = plt.subplots(3, 3, figsize=(9.5,9), sharex=True, sharey=True, layout='constrained', subplot_kw=dict(box_aspect=1))

plot_contourf(x1, y1, tilt_y1, 'zvort', ax[0,0], levels=levs1, datalims=lims1, xlims=xl1, ylims=yl1, cmap=cm, cbar=False)
# ax[0,0].scatter(xp1, yp1, s=30, color='k', marker='.')
ax[0,0].scatter(xp1, yp1, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[0,0].set_ylabel('y (km)', fontsize=14)

plot_contourf(x1, y1, stretch_y1, 'zvort', ax[0,1], levels=levs1, datalims=lims1, xlims=xl1, ylims=yl1, cmap=cm, cbar=False)
# ax[0,1].scatter(xp1, yp1, s=30, color='k', marker='.')
ax[0,1].scatter(xp1, yp1, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[0,1].set_title(f"{times1[0]}-{times1[-1]} min", fontsize=16)

c1 = plot_contourf(x1, y1, bcl_y1, 'zvort', ax[0,2], levels=levs1, datalims=lims1, xlims=xl1, ylims=yl1, cmap=cm, cbar=False)
cb1 = plt.colorbar(c1, ax=ax[0,2], extend='both')
cb1.set_label("d\u03B7/dt (s$^{-2}$)", fontsize=12)
cb1.formatter.set_powerlimits((0,0))
# ax[0,2].scatter(xp1, yp1, s=30, color='k', marker='.')
ax[0,2].scatter(xp1, yp1, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)


plot_contourf(x2, y2, tilt_y2, 'zvort', ax[1,0], levels=levs2, datalims=lims2, xlims=xl2, ylims=yl2, cmap=cm, cbar=False)
# ax[1,0].scatter(xp2, yp2, s=30, color='k', marker='.')
ax[1,0].scatter(xp2, yp2, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[1,0].set_ylabel('y (km)', fontsize=14)

plot_contourf(x2, y2, stretch_y2, 'zvort', ax[1,1], levels=levs2, datalims=lims2, xlims=xl2, ylims=yl2, cmap=cm, cbar=False)
# ax[1,1].scatter(xp2, yp2, s=30, color='k', marker='.')
ax[1,1].scatter(xp2, yp2, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[1,1].set_title(f"{times2[0]}-{times2[-1]} min", fontsize=16)

c2 = plot_contourf(x2, y2, bcl_y2, 'zvort', ax[1,2], levels=levs2, datalims=lims2, xlims=xl2, ylims=yl2, cmap=cm, cbar=False)
cb2 = plt.colorbar(c2, ax=ax[1,2], extend='both')
cb2.set_label("d\u03B7/dt (s$^{-2}$)", fontsize=12)
cb2.formatter.set_powerlimits((0,0))
# ax[1,2].scatter(xp2, yp2, s=30, color='k', marker='.')
ax[1,2].scatter(xp2, yp2, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)


plot_contourf(x3, y3, tilt_y3, 'zvort', ax[2,0], levels=levs3, datalims=lims3, xlims=xl3, ylims=yl3, cmap=cm, cbar=False)
# ax[2,0].scatter(xp3, yp3, s=30, color='k', marker='.')
ax[2,0].scatter(xp3, yp3, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[2,0].set_xlabel('x (km)', fontsize=14)
ax[2,0].set_ylabel('y (km)', fontsize=14)

plot_contourf(x3, y3, stretch_y3, 'zvort', ax[2,1], levels=levs3, datalims=lims3, xlims=xl3, ylims=yl3, cmap=cm, cbar=False)
# ax[2,1].scatter(xp3, yp3, s=30, color='k', marker='.')
ax[2,1].scatter(xp3, yp3, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[2,1].set_title(f"{times3[0]}-{times3[-1]} min", fontsize=16)
ax[2,1].set_xlabel('x (km)', fontsize=14)

c3 = plot_contourf(x3, y3, bcl_y3, 'zvort', ax[2,2], levels=levs3, datalims=lims3, xlims=xl3, ylims=yl3, cmap=cm, cbar=False)
cb3 = plt.colorbar(c3, ax=ax[2,2], extend='both')
cb3.set_label("d\u03B7/dt (s$^{-2}$)", fontsize=12)
cb3.formatter.set_powerlimits((0,0))
# ax[2,2].scatter(xp3, yp3, s=30, color='k', marker='.')
ax[2,2].scatter(xp3, yp3, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[2,2].set_xlabel('x (km)', fontsize=14)

plt.suptitle(f"Composite \u03B7 tendency (parcel-centered)", fontsize=18)

if figsave:
    plt.savefig(f"/Users/morgan.schneider/Documents/merger/merger-125m/vort_tendency/yvort_composite.png", dpi=300)





### Streamwise ###
fig,ax = plt.subplots(3, 3, figsize=(9.5,9), sharex=True, sharey=True, layout='constrained', subplot_kw=dict(box_aspect=1))

plot_contourf(x1, y1, tilt_sw1, 'zvort', ax[0,0], levels=levs1, datalims=lims1, xlims=xl1, ylims=yl1, cmap=cm, cbar=False)
# ax[0,0].scatter(xp1, yp1, s=30, color='k', marker='.')
ax[0,0].scatter(xp1, yp1, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[0,0].set_ylabel('y (km)', fontsize=14)

plot_contourf(x1, y1, stretch_sw1, 'zvort', ax[0,1], levels=levs1, datalims=lims1, xlims=xl1, ylims=yl1, cmap=cm, cbar=False)
# ax[0,1].scatter(xp1, yp1, s=30, color='k', marker='.')
ax[0,1].scatter(xp1, yp1, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[0,1].set_title(f"{times1[0]}-{times1[-1]} min", fontsize=16)

c1 = plot_contourf(x1, y1, bcl_sw1, 'zvort', ax[0,2], levels=levs1, datalims=lims1, xlims=xl1, ylims=yl1, cmap=cm, cbar=False)
cb1 = plt.colorbar(c1, ax=ax[0,2], extend='both')
cb1.set_label("d\u03c9$_{SW}$/dt (s$^{-2}$)", fontsize=12)
cb1.formatter.set_powerlimits((0,0))
# ax[0,2].scatter(xp1, yp1, s=30, color='k', marker='.')
ax[0,2].scatter(xp1, yp1, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)


plot_contourf(x2, y2, tilt_sw2, 'zvort', ax[1,0], levels=levs2, datalims=lims2, xlims=xl2, ylims=yl2, cmap=cm, cbar=False)
# ax[1,0].scatter(xp2, yp2, s=30, color='k', marker='.')
ax[1,0].scatter(xp2, yp2, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[1,0].set_ylabel('y (km)', fontsize=14)

plot_contourf(x2, y2, stretch_sw2, 'zvort', ax[1,1], levels=levs2, datalims=lims2, xlims=xl2, ylims=yl2, cmap=cm, cbar=False)
# ax[1,1].scatter(xp2, yp2, s=30, color='k', marker='.')
ax[1,1].scatter(xp2, yp2, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[1,1].set_title(f"{times2[0]}-{times2[-1]} min", fontsize=16)

c2 = plot_contourf(x2, y2, bcl_sw2, 'zvort', ax[1,2], levels=levs2, datalims=lims2, xlims=xl2, ylims=yl2, cmap=cm, cbar=False)
cb2 = plt.colorbar(c2, ax=ax[1,2], extend='both')
cb2.set_label("d\u03c9$_{SW}$/dt (s$^{-2}$)", fontsize=12)
cb2.formatter.set_powerlimits((0,0))
# ax[1,2].scatter(xp2, yp2, s=30, color='k', marker='.')
ax[1,2].scatter(xp2, yp2, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)


plot_contourf(x3, y3, tilt_sw3, 'zvort', ax[2,0], levels=levs3, datalims=lims3, xlims=xl3, ylims=yl3, cmap=cm, cbar=False)
# ax[2,0].scatter(xp3, yp3, s=30, color='k', marker='.')
ax[2,0].scatter(xp3, yp3, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[2,0].set_xlabel('x (km)', fontsize=14)
ax[2,0].set_ylabel('y (km)', fontsize=14)

plot_contourf(x3, y3, stretch_sw3, 'zvort', ax[2,1], levels=levs3, datalims=lims3, xlims=xl3, ylims=yl3, cmap=cm, cbar=False)
# ax[2,1].scatter(xp3, yp3, s=30, color='k', marker='.')
ax[2,1].scatter(xp3, yp3, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[2,1].set_title(f"{times3[0]}-{times3[-1]} min", fontsize=16)
ax[2,1].set_xlabel('x (km)', fontsize=14)

c3 = plot_contourf(x3, y3, bcl_sw3, 'zvort', ax[2,2], levels=levs3, datalims=lims3, xlims=xl3, ylims=yl3, cmap=cm, cbar=False)
cb3 = plt.colorbar(c3, ax=ax[2,2], extend='both')
cb3.set_label("d\u03c9$_{SW}$/dt (s$^{-2}$)", fontsize=12)
cb3.formatter.set_powerlimits((0,0))
# ax[2,2].scatter(xp3, yp3, s=30, color='k', marker='.')
ax[2,2].scatter(xp3, yp3, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[2,2].set_xlabel('x (km)', fontsize=14)

plt.suptitle(f"Composite streamwise \u03c9 tendency (parcel-centered)", fontsize=18)

if figsave:
    plt.savefig(f"/Users/morgan.schneider/Documents/merger/merger-125m/vort_tendency/swvort_composite.png", dpi=300)




### Crosswise ###
fig,ax = plt.subplots(3, 3, figsize=(9.5,9), sharex=True, sharey=True, layout='constrained', subplot_kw=dict(box_aspect=1))

plot_contourf(x1, y1, tilt_cw1, 'zvort', ax[0,0], levels=levs1, datalims=lims1, xlims=xl1, ylims=yl1, cmap=cm, cbar=False)
# ax[0,0].scatter(xp1, yp1, s=30, color='k', marker='.')
ax[0,0].scatter(xp1, yp1, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[0,0].set_ylabel('y (km)', fontsize=14)

plot_contourf(x1, y1, stretch_cw1, 'zvort', ax[0,1], levels=levs1, datalims=lims1, xlims=xl1, ylims=yl1, cmap=cm, cbar=False)
# ax[0,1].scatter(xp1, yp1, s=30, color='k', marker='.')
ax[0,1].scatter(xp1, yp1, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[0,1].set_title(f"{times1[0]}-{times1[-1]} min", fontsize=16)

c1 = plot_contourf(x1, y1, bcl_cw1, 'zvort', ax[0,2], levels=levs1, datalims=lims1, xlims=xl1, ylims=yl1, cmap=cm, cbar=False)
cb1 = plt.colorbar(c1, ax=ax[0,2], extend='both')
cb1.set_label("d\u03c9$_{CW}$/dt (s$^{-2}$)", fontsize=12)
cb1.formatter.set_powerlimits((0,0))
# ax[0,2].scatter(xp1, yp1, s=30, color='k', marker='.')
ax[0,2].scatter(xp1, yp1, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)


plot_contourf(x2, y2, tilt_cw2, 'zvort', ax[1,0], levels=levs2, datalims=lims2, xlims=xl2, ylims=yl2, cmap=cm, cbar=False)
# ax[1,0].scatter(xp2, yp2, s=30, color='k', marker='.')
ax[1,0].scatter(xp2, yp2, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[1,0].set_ylabel('y (km)', fontsize=14)

plot_contourf(x2, y2, stretch_cw2, 'zvort', ax[1,1], levels=levs2, datalims=lims2, xlims=xl2, ylims=yl2, cmap=cm, cbar=False)
# ax[1,1].scatter(xp2, yp2, s=30, color='k', marker='.')
ax[1,1].scatter(xp2, yp2, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[1,1].set_title(f"{times2[0]}-{times2[-1]} min", fontsize=16)

c2 = plot_contourf(x2, y2, bcl_cw2, 'zvort', ax[1,2], levels=levs2, datalims=lims2, xlims=xl2, ylims=yl2, cmap=cm, cbar=False)
cb2 = plt.colorbar(c2, ax=ax[1,2], extend='both')
cb2.set_label("d\u03c9$_{CW}$/dt (s$^{-2}$)", fontsize=12)
cb2.formatter.set_powerlimits((0,0))
# ax[1,2].scatter(xp2, yp2, s=30, color='k', marker='.')
ax[1,2].scatter(xp2, yp2, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)


plot_contourf(x3, y3, tilt_cw3, 'zvort', ax[2,0], levels=levs3, datalims=lims3, xlims=xl3, ylims=yl3, cmap=cm, cbar=False)
# ax[2,0].scatter(xp3, yp3, s=30, color='k', marker='.')
ax[2,0].scatter(xp3, yp3, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[2,0].set_xlabel('x (km)', fontsize=14)
ax[2,0].set_ylabel('y (km)', fontsize=14)

plot_contourf(x3, y3, stretch_cw3, 'zvort', ax[2,1], levels=levs3, datalims=lims3, xlims=xl3, ylims=yl3, cmap=cm, cbar=False)
# ax[2,1].scatter(xp3, yp3, s=30, color='k', marker='.')
ax[2,1].scatter(xp3, yp3, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[2,1].set_title(f"{times3[0]}-{times3[-1]} min", fontsize=16)
ax[2,1].set_xlabel('x (km)', fontsize=14)

c3 = plot_contourf(x3, y3, bcl_cw3, 'zvort', ax[2,2], levels=levs3, datalims=lims3, xlims=xl3, ylims=yl3, cmap=cm, cbar=False)
cb3 = plt.colorbar(c3, ax=ax[2,2], extend='both')
cb3.set_label("d\u03c9$_{CW}$/dt (s$^{-2}$)", fontsize=12)
cb3.formatter.set_powerlimits((0,0))
# ax[2,2].scatter(xp3, yp3, s=30, color='k', marker='.')
ax[2,2].scatter(xp3, yp3, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[2,2].set_xlabel('x (km)', fontsize=14)

plt.suptitle(f"Composite crosswise \u03c9 tendency (parcel-centered)", fontsize=18)

if figsave:
    plt.savefig(f"/Users/morgan.schneider/Documents/merger/merger-125m/vort_tendency/cwvort_composite.png", dpi=300)



#%% Plan views of my chosen times ***PAPER FIGS?***

fp = '/Users/mschne28/Documents/merger/merger-125m/'

time = 220

### Choose averaging times ###
times = np.arange(218,221)

figsave = False


# Load coordinate data
dbfile = open(fp+'coords.pkl', 'rb')
coords = pickle.load(dbfile)
xh = coords['xh']
yh = coords['yh']
zh = coords['zh']
ptime = coords['ptime']
dbfile.close()


# Load filtered parcel data
dbfile = open(fp+'traj_MV1.pkl', 'rb')
traj = pickle.load(dbfile)
dbfile.close()

# Load source regions and pull mid-level source
dbfile = open(fp+f"traj_clusters_{time}min_v2.pkl", 'rb')
tmp = pickle.load(dbfile)
cc = tmp['mv1']
dbfile.close()

z_ml = traj[f"{time}min"]['z'][:,(cc==1)]
z_median = np.median(z_ml, axis=1)
u_ml = traj[f"{time}min"]['u'][:,(cc==1)]
v_ml = traj[f"{time}min"]['v'][:,(cc==1)]


# Load vorticity tendencies
dbfile = open(fp+f"vten_traj_{time}min_parcels.pkl", 'rb')
vten = pickle.load(dbfile)
stime = vten['time']
dbfile.close()

# Choose parcel ending layer
dz = 1000
itmv = np.where(ptime == time*60)[0][0]
iz = ((z_ml[itmv,:] >= z_median[itmv]-dz) & (z_ml[itmv,:] <= z_median[itmv]+dz))







it0 = np.where(stime == times[0])[0][0]
itf = np.where(stime == times[-1])[0][0]
it = slice(it0,itf+1)


stretch_x = np.mean(vten['stretch_x'][iz,it,:,:], axis=(0,1))
stretch_y = np.mean(vten['stretch_y'][iz,it,:,:], axis=(0,1))
stretch_z = np.mean(vten['stretch_z'][iz,it,:,:], axis=(0,1))
stretch_h = np.mean(vten['stretch_h'][iz,it,:,:], axis=(0,1))
stretch_sw = np.mean(vten['stretch_sw'][iz,it,:,:], axis=(0,1))
stretch_cw = np.mean(vten['stretch_cw'][iz,it,:,:], axis=(0,1))

tilt_x = np.mean(vten['tilt_x'][iz,it,:,:], axis=(0,1))
tilt_y = np.mean(vten['tilt_y'][iz,it,:,:], axis=(0,1))
tilt_z = np.mean(vten['tilt_z'][iz,it,:,:], axis=(0,1))
tilt_h = np.mean(vten['tilt_h'][iz,it,:,:], axis=(0,1))
tilt_sw = np.mean(vten['tilt_sw'][iz,it,:,:], axis=(0,1))
tilt_cw = np.mean(vten['tilt_cw'][iz,it,:,:], axis=(0,1))

bcl_x = np.mean(vten['bcl_x'][iz,it,:,:], axis=(0,1))
bcl_y = np.mean(vten['bcl_y'][iz,it,:,:], axis=(0,1))
bcl_h = np.mean(vten['bcl_h'][iz,it,:,:], axis=(0,1))
bcl_sw = np.mean(vten['bcl_sw'][iz,it,:,:], axis=(0,1))
bcl_cw = np.mean(vten['bcl_cw'][iz,it,:,:], axis=(0,1))


# does each term contribute positively or negatively to crosswise vorticity in its current direction (AKA magnitude)
dbfile = open(fp+f"hvort_traj_{time}min.pkl", 'rb')
vort_traj = pickle.load(dbfile)
vort_x = vort_traj['xvort_ml']
vort_y = vort_traj['yvort_ml']
vort_sw = vort_traj['vort_sw_ml']
vort_cw = vort_traj['vort_cw_ml']
dbfile.close()

itp0 = np.where(ptime == times[0]*60)[0][0]
itpf = np.where(ptime == times[-1]*60)[0][0]
itp = slice(itp0,itpf+1)
vort_sw = vort_sw[itp,iz]
vort_cw = vort_cw[itp,iz]

cw_sign = vort_cw / np.abs(vort_cw)

scw = vten['stretch_cw'][iz,it,:,:]
tcw = vten['tilt_cw'][iz,it,:,:]
bcw = vten['bcl_cw'][iz,it,:,:]

for i in range(len(times)):
    for j in range(len(iz[(iz)])):
        scw[j,i,:,:] = scw[j,i,:,:] * cw_sign[i,j]
        tcw[j,i,:,:] = tcw[j,i,:,:] * cw_sign[i,j]
        bcw[j,i,:,:] = bcw[j,i,:,:] * cw_sign[i,j]

stretch_cw2 = np.mean(scw, axis=(0,1))
tilt_cw2 = np.mean(tcw, axis=(0,1))
bcl_cw2 = np.mean(bcw, axis=(0,1))

xvort = np.mean(np.median(vort_x[itp,iz], axis=0))
yvort = np.mean(np.median(vort_y[itp,iz], axis=0))

# Load storm motion
dbfile = open(fp+'storm_motion.pkl', 'rb')
sm = pickle.load(dbfile)
u_storm = sm['u_storm']
v_storm = sm['v_storm']
ws_storm = np.sqrt(u_storm**2 + v_storm**2)
dbfile.close()


u_sr = np.mean(np.median(u_ml[itp,iz], axis=0)) - np.mean(u_storm[times[0]-180:times[-1]-179])
v_sr = np.mean(np.median(v_ml[itp,iz], axis=0)) - np.mean(v_storm[times[0]-180:times[-1]-179])
ws_sr = np.mean(ws_storm[times[0]-180:times[-1]-179])
swvort = np.mean(np.median(vort_sw, axis=0))
cwvort = np.mean(np.median(vort_cw, axis=0))
swvort_x = u_sr/ws_sr * swvort
swvort_y = v_sr/ws_sr * swvort
cwvort_x = -v_sr/ws_sr * cwvort
cwvort_y = u_sr/ws_sr * cwvort

#% Make the plots

from matplotlib.ticker import MultipleLocator

from matplotlib.patches import FancyArrowPatch

def add_vectors(ax, x, y, dx, dy, lengthscale=10, *args, **kwargs):
    # SMALLER lengthscale for SHORTER arrows -- this is REVERSED from scale kwarg in quiver
    for i in range(len(x)):
        x1 = x[i]
        y1 = y[i]
        dx1 = dx[i]*lengthscale
        dy1 = dy[i]*lengthscale
        x2 = x1 + dx1
        y2 = y1 + dy1
        arrow = FancyArrowPatch((x1, y1), (x2, y2), *args, **kwargs)
        ax.add_patch(arrow)
        
    return arrow


# figsave = False


xx = np.linspace(-2, 2, 33)
yy = np.linspace(-2, 2, 33)
lims = [-2e-4, 2e-4]
levs = np.linspace(lims[0], lims[1], 41)
levs = np.append(np.append([-0.01], levs), [0.01])
xl = [-1, 1]
yl = [-1, 1]

cm = 'balance'


### Vertical ###
fig,ax = plt.subplots(1, 2, figsize=(7.5,3.5), sharex=True, sharey=True, layout='constrained', subplot_kw=dict(box_aspect=1))

plot_contourf(xx, yy, tilt_z, 'zvort', ax[0], levels=levs, datalims=lims, xlims=xl, ylims=yl, cmap=cm, cbar=False)
ax[0].scatter(0, 0, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[0].set_xlabel('x (km)', fontsize=14)
ax[0].set_ylabel('y (km)', fontsize=14)
ax[0].text(0.2, -0.9, 'Tilting', color='k', fontsize=18, fontweight='bold')
# ax[0].quiver(0, 0, xvort, yvort, color='k', scale=0.1, width=0.02, pivot='tail')
a_wind = add_vectors(ax[0], [0], [0], [u_sr], [v_sr],
            lengthscale=0.03, arrowstyle='simple', mutation_scale=15, ec='k', fc='k', lw=1)
a_vort = add_vectors(ax[0], [0], [0], [xvort], [yvort],
            lengthscale=20, arrowstyle='simple', mutation_scale=15, ec='k', fc='lightgray', lw=1)

c = plot_contourf(xx, yy, stretch_z, 'zvort', ax[1], levels=levs, datalims=lims, xlims=xl, ylims=yl, cmap=cm, cbar=False)
cb = plt.colorbar(c, ax=ax[1], extend='both')
cb.set_label("d\u03B6/dt (s$^{-2}$)", fontsize=12)
# cb.set_ticks(np.linspace(-3e-4, 3e-4, 7))
cb.formatter.set_powerlimits((0,0))
ax[1].scatter(0, 0, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[1].set_xlabel('x (km)', fontsize=14)
ax[1].text(-0.2, -0.9, 'Stretching', color='k', fontsize=18, fontweight='bold')
# ax[1].quiver(0, 0, xvort, yvort, color='k', scale=0.1, width=0.02, pivot='tail')
a_wind = add_vectors(ax[1], [0], [0], [u_sr], [v_sr],
            lengthscale=0.03, arrowstyle='simple', mutation_scale=15, ec='k', fc='k', lw=1)
a_vort = add_vectors(ax[1], [0], [0], [xvort], [yvort],
            lengthscale=20, arrowstyle='simple', mutation_scale=15, ec='k', fc='lightgray', lw=1)

for i in range(2):
    ax[i].xaxis.set_major_locator(MultipleLocator(0.5))
    ax[i].xaxis.set_minor_locator(MultipleLocator(0.25))
    ax[i].yaxis.set_major_locator(MultipleLocator(0.5))
    ax[i].yaxis.set_minor_locator(MultipleLocator(0.25))
    ax[i].tick_params(axis='both', labelsize=12)

# plt.suptitle(f" Composite \u03B6 tendency (parcel-centered) \n {times[0]}-{times[-1]} min ", fontsize=16)
plt.suptitle(f"{times[0]}-{times[-1]} min", fontsize=16)

if figsave:
    plt.savefig(fp+f"figs/zvort_composite_{times[0]}-{times[-1]}min.png", dpi=300)



### Horizontal ###
fig,ax = plt.subplots(1, 3, figsize=(9.5,3), sharex=True, sharey=True, layout='constrained', subplot_kw=dict(box_aspect=1))

plot_contourf(xx, yy, tilt_h, 'zvort', ax[0], levels=levs, datalims=lims, xlims=xl, ylims=yl, cmap=cm, cbar=False)
ax[0].scatter(0, 0, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[0].set_xlabel('x (km)', fontsize=14)
ax[0].set_ylabel('y (km)', fontsize=14)
ax[0].text(0.15, -0.9, 'Tilting', color='k', fontsize=18, fontweight='bold')
a_wind = add_vectors(ax[0], [0], [0], [u_sr], [v_sr],
            lengthscale=0.03, arrowstyle='simple', mutation_scale=15, ec='k', fc='k', lw=1)
a_vort = add_vectors(ax[0], [0], [0], [xvort], [yvort],
            lengthscale=20, arrowstyle='simple', mutation_scale=15, ec='k', fc='lightgray', lw=1)
# a_swvort = add_vectors(ax[0], [0], [0], [swvort_x], [swvort_y],
#             lengthscale=20, arrowstyle='simple', mutation_scale=12, ec='k', fc='magenta', lw=1)
# a_cwvort = add_vectors(ax[0], [0], [0], [cwvort_x], [cwvort_y],
#             lengthscale=20, arrowstyle='simple', mutation_scale=12, ec='k', fc='blue', lw=1)
ax[0].legend(handles=[a_vort,a_wind], labels=["\u03c9$_H$","V$_{SR}$"], loc='lower left', fontsize=10)

plot_contourf(xx, yy, stretch_h, 'zvort', ax[1], levels=levs, datalims=lims, xlims=xl, ylims=yl, cmap=cm, cbar=False)
ax[1].scatter(0, 0, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[1].set_title(f"{times[0]}-{times[-1]} min", fontsize=16)
ax[1].set_xlabel('x (km)', fontsize=14)
ax[1].text(-0.4, -0.9, 'Stretching', color='k', fontsize=18, fontweight='bold')
a_wind = add_vectors(ax[1], [0], [0], [u_sr], [v_sr],
            lengthscale=0.03, arrowstyle='simple', mutation_scale=15, ec='k', fc='k', lw=1)
a_vort = add_vectors(ax[1], [0], [0], [xvort], [yvort],
            lengthscale=20, arrowstyle='simple', mutation_scale=15, ec='k', fc='lightgray', lw=1)

c = plot_contourf(xx, yy, bcl_h, 'zvort', ax[2], levels=levs, datalims=lims, xlims=xl, ylims=yl, cmap=cm, cbar=False)
cb = plt.colorbar(c, ax=ax[2], extend='both')
cb.set_label("d\u03c9$_H$/dt (s$^{-2}$)", fontsize=12)
# cb.set_ticks(np.linspace(-2.5e-4, 2.5e-4, 5))
cb.formatter.set_powerlimits((0,0))
ax[2].scatter(0, 0, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[2].set_xlabel('x (km)', fontsize=14)
ax[2].text(-0.3, -0.9, 'Baroclinic', color='k', fontsize=18, fontweight='bold')
a_wind = add_vectors(ax[2], [0], [0], [u_sr], [v_sr],
            lengthscale=0.03, arrowstyle='simple', mutation_scale=15, ec='k', fc='k', lw=1)
a_vort = add_vectors(ax[2], [0], [0], [xvort], [yvort],
            lengthscale=20, arrowstyle='simple', mutation_scale=15, ec='k', fc='lightgray', lw=1)

for i in range(3):
    ax[i].xaxis.set_major_locator(MultipleLocator(0.5))
    ax[i].xaxis.set_minor_locator(MultipleLocator(0.25))
    ax[i].yaxis.set_major_locator(MultipleLocator(0.5))
    ax[i].yaxis.set_minor_locator(MultipleLocator(0.25))
    ax[i].tick_params(axis='both', labelsize=12)

# plt.suptitle(f" Composite |\u03c9$_H$| tendency (parcel-centered) \n {times[0]}-{times[-1]} min ", fontsize=16)
# plt.suptitle(f"{times[0]}-{times[-1]} min", fontsize=16)

if figsave:
    plt.savefig(fp+f"figs/hvort_composite_{times[0]}-{times[-1]}min.png", dpi=300)


#% Streamwise and crosswise plots

### Streamwise ###
fig,ax = plt.subplots(1, 3, figsize=(9.5,3), sharex=True, sharey=True, layout='constrained', subplot_kw=dict(box_aspect=1))

plot_contourf(xx, yy, tilt_sw, 'zvort', ax[0], levels=levs, datalims=lims, xlims=xl, ylims=yl, cmap=cm, cbar=False)
ax[0].scatter(0, 0, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[0].set_xlabel('x (km)', fontsize=14)
ax[0].set_ylabel('y (km)', fontsize=14)
ax[0].text(0.15, -0.9, 'Tilting', color='k', fontsize=18, fontweight='bold')
a_wind = add_vectors(ax[0], [0], [0], [u_sr], [v_sr],
            lengthscale=0.03, arrowstyle='simple', mutation_scale=15, ec='k', fc='k', lw=1)
a_vort = add_vectors(ax[0], [0], [0], [xvort], [yvort],
            lengthscale=20, arrowstyle='simple', mutation_scale=15, ec='k', fc='lightgray', lw=1)


plot_contourf(xx, yy, stretch_sw, 'zvort', ax[1], levels=levs, datalims=lims, xlims=xl, ylims=yl, cmap=cm, cbar=False)
ax[1].scatter(0, 0, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[1].set_title(f"{times[0]}-{times[-1]} min", fontsize=16)
ax[1].set_xlabel('x (km)', fontsize=14)
ax[1].text(-0.4, -0.9, 'Stretching', color='k', fontsize=18, fontweight='bold')
a_wind = add_vectors(ax[1], [0], [0], [u_sr], [v_sr],
            lengthscale=0.03, arrowstyle='simple', mutation_scale=15, ec='k', fc='k', lw=1)
a_vort = add_vectors(ax[1], [0], [0], [xvort], [yvort],
            lengthscale=20, arrowstyle='simple', mutation_scale=15, ec='k', fc='lightgray', lw=1)


c = plot_contourf(xx, yy, bcl_sw, 'zvort', ax[2], levels=levs, datalims=lims, xlims=xl, ylims=yl, cmap=cm, cbar=False)
cb = plt.colorbar(c, ax=ax[2], extend='both')
cb.set_label("d\u03c9$_{{sw}}$/dt (s$^{-2}$)", fontsize=12)
# cb.set_ticks(np.linspace(-2.5e-4, 2.5e-4, 5))
cb.formatter.set_powerlimits((0,0))
ax[2].scatter(0, 0, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[2].set_xlabel('x (km)', fontsize=14)
ax[2].text(-0.3, -0.9, 'Baroclinic', color='k', fontsize=18, fontweight='bold')
a_wind = add_vectors(ax[2], [0], [0], [u_sr], [v_sr],
            lengthscale=0.03, arrowstyle='simple', mutation_scale=15, ec='k', fc='k', lw=1)
a_vort = add_vectors(ax[2], [0], [0], [xvort], [yvort],
            lengthscale=20, arrowstyle='simple', mutation_scale=15, ec='k', fc='lightgray', lw=1)


for i in range(3):
    ax[i].xaxis.set_major_locator(MultipleLocator(0.5))
    ax[i].xaxis.set_minor_locator(MultipleLocator(0.25))
    ax[i].yaxis.set_major_locator(MultipleLocator(0.5))
    ax[i].yaxis.set_minor_locator(MultipleLocator(0.25))
    ax[i].tick_params(axis='both', labelsize=12)

# plt.suptitle(f" Composite \u03c9$_{{sw}}$ tendency (parcel-centered) \n {times[0]}-{times[-1]} min ", fontsize=16)

if figsave:
    plt.savefig(fp+f"figs/swvort_composite_{times[0]}-{times[-1]}min.png", dpi=300)



### Crosswise ###
fig,ax = plt.subplots(1, 3, figsize=(9.5,3), sharex=True, sharey=True, layout='constrained', subplot_kw=dict(box_aspect=1))

plot_contourf(xx, yy, tilt_cw, 'zvort', ax[0], levels=levs, datalims=lims, xlims=xl, ylims=yl, cmap=cm, cbar=False)
ax[0].scatter(0, 0, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[0].set_xlabel('x (km)', fontsize=14)
ax[0].set_ylabel('y (km)', fontsize=14)
ax[0].text(0.15, -0.9, 'Tilting', color='k', fontsize=18, fontweight='bold')
a_wind = add_vectors(ax[0], [0], [0], [u_sr], [v_sr],
            lengthscale=0.03, arrowstyle='simple', mutation_scale=15, ec='k', fc='k', lw=1)
a_vort = add_vectors(ax[0], [0], [0], [xvort], [yvort],
            lengthscale=20, arrowstyle='simple', mutation_scale=15, ec='k', fc='lightgray', lw=1)


plot_contourf(xx, yy, stretch_cw, 'zvort', ax[1], levels=levs, datalims=lims, xlims=xl, ylims=yl, cmap=cm, cbar=False)
ax[1].scatter(0, 0, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[1].set_title(f"{times[0]}-{times[-1]} min", fontsize=16)
ax[1].set_xlabel('x (km)', fontsize=14)
ax[1].text(-0.4, -0.9, 'Stretching', color='k', fontsize=18, fontweight='bold')
a_wind = add_vectors(ax[1], [0], [0], [u_sr], [v_sr],
            lengthscale=0.03, arrowstyle='simple', mutation_scale=15, ec='k', fc='k', lw=1)
a_vort = add_vectors(ax[1], [0], [0], [xvort], [yvort],
            lengthscale=20, arrowstyle='simple', mutation_scale=15, ec='k', fc='lightgray', lw=1)


c = plot_contourf(xx, yy, bcl_cw, 'zvort', ax[2], levels=levs, datalims=lims, xlims=xl, ylims=yl, cmap=cm, cbar=False)
cb = plt.colorbar(c, ax=ax[2], extend='both')
cb.set_label("d\u03c9$_{{sw}}$/dt (s$^{-2}$)", fontsize=12)
# cb.set_ticks(np.linspace(-2.5e-4, 2.5e-4, 5))
cb.formatter.set_powerlimits((0,0))
ax[2].scatter(0, 0, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[2].set_xlabel('x (km)', fontsize=14)
ax[2].text(-0.3, -0.9, 'Baroclinic', color='k', fontsize=18, fontweight='bold')
a_wind = add_vectors(ax[2], [0], [0], [u_sr], [v_sr],
            lengthscale=0.03, arrowstyle='simple', mutation_scale=15, ec='k', fc='k', lw=1)
a_vort = add_vectors(ax[2], [0], [0], [xvort], [yvort],
            lengthscale=20, arrowstyle='simple', mutation_scale=15, ec='k', fc='lightgray', lw=1)


for i in range(3):
    ax[i].xaxis.set_major_locator(MultipleLocator(0.5))
    ax[i].xaxis.set_minor_locator(MultipleLocator(0.25))
    ax[i].yaxis.set_major_locator(MultipleLocator(0.5))
    ax[i].yaxis.set_minor_locator(MultipleLocator(0.25))
    ax[i].tick_params(axis='both', labelsize=12)

# plt.suptitle(f" Composite \u03c9$_{{cw}}$ tendency (parcel-centered) \n {times[0]}-{times[-1]} min ", fontsize=16)

if figsave:
    plt.savefig(fp+f"figs/cwvort_composite_{times[0]}-{times[-1]}min.png", dpi=300)


#%% Composite tendency animations?

ip = '/Users/morgan.schneider/Documents/merger/merger-125m/cross_sections/MV1_vten/'


# Load parcel time
ds = nc.Dataset('/Volumes/Promise_Pegasus_70TB/merger/merger-125m/cm1out_pdata.nc')
ptime = ds.variables['time'][:].data
ds.close()

# Load filtered parcel data
dbfile = open('/Users/morgan.schneider/Documents/merger/traj_MV1.pkl', 'rb')
traj = pickle.load(dbfile)
dbfile.close()

# Load 210 min source regions and pull mid-level source
dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/traj_clusters/traj_clusters_210min_v2.pkl", 'rb')
tmp = pickle.load(dbfile)
cc = tmp['mv1']
dbfile.close()

x1_ml = traj[f"210min"]['x'][:,(cc==1)]/1000
y1_ml = traj[f"210min"]['y'][:,(cc==1)]/1000
z1_ml = traj[f"210min"]['z'][:,(cc==1)]/1000

# Load 220 min source regions and pull mid-level source
dbfile = open(f"/Users/morgan.schneider/Documents/merger/merger-125m/traj_clusters/traj_clusters_220min_v2.pkl", 'rb')
tmp = pickle.load(dbfile)
cc = tmp['mv1']
dbfile.close()

x2_ml = traj[f"220min"]['x'][:,(cc==1)]/1000
y2_ml = traj[f"220min"]['y'][:,(cc==1)]/1000
z2_ml = traj[f"220min"]['z'][:,(cc==1)]/1000

x1_median = np.median(x1_ml, axis=1)
y1_median = np.median(y1_ml, axis=1)
z1_median = np.median(z1_ml, axis=1)
x2_median = np.median(x2_ml, axis=1)
y2_median = np.median(y2_ml, axis=1)
z2_median = np.median(z2_ml, axis=1)


# Load model grid
ds = nc.Dataset('/Volumes/Promise_Pegasus_70TB/merger/merger-125m/cm1out_000014.nc')
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
zh = ds.variables['z'][:].data
ds.close()

# Load vorticity tendencies
dbfile = open(ip+f"vten_traj_210min.pkl", 'rb')
vten1 = pickle.load(dbfile)
stime1 = vten1['time']
dbfile.close()

dbfile = open(ip+f"vten_traj_220min.pkl", 'rb')
vten2 = pickle.load(dbfile)
stime2 = vten2['time']
dbfile.close()



### Choose averaging times ###
times1 = np.arange(195,211)
times2 = np.arange(205,221)


### Time 1 ###
it0 = np.where(stime1 == times1[0])[0][0]
itf = np.where(stime1 == times1[-1])[0][0]
it1 = slice(it0,itf)

stretch_z1 = vten1['stretch_z']
stretch_h1 = vten1['stretch_h']
tilt_z1 = vten1['tilt_z']
tilt_h1 = vten1['tilt_h']
bcl_h1 = vten1['bcl_h']


### Time 2 ###
it0 = np.where(stime2 == times2[0])[0][0]
itf = np.where(stime2 == times2[-1])[0][0]
it2 = slice(it0,itf)

stretch_z2 = vten2['stretch_z']
stretch_h2 = vten2['stretch_h']
tilt_z2 = vten2['tilt_z']
tilt_h2 = vten2['tilt_h']
bcl_h2 = vten2['bcl_h']


itp1 = np.where(ptime/60 == times1[-1])[0][0]
xp1 = x1_ml[itp1,:]
yp1 = y1_ml[itp1,:]
ixp = np.abs(xh - x1_median[itp1]).argmin()
iyp = np.abs(yh - y1_median[itp1]).argmin()
ix1 = slice(ixp-8,ixp+9)
iy1 = slice(iyp-8,iyp+9)

itp2 = np.where(ptime/60 == times2[-1])[0][0]
xp2 = x2_ml[itp2,:]
yp2 = y2_ml[itp2,:]
ixp = np.abs(xh - x2_median[itp2]).argmin()
iyp = np.abs(yh - y2_median[itp2]).argmin()
ix2 = slice(ixp-8,ixp+9)
iy2 = slice(iyp-8,iyp+9)




# if using parcel-centered composites
xh1 = xh; yh1 = yh; xh2 = xh; yh2 = yh

xl1 = [xh1[ix1][0], xh1[ix1][-1]] # [-19,-14] at 205-210 / [-11,-5] at 220-225
yl1 = [yh1[iy1][0], yh1[iy1][-1]] # [-72,-67]+0.75 at 205-210 / ?[-55,-49]+1 at 220-225
xl2 = [xh2[ix2][0], xh2[ix2][-1]]
yl2 = [yh2[iy2][0], yh2[iy2][-1]]


lims = [-1e-4, 1e-4]; levs = np.linspace(lims[0], lims[1], 41)
levs = np.append(np.append([-0.01], levs), [0.01])

cm = 'balance'


plot_time = 210

if plot_time == 210:
    x = xh1[ix1]
    y = yh1[iy1]
    tilt_z = tilt_z1
    stretch_z = stretch_z1
    tilt_h = tilt_h1
    stretch_h = stretch_h1
    bcl_h = bcl_h1
    xp = xp1
    yp = yp1
    xl = xl1
    yl = yl1
    times = times1
elif plot_time == 220:
    x = xh1[ix2]
    y = yh2[iy2]
    tilt_z = tilt_z2
    stretch_z = stretch_z2
    tilt_h = tilt_h2
    stretch_h = stretch_h2
    bcl_h = bcl_h2
    xp = xp2
    yp = yp2
    xl = xl2
    yl = yl2
    times = times2


# zvort tendency animation
if False:
    fig,ax = plt.subplots(1, 2, figsize=(9,4), sharey=True, subplot_kw=dict(box_aspect=1), layout='constrained')
    plot_contourf(x, y, tilt_z[0,:,:], 'zvort', ax[0], levels=levs, datalims=lims, xlims=xl, ylims=yl, cmap=cm, cbar=False)
    # ax[0,0].scatter(xp1, yp1, s=30, color='k', marker='.')
    ax[0].scatter(np.median(xp), np.median(yp), s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
    ax[0].set_title(f"                                                  {times[0]} min", fontsize=16)
    ax[0].set_ylabel('y (km)', fontsize=12)
    ax[0].set_xlabel('x (km)', fontsize=12)

    c1 = plot_contourf(x, y, stretch_z[0,:,:], 'zvort', ax[1], levels=levs, datalims=lims, xlims=xl, ylims=yl, cmap=cm, cbar=False)
    cb1 = plt.colorbar(c1, ax=ax[1], extend='both')
    cb1.set_label("d\u03B6/dt (s$^{-2}$)", fontsize=12)
    cb1.formatter.set_powerlimits((0,0))
    # ax[1].scatter(xp1, yp1, s=30, color='k', marker='.')
    ax[1].scatter(np.median(xp), np.median(yp), s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
    ax[1].set_xlabel('x (km)', fontsize=12)
    # plt.tight_layout()
    
    def animate_zten(i):
        global ax
        ax[0].clear()
        ax[1].clear()
        
        plot_contourf(x, y, tilt_z[i,:,:], 'zvort', ax[0], levels=levs, datalims=lims, xlims=xl, ylims=yl, cmap=cm, cbar=False)
        ax[0].scatter(np.median(xp), np.median(yp), s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
        ax[0].set_title(f"                                                  {times[i]} min", fontsize=16)
        ax[0].set_xlabel('x (km)', fontsize=12)
        ax[0].set_ylabel('y (km)', fontsize=12)
        
        plot_contourf(x, y, stretch_z[i,:,:], 'zvort', ax[1], levels=levs, datalims=lims, xlims=xl, ylims=yl, cmap=cm, cbar=False)
        ax[1].scatter(np.median(xp), np.median(yp), s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
        ax[1].set_xlabel('x (km)', fontsize=12)
        
    
    figsave = True
    
    anim = FuncAnimation(fig, animate_zten, frames=16, interval=500, repeat=False, blit=False)
    if figsave:
        anim.save(f"/Users/morgan.schneider/Documents/merger/merger-125m/vort_tendency/composite_zvort_{plot_time}min.gif", dpi=300)
    plt.show()



# hvort tendency animation
if True:
    fig,ax = plt.subplots(1, 3, figsize=(9.25,3), sharey=True, subplot_kw=dict(box_aspect=1), layout='constrained')
    plot_contourf(x, y, tilt_h[0,:,:], 'zvort', ax[0], levels=levs, datalims=lims, xlims=xl, ylims=yl, cmap=cm, cbar=False)
    # ax[0,0].scatter(xp1, yp1, s=30, color='k', marker='.')
    ax[0].scatter(np.median(xp), np.median(yp), s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
    ax[0].set_ylabel('y (km)', fontsize=12)
    ax[0].set_xlabel('x (km)', fontsize=12)
    
    plot_contourf(x, y, stretch_h[0,:,:], 'zvort', ax[1], levels=levs, datalims=lims, xlims=xl, ylims=yl, cmap=cm, cbar=False)
    # ax[1].scatter(xp1, yp1, s=30, color='k', marker='.')
    ax[1].scatter(np.median(xp), np.median(yp), s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
    ax[1].set_xlabel('x (km)', fontsize=12)
    ax[1].set_title(f"{times[0]} min", fontsize=16)

    c1 = plot_contourf(x, y, bcl_h[0,:,:], 'zvort', ax[2], levels=levs, datalims=lims, xlims=xl, ylims=yl, cmap=cm, cbar=False)
    cb1 = plt.colorbar(c1, ax=ax[2], extend='both')
    cb1.set_label("d|\u03c9$_H$|/dt (s$^{-2}$)", fontsize=12)
    cb1.formatter.set_powerlimits((0,0))
    # ax[1].scatter(xp1, yp1, s=30, color='k', marker='.')
    ax[2].scatter(np.median(xp), np.median(yp), s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
    ax[2].set_xlabel('x (km)', fontsize=12)
    # plt.tight_layout()
    
    def animate_hten(i):
        global ax
        ax[0].clear()
        ax[1].clear()
        ax[2].clear()
        
        plot_contourf(x, y, tilt_h[i,:,:], 'zvort', ax[0], levels=levs, datalims=lims, xlims=xl, ylims=yl, cmap=cm, cbar=False)
        ax[0].scatter(np.median(xp), np.median(yp), s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
        ax[0].set_xlabel('x (km)', fontsize=12)
        ax[0].set_ylabel('y (km)', fontsize=12)
        
        plot_contourf(x, y, stretch_h[i,:,:], 'zvort', ax[1], levels=levs, datalims=lims, xlims=xl, ylims=yl, cmap=cm, cbar=False)
        ax[1].scatter(np.median(xp), np.median(yp), s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
        ax[1].set_xlabel('x (km)', fontsize=12)
        ax[1].set_title(f"{times[i]} min", fontsize=16)
        
        plot_contourf(x, y, bcl_h[i,:,:], 'zvort', ax[2], levels=levs, datalims=lims, xlims=xl, ylims=yl, cmap=cm, cbar=False)
        ax[2].scatter(np.median(xp), np.median(yp), s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
        ax[2].set_xlabel('x (km)', fontsize=12)
        # plt.tight_layout()
    
    figsave = True
    
    anim = FuncAnimation(fig, animate_hten, frames=16, interval=500, repeat=False, blit=False)
    if figsave:
        anim.save(f"/Users/morgan.schneider/Documents/merger/merger-125m/vort_tendency/composite_hvort_{plot_time}min.gif", dpi=300)
    plt.show()

