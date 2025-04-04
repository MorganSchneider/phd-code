#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 12:27:13 2023

@author: morgan.schneider
"""

import plotly.graph_objects as go
from CM1utils import *
from RaxpolUtils import *
from scipy.ndimage import gaussian_filter


#%% Load CM1 data

fp = '/Volumes/Promise_Pegasus_70TB/2021-SupSims-final/R2_rerun/'
fn = 'cm1out_000079.nc'

ds = nc.Dataset(fp+'cm1out_000001.nc')
thr0 = ds.variables['th'][:].data * (1 + 0.61*ds.variables['qv'][:].data)
prs0 = ds.variables['prs'][:].data
ds.close()

ds = nc.Dataset(fp+fn)
xh = ds.variables['xh'][:].data
yh = ds.variables['yh'][:].data
z = ds.variables['z'][:].data
zvort = ds.variables['zvort'][:].data
# th = ds.variables['th'][:].data
# qv = ds.variables['qv'][:].data
winterp = ds.variables['winterp'][:].data
prs = ds.variables['prs'][:].data

thr = ds.variables['th'][:].data * (1 + 0.61*ds.variables['qv'][:].data - (
    ds.variables['qc'][:].data + ds.variables['qr'][:].data + ds.variables['qi'][:].data +
    ds.variables['qs'][:].data + ds.variables['qg'][:].data + ds.variables['qhl'][:].data))

thrpert = thr - thr0
prspert = prs - prs0
ds.close()

#%% Plot CM1 isosurfaces(zvort, prspert, thrpert, w)

ix = np.where((xh>=-2) & (xh<=3.0))[0]
iy = np.where((yh>=0.0) & (yh<=5.0))[0]
iz = np.where(z<=3.1)[0]

Z,Y,X = np.meshgrid(z[0:iz[-1]], yh[iy[0]:iy[-1]], xh[ix[0]:ix[-1]], indexing='ij')
vort = zvort[0, 0:iz[-1], iy[0]:iy[-1], ix[0]:ix[-1]]
prsp = prspert[0, 0:iz[-1], iy[0]:iy[-1], ix[0]:ix[-1]] / 100
thrp = gaussian_filter(thrpert[0, 0:iz[-1], iy[0]:iy[-1], ix[0]:ix[-1]], 1)
w = gaussian_filter(winterp[0, 0:iz[-1], iy[0]:iy[-1], ix[0]:ix[-1]], 2)

fig = go.Figure(data=[
    go.Isosurface(
    x=X.flatten(),
    y=Y.flatten(),
    z=Z.flatten(),
    value=vort.flatten(),
    colorscale='Reds',
    showscale=False,
    opacity=0.9,
    isomin=0.1,
    isomax=0.3,
    surface_count=3,
    # caps=dict(x_show=False, y_show=False),
    ),
    
    # go.Isosurface(
    # x=X.flatten(),
    # y=Y.flatten(),
    # z=Z.flatten(),
    # value=thrp.flatten(),
    # colorscale='plotly3',
    # showscale=False,
    # opacity=0.25,
    # isomin=-2,
    # isomax=100,
    # surface_count=2,
    # caps=dict(x_show=False, y_show=False),
    # ),
    
    # go.Isosurface(
    # x=X.flatten(),
    # y=Y.flatten(),
    # z=Z.flatten(),
    # value=w.flatten(),
    # colorscale='viridis_r',
    # showscale=False,
    # opacity=0.15,
    # isomin=10,
    # isomax=30,
    # surface_count=3,
    # caps=dict(x_show=False, y_show=False),
    # ),
    
    # go.Isosurface(
    # x=X.flatten(),
    # y=Y.flatten(),
    # z=Z.flatten(),
    # value=prsp.flatten(),
    # colorscale='Blues_r',
    # showscale=False,
    # opacity=0.9,
    # isomin=-10,
    # isomax=-10,
    # surface_count=1,
    # # caps=dict(x_show=False, y_show=False),
    # ),
    
    ])

fig.update_layout(
    margin=dict(t=0, l=0, b=0),
    scene_camera_eye=dict(x=-0.3, y=-2.2, z=0.2))
fig.show(renderer='browser')



#%% Load RaXPol data


dbfile= open('/Users/morgan.schneider/Documents/perils2023/iop2/circuit_data.pkl', 'rb')
tmp = pickle.load(dbfile)
vol = tmp['Rax']
dbfile.close()

vi = 7


#%%
x = vol[vi]['xx']
y = vol[vi]['yy']
z = vol[vi]['zz']
az = vol[vi]['az'][0,:].round(0)
r = (x[0,0,:]**2 + y[0,0,:]**2)**0.5

az_lims = [270, 0]
azi1 = np.where(az >= az_lims[0])[0][0]
azi2 = np.where(az >= az_lims[1])[0][0]
ri1 = np.where(r >= 0)[0][0]
ri2 = np.where(r >= 6)[0][1]
ri = slice(ri1,ri2)

if az_lims[1] < az_lims[0]:
    xx = np.append(x[:,azi1:,ri], x[:,:azi2,ri], axis=1)
    yy = np.append(y[:,azi1:,ri], y[:,:azi2,ri], axis=1)
    zz = np.append(z[:,azi1:,ri], z[:,:azi2,ri], axis=1)
    
    dbz = np.append(vol[vi]['dbz'][1,azi1:,ri], vol[vi]['dbz'][1,:azi2,ri], axis=0)
    vel = np.append(vol[vi]['vel'][1,azi1:,ri], vol[vi]['vel'][1,:azi2,ri], axis=0)
    zvort = np.append(vol[vi]['zvort'][:,azi1:,ri], vol[vi]['zvort'][:,:azi2,ri], axis=1)
    hvort = np.append(vol[vi]['hvort'][:,azi1:,ri], vol[vi]['hvort'][:,:azi2,ri], axis=1)
else:
    xx = x[:,slice(azi1,azi2),ri]
    yy = y[:,slice(azi1,azi2),ri]
    zz = z[:,slice(azi1,azi2),ri]
    dbz = vol[vi]['dbz'][1,slice(azi1,azi2),ri]
    vel = vol[vi]['vel'][1,slice(azi1,azi2),ri]
    zvort = vol[vi]['zvort'][:,slice(azi1,azi2),ri]
    hvort = vol[vi]['hvort'][:,slice(azi1,azi2),ri]

zvort = gaussian_filter(zvort, 2)
hvort = gaussian_filter(hvort, 2)


#%%

from scipy.interpolate import griddata



xx = vol[vi]['xx'][:,:,:230]
yy = vol[vi]['yy'][:,:,:230]
zz = vol[vi]['zz'][:,:,:230]
dbz = vol[vi]['dbz'][1,:,:230]
vel = vol[vi]['vel'][1,:,:230]
zvort = vol[vi]['zvort'][:,:,:230]
hvort = vol[vi]['hvort'][:,:,:230]

pts = np.array( (zz.flatten(), yy.flatten(), xx.flatten()) ).T
Z,Y,X = np.meshgrid(np.arange(0,1.5,0.1), np.arange(0,4.503,0.03), np.arange(-4.5,0.03,0.03), indexing='ij')
pts0 = np.array( (yy[1,:,:].flatten(), xx[1,:,:].flatten()) ).T

dbz_vals = dbz.flatten()
vel_vals = vel.flatten()
zv_vals = zvort.flatten()
hv_vals = hvort.flatten()

dbz_grid = griddata(pts0, dbz_vals, (Y[1,:,:],X[1,:,:]))
vel_grid = griddata(pts0, vel_vals, (Y[1,:,:],X[1,:,:]))
zvort_grid = griddata(pts, zv_vals, (Z,Y,X))
hvort_grid = griddata(pts, hv_vals, (Z,Y,X))

# dbz = dbz_grid[1,:,:]
# vel = vel_grid[1,:,:]
zvort_grid = np.ma.masked_array(zvort_grid, np.isnan(zvort_grid))
hvort_grid = np.ma.masked_array(hvort_grid, np.isnan(hvort_grid))
#%%
zvort = gaussian_filter(zvort_grid[1:,:,:], 1.5)
hvort = gaussian_filter(hvort_grid[1:,:,:], 1.5)

#%%

fig,(ax1,ax2) = plt.subplots(1,2,figsize=(10,4), sharex=True, sharey=True, subplot_kw=dict(box_aspect=1), layout='constrained')

plot_cfill(X[1,:,:], Y[1,:,:], np.nanmax(zvort[2:,:,:],axis=0), 'vort', ax1, datalims=[0,0.001], xlims=[-6,0], ylims=[0,6])
ax1.set_title(f"{filetime} UTC Max vertical pseudovorticity", fontsize=14)
ax1.set_xlabel('E-W distance from radar (km)', fontsize=12)
ax1.set_ylabel('N-S distance from radar (km)', fontsize=12)
ax1.plot(vol[vi]['xx'][eli,azi,:], vol[vi]['yy'][eli,azi,:], '--k', linewidth=1.25)

plot_cfill(X[1,:,:], Y[1,:,:], np.nanmax(np.abs(hvort[3:,:,:]),axis=0), 'vort', ax2, datalims=[0,0.001], xlims=[-6,0], ylims=[0,6])
ax2.set_title(f"{filetime} UTC Max horizontal pseudovorticity", fontsize=14)
ax2.set_xlabel('E-W distance from radar (km)', fontsize=12)
ax2.set_ylabel('N-S distance from radar (km)', fontsize=12)


plt.show()



fig,(ax1,ax2) = plt.subplots(1,2,figsize=(10,4), sharex=True, sharey=True, subplot_kw=dict(box_aspect=1), layout='constrained')

plot_cfill(Y[1:,:,83], Z[1:,:,83], zvort[:,:,83], 'vort', ax1, datalims=[0,0.001], xlims=[0,4.5], ylims=[0,1.5])
ax1.set_title(f"{filetime} UTC Max vertical pseudovorticity", fontsize=14)
ax1.set_xlabel('E-W distance from radar (km)', fontsize=12)
ax1.set_ylabel('N-S distance from radar (km)', fontsize=12)
# ax1.plot(vol[vi]['xx'][eli,azi,:], vol[vi]['yy'][eli,azi,:], '--k', linewidth=1.25)

plot_cfill(Y[1:,:,83], Z[1:,:,83], hvort[:,:,83], 'vort', ax2, datalims=[0,0.001], xlims=[0,4.5], ylims=[0,1.5])
ax2.set_title(f"{filetime} UTC Max horizontal pseudovorticity", fontsize=14)
ax2.set_xlabel('E-W distance from radar (km)', fontsize=12)
ax2.set_ylabel('N-S distance from radar (km)', fontsize=12)

plt.show()

#%%


fig1 = go.Figure(data=[
    go.Isosurface(
    x=X[:-2, :-17, 17:].flatten(),
    y=Y[:-2, :-17, 17:].flatten(),
    z=Z[:-2, :-17, 17:].flatten(),
    value=zvort[:-2, :-17, 17:].flatten(),
    colorscale='reds',
    showscale=False,
    opacity=0.9,
    isomin=0.0004,
    isomax=0.0006,
    surface_count=3,
    caps=dict(x_show=False, y_show=False, z_show=False),
    slices_z=dict(show=True, locations=[0])
    ),
    ])

fig1.update_layout(
    margin=dict(t=0, l=0, b=0),
    scene_camera_eye=dict(x=1.5, y=0, z=0.05)) # 2, 1, 0.1
fig1.show(renderer='browser')



fig2 = go.Figure(data=[
    go.Isosurface(
    x=X[:-2, :-17, 17:].flatten(),
    y=Y[:-2, :-17, 17:].flatten(),
    z=Z[:-2, :-17, 17:].flatten(),
    value=-1*hvort[:-2, :-17, 17:].flatten(),
    colorscale='blues',
    showscale=False,
    opacity=0.9,
    isomin=0.0004,
    isomax=0.0006,
    surface_count=3,
    caps=dict(x_show=False, y_show=False, z_show=False),
    slices_z=dict(show=True, locations=[0]),
    ),
    
    go.Isosurface(
    x=X[:-2, :-17, 17:].flatten(),
    y=Y[:-2, :-17, 17:].flatten(),
    z=Z[:-2, :-17, 17:].flatten(),
    value=hvort[:-2, :-17, 17:].flatten(),
    colorscale='Reds',
    showscale=False,
    opacity=0.9,
    isomin=0.0005,
    isomax=0.0008,
    surface_count=2,
    caps=dict(x_show=False, y_show=False, z_show=False),
    ),
    
    # go.Surface(
    # x=X[1, :-17, 17:],
    # y=Y[1, :-17, 17:],
    # z=0.001*dbz_grid,
    # colorscale=[[0.0,'green'], [0.33,'gold'], [0.67,'orange'], [1.0,'red']],
    # showscale=False,
    # opacity=0.9,
    # ),
    ])

fig2.update_layout(
    margin=dict(t=0, l=0, b=0),
    scene_camera_eye=dict(x=1.5, y=0, z=0.05)) # 2, 0.6, 1
fig2.show(renderer='browser')































