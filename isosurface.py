#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 12:27:13 2023

@author: morgan.schneider
"""

import plotly.graph_objects as go
from CM1utils import *
from scipy.ndimage import gaussian_filter


#%%

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

#%%

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
    ])

fig.update_layout(
    margin=dict(t=0, l=0, b=0),
    scene_camera_eye=dict(x=-0.3, y=-2.2, z=0.2))
fig.show(renderer='browser')

#%%

fig = go.Figure(data=[
    go.Isosurface(
    x=X.flatten(),
    y=Y.flatten(),
    z=Z.flatten(),
    value=prsp.flatten(),
    colorscale='Blues_r',
    showscale=False,
    opacity=0.9,
    isomin=-10,
    isomax=-10,
    surface_count=1,
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
    ])

fig.update_layout(
    margin=dict(t=0, l=0, b=0),
    scene_camera_eye=dict(x=-0.3, y=-2.2, z=0.2))
fig.show(renderer='browser')
