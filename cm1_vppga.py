#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 14:26:25 2024

@author: morgan.schneider

Calculate layer VPPGA
"""

####################
### Load modules ###
####################

from CM1utils import *

#%%

fp = '/Volumes/Promise_Pegasus_70TB/merger/merger-125m/'

ds = nc.Dataset(fp+f"base/cm1out_000001.nc")
z = ds.variables['z'][:].data
iz1 = np.where(z >= 1)[0][0]
dz = z[iz1] - z[0]
rho0 = np.mean(ds.variables['rho'][:].data[0,0:iz1,-1,-1])
ds.close()

vppga_dl = {}
vppga_dn = {}
vppga_b = {}
for fn in np.arange(14,74):
    print(f"dyn_{fn:06d}")
    ds = nc.Dataset(fp+f"pp/dyn_{fn:06d}.nc")
    time = ds.variables['time'][:].data[0]/60
    dp_dl = ds.variables['p_dl'][:].data[iz1,:,:] - ds.variables['p_dl'][:].data[0,:,:]
    dp_dn = ds.variables['p_dn'][:].data[iz1,:,:] - ds.variables['p_dn'][:].data[0,:,:]
    dp_b = ds.variables['p_b'][:].data[iz1,:,:] - ds.variables['p_b'][:].data[0,:,:]
    ds.close()
    
    pgf_dl = -1/rho0 * dp_dl / dz
    pgf_dn = -1/rho0 * dp_dn / dz
    pgf_b = -1/rho0 * dp_b / dz
    
    vppga_dl.update({f"f{fn:06d}":pgf_dl})
    vppga_dn.update({f"f{fn:06d}":pgf_dn})
    vppga_b.update({f"f{fn:06d}":pgf_b})
    
vppga = {'dl':vppga_dl, 'dn':vppga_dn, 'b':vppga_b}

dbfile = open(fp+'pp/vppga_1km.pkl', 'wb')
pickle.dump(vppga, dbfile)
dbfile.close()
    
    




