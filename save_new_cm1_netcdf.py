#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  1 14:36:43 2025

@author: morgan.schneider
"""

import netCDF4 as nc
import numpy as np


#%%

def del_var_from_netcdf(in_fp, out_fp, exclude_vars):
    with nc.Dataset(in_fp,'r') as src, nc.Dataset(out_fp,'w') as ds:
        ds.setncatts(src.__dict__)
        for name,dimension in src.dimensions.items():
            ds.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))
        
        for name,variable in src.variables.items():
            if name not in exclude_vars:
                x = ds.createVariable(name, variable.datatype, variable.dimensions)
                ds[name].setncatts(src[name].__dict__)
                ds[name][:] = src[name][:]


fp = '/Volumes/Promise_Pegasus_70TB/merger/merger-125m/'

ds = nc.Dataset(fp+f"base/cm1out_000001.nc")
z = ds.variables['z'][:].data
th0 = ds.variables['th0'][:].data[0,:,0,0]
prs0 = ds.variables['prs0'][:].data[0,:,0,0]
pi0 = ds.variables['pi0'][:].data[0,:,0,0]
qv0 = ds.variables['qv0'][:].data[0,:,0,0]
u0 = ds.variables['u0'][:].data[0,:,0,0]
v0 = ds.variables['v0'][:].data[0,:,0,0]
thr0 = th0 * (1 + 0.61*qv0)
ds.close()


ds = nc.Dataset(fp+f"cm1out_basestate.nc", 'w')
ds.createDimension('nk', len(z))

th = ds.createVariable('th', 'f4', ('nk'))
th.units = 'K'
th.long_name = 'base state potential temperature'
th[:] = th0[:]

prs = ds.createVariable('prs', 'f4', ('nk'))
prs.units = 'Pa'
prs.long_name = 'base state pressure'
prs[:] = prs0[:]

pi = ds.createVariable('pi', 'f4', ('nk'))
pi.units = 'nondimensional'
pi.long_name = 'base state nondimensional pressure'
pi[:] = pi0[:]

thr = ds.createVariable('thr', 'f4', ('nk'))
thr.units = 'K'
thr.long_name = 'base state density potential temperature'
thr[:] = thr0[:]

qv = ds.createVariable('qv', 'f4', ('nk'))
qv.units = 'kg/kg'
qv.long_name = 'base state water vapor mixing ratio'
qv[:] = qv0[:]

u = ds.createVariable('u', 'f4', ('nk'))
u.units = 'm/s'
u.long_name = 'base state zonal wind'
u[:] = u0[:]

v = ds.createVariable('v', 'f4', ('nk'))
v.units = 'm/s'
v.long_name = 'base state meridional wind'
v[:] = v0[:]

ds.close()


for fn in np.arange(13,15):
    if fn <= 13:
        exclude_vars = ['f_cor', 'ztop', 'sws', 'svs', 'sps', 'srs', 'sgs', 'sus', 'shs',
                        'sws2', 'svs2', 'sps2', 'srs2', 'sgs2', 'sus2', 'shs2', 'cref',
                        'pipert', 'ccn', 'ccw', 'crw', 'cci', 'csw', 'chw', 'chl', 'vhw',
                        'vhl', 'dbz', 'dbz2', 'pi0', 'th0', 'prs0', 'qv0', 'u0', 'v0']
    else:
        exclude_vars = ['f_cor', 'ztop', 'cref', 'ccn', 'ccw', 'crw', 'cci', 'csw', 'chw',
                        'chl', 'vhw', 'vhl']
    
    del_var_from_netcdf(fp+f"cm1out_{fn:06d}.nc", fp+f"cm1out_{fn:06d}_v2.nc", exclude_vars)
    
    if fn <= 13:
        ds = nc.Dataset(fp+f"cm1out_{fn:06d}.nc")
        dbz2 = ds.variables['dbz2'][:].data
        thrho = ds.variables['th'][:].data * (1 + 0.61*ds.variables['qv'][:].data - 
                    (ds.variables['qc'][:].data + ds.variables['qr'][:].data + 
                     ds.variables['qi'][:].data + ds.variables['qs'][:].data + 
                     ds.variables['qg'][:].data + ds.variables['qhl'][:].data))
        ds.close()
        
        
        ds = nc.Dataset(fp+f"cm1out_{fn:06d}_v2.nc", 'r+')
        
        dbz = ds.createVariable('dbz', 'f4', ('time', 'nk', 'nj', 'ni'))
        dbz.units = 'DBZ'
        dbz.long_name = 'reflectivity'
        dbz[:,:,:,:] = dbz2[:,:,:,:]
        del dbz,dbz2
        
        thr = ds.createVariable('thr', 'f4', ('nk', 'nj', 'ni'))
        thr.units = 'K'
        thr.long_name = 'perturbation density potential temperature'
        thr[:,:,:] = thrho[:,:,:]
        del thr,thrho
        
        ds.close()
        
    else:
        ds = nc.Dataset(fp+"pp/dyn_{fn:06d}.nc")
        thrho = ds.variables['thrpert'][:].data + np.moveaxis(np.tile(thr0, (2880,2880,1)), -1, 0)
        ds.close()
        
        ds = nc.Dataset(fp+f"cm1out_{fn:06d}_v2.nc", 'r+')
        
        thr = ds.createVariable('thr', 'f4', ('nk', 'nj', 'ni'))
        thr.units = 'K'
        thr.long_name = 'perturbation density potential temperature'
        thr[:,:,:] = thrho[:,:,:]
        del thr,thrho
        
        ds.close()














