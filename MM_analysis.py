#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  4 14:52:22 2025

@author: morgan.schneider
"""

from RaxpolUtils import *
from MM_functions import *

#%% Load MM data and vortex locations

# Probe 1 - Tyler's bias-corrected files
ds = nc.Dataset('/Users/morgan.schneider/Documents/perils2023/iop2/mesonet/Probe_1_IOP2_QC_all.nc')
P1 = dict(time=datetime.fromtimestamp(ds.variables['time'][:].data+21600),
          lat=ds.variables['lat'][:].data,
          lon=ds.variables['lon'][:].data,
          elev=ds.variables['elev'][:].data,
          temp_raw=ds.variables['temp'][:].data,
          pres_raw=ds.variables['pres'][:].data,
          rh_raw=ds.variables['rh'][:].data,
          wspd_raw=ds.variables['wspd'][:].data,
          wdir_raw=ds.variables['wdir'][:].data,
          wspd_corr=ds.variables['wspd_corr'][:].data,
          wdir_corr=ds.variables['wdir_corr'][:].data,
          u_corr=ds.variables['u_corr'][:].data,
          v_corr=ds.variables['v_corr'][:].data,
          temp_corr=ds.variables['temp_unbiased'][:].data,
          rh_corr=ds.variables['rh_unbiased'][:].data,
          pres_corr=ds.variables['pres_unbiased'][:].data)
ds.close()

# Probe 2 - Tyler's bias-corrected files
ds = nc.Dataset('/Users/morgan.schneider/Documents/perils2023/iop2/mesonet/Probe_2_IOP2_QC_all.nc')
P2 = dict(time=datetime.fromtimestamp(ds.variables['time'][:].data+21600),
          lat=ds.variables['lat'][:].data,
          lon=ds.variables['lon'][:].data,
          elev=ds.variables['elev'][:].data,
          temp_raw=ds.variables['temp'][:].data,
          pres_raw=ds.variables['pres'][:].data,
          rh_raw=ds.variables['rh'][:].data,
          wspd_raw=ds.variables['wspd'][:].data,
          wdir_raw=ds.variables['wdir'][:].data,
          wspd_corr=ds.variables['wspd_corr'][:].data,
          wdir_corr=ds.variables['wdir_corr'][:].data,
          u_corr=ds.variables['u_corr'][:].data,
          v_corr=ds.variables['v_corr'][:].data,
          temp_corr=ds.variables['temp_unbiased'][:].data,
          rh_corr=ds.variables['rh_unbiased'][:].data,
          pres_corr=ds.variables['pres_unbiased'][:].data)
ds.close()




fn = glob(f"/Users/morgan.schneider/Documents/perils2023/iop2/raxpol/CFradial_cal/*_081928_*.nc")[0]
rax = pyart.io.read(fn)
rax_lat = rax.latitude['data'][0]
rax_lon = rax.longitude['data'][0]

dbfile = open('/Users/morgan.schneider/Documents/perils2023/iop2/raxpol_vortex_locs.pkl', 'rb')
vortex_locs = pickle.load(dbfile)
dbfile.close()



# convert all the vortex loc x/y to lat/lons, then figure out how to find and track a center point
# --> center on the visible couplet (vortex 3)? use the scans where it's visible to get couplet motion
# and then back out the position for the earlier times

meso_times = [datetime(2023,3,3,8,20,0), datetime(2023,3,3,8,20,30), datetime(2023,3,3,8,21,0),
              datetime(2023,3,3,8,21,30), datetime(2023,3,3,8,22,30), datetime(2023,3,3,8,23,0),
              datetime(2023,3,3,8,23,30), datetime(2023,3,3,8,24,0), datetime(2023,3,3,8,24,30),
              datetime(2023,3,3,8,25,0)]












