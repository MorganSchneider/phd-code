#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:43:24 2023

@author: morgan.schneider

What it says on the package.
"""

####################
### Load modules ###
####################

from CM1utils import *
from metpy.plots import SkewT, Hodograph
import metpy.calc as mc
from metpy.units import units

#%% Define functions for things

# Effective Shear Algorithm
def effective_layer(p, T, Td, h, height_layer=False):
    '''This function determines the effective inflow layer for a convective sounding.
    
    Input:
      - p: sounding pressure with units
      - T: sounding temperature with units
      - Td: sounding dewpoint temperature with units
      - h: sounding heights with units
      
    Returns:
      - pbot/hbot, ptop/htop: pressure/height of the bottom level, pressure/height of the top level
    '''
    
    pbot = None
    
    for i in range(p.shape[0]):
        prof = mc.parcel_profile(p[i:], T[i], Td[i])
        sbcape, sbcin = mc.cape_cin(p[i:], T[i:], Td[i:], prof)
        if sbcape >= 100 * units('J/kg') and sbcin > -250 * units('J/kg'):
            pbot = p[i]
            hbot = h[i]
            bot_idx = i
            break
    if not pbot:
        return None, None
    
    for i in range(bot_idx+1, p.shape[0]):
        prof = mc.parcel_profile(p[i:], T[i], Td[i])
        sbcape, sbcin = mc.cape_cin(p[i:], T[i:], Td[i:], prof)
        if sbcape < 100 * units('J/kg') or sbcin < -250 * units('J/kg'):
            ptop = p[i]
            htop = h[i]
            break
            
    if height_layer:
        return hbot, htop
    else:
        return pbot, ptop

#%% Load data

fp = '/Volumes/Promise_Pegasus_70TB/merger/merger-125m/base/'
fn = 'cm1out_000001.nc'

ds = nc.Dataset(fp+fn)
z = ds.variables['z'][:].data # km
u = ds.variables['uinterp'][:].data[0,:,-1,-1] + 6 # m/s
v = ds.variables['vinterp'][:].data[0,:,-1,-1] # m/s
th = ds.variables['th'][:].data[0,:,-1,-1] # K
qv = ds.variables['qv'][:].data[0,:,-1,-1] # kg/kg
prs = ds.variables['prs'][:].data[0,:,-1,-1] # Pa
ds.close()

# Calculate T, Td
T = th * (prs/100000.)**0.286
e = (qv * prs/100) / (0.622+qv)
Td = 243.5 / ((17.67/(np.log(e/6.112)))-1) + 273.15


#%% Calculate things!

# Effective inflow layer
eff_cm1 = effective_layer(prs*units.Pa, T*units.K, Td*units.K, z*1000*units.m, height_layer=True)
ebot = eff_cm1[0].magnitude
etop = eff_cm1[1].magnitude

# Storm motion estimates
bwnd = mc.bunkers_storm_motion(prs*units.Pa, u*units('m/s'), v*units('m/s'), z*units.km)
u_RM = bwnd[0].magnitude[0] # Bunkers RM
v_RM = bwnd[0].magnitude[1]
u_LM = bwnd[1].magnitude[0] # Bunkers LM
v_LM = bwnd[1].magnitude[1]
u_06 = bwnd[2].magnitude[0] # 0-6 km mean wind
v_06 = bwnd[2].magnitude[1]
wspd_RM = np.sqrt(u_RM**2 + v_RM**2)
wdir_RM = 180 + np.arctan2(v_RM, u_RM)*180/np.pi
wspd_06 = np.sqrt(u_06**2 + v_06**2)
wdir_06 = 180 + np.arctan2(v_06, u_06)*180/np.pi


T_parcel = mc.parcel_profile(prs*units.Pa, T[0]*units.K, Td[0]*units.K)
# Regular CAPE/CIN
CC = mc.cape_cin(prs*units.Pa, T*units.K, Td*units.K, T_parcel)
cape = CC[0].magnitude
cin = CC[1].magnitude
# MUCAPE/MUCIN
MU = mc.most_unstable_cape_cin(prs*units.Pa, T*units.K, Td*units.K, height=z*units.km)
mucape = MU[0].magnitude
mucin = MU[1].magnitude
# SBCAPE/SBCIN
SB = mc.surface_based_cape_cin(prs*units.Pa, T*units.K, Td*units.K)
sbcape = SB[0].magnitude
sbcin = SB[1].magnitude
# MLCAPE/MLCIN
ML = mc.mixed_layer_cape_cin(prs*units.Pa, T*units.K, Td*units.K, height=z*units.km)
mlcape = ML[0].magnitude
mlcin = ML[1].magnitude

# LCL height
lcl = mc.lcl(prs[0]*units.Pa, T[0]*units.K, Td[0]*units.K)
p_lcl = lcl[0].magnitude # Pa
i1 = np.where(prs >= p_lcl)[0][-1] # below LCL
i2 = np.where(prs <= p_lcl)[0][0] # above LCL
z_lcl = z[i2] - (z[i2]-z[i1]) * (p_lcl-prs[i2])/(prs[i1]-prs[i2]) # km
# LFC height
lfc = mc.lfc(prs*units.Pa, T*units.K, Td*units.K)
p_lfc = lfc[0].magnitude # Pa
i1 = np.where(prs >= p_lfc)[0][-1] # below LFC
i2 = np.where(prs <= p_lfc)[0][0] # above LFC
z_lfc = z[i2] - (z[i2]-z[i1]) * (p_lfc-prs[i2])/(prs[i1]-prs[i2]) # km
# EL height
el = mc.el(prs*units.Pa, T*units.K, Td*units.K)
p_el = el[0].magnitude # Pa
i1 = np.where(prs >= p_el)[0][-1] # below EL
i2 = np.where(prs <= p_el)[0][0] # above EL
z_el = z[i2] - (z[i2]-z[i1]) * (p_el-prs[i2])/(prs[i1]-prs[i2]) # km

# 0-1 km bulk shear
shr = mc.bulk_shear(prs*units.Pa, u*units('m/s'), v*units('m/s'), height=z*1000*units.m, bottom=z[0]*1000*units.m, depth=1000*units.m)
u_shr = shr[0].magnitude
v_shr = shr[1].magnitude
shr_01 = np.sqrt(u_shr**2 + v_shr**2)
shr01_dir = 180 + np.arctan2(v_shr, u_shr)*180/np.pi
# 0-3 km bulk shear
shr = mc.bulk_shear(prs*units.Pa, u*units('m/s'), v*units('m/s'), height=z*1000*units.m, bottom=z[0]*1000*units.m, depth=3000*units.m)
u_shr = shr[0].magnitude
v_shr = shr[1].magnitude
shr_03 = np.sqrt(u_shr**2 + v_shr**2)
shr03_dir = 180 + np.arctan2(v_shr, u_shr)*180/np.pi
# 0-6 km bulk shear
shr = mc.bulk_shear(prs*units.Pa, u*units('m/s'), v*units('m/s'), height=z*1000*units.m, bottom=z[0]*1000*units.m, depth=6000*units.m)
u_shr = shr[0].magnitude
v_shr = shr[1].magnitude
shr_06 = np.sqrt(u_shr**2 + v_shr**2)
shr06_dir = 180 + np.arctan2(v_shr, u_shr)*180/np.pi
# Effective layer bulk shear
shr = mc.bulk_shear(prs*units.Pa, u*units('m/s'), v*units('m/s'), height=z*1000*units.m, bottom=ebot*units.m, depth=(z_el*1000-ebot)*units.m)
u_shr = shr[0].magnitude
v_shr = shr[1].magnitude
EBS = np.sqrt(u_shr**2 + v_shr**2)
EBS_dir = 180 + np.arctan2(v_shr, u_shr)*180/np.pi

# Storm-relative helicity (0-500 m, 0-1 km, 0-3 km, effective layer)
srh = mc.storm_relative_helicity(z*1000*units.m, u*units('m/s'), v*units('m/s'), depth=500*units.m, storm_u=u_RM*units('m/s'), storm_v=v_RM*units('m/s'))
SRH_500 = srh[2].magnitude
srh = mc.storm_relative_helicity(z*1000*units.m, u*units('m/s'), v*units('m/s'), depth=1000*units.m, storm_u=u_RM*units('m/s'), storm_v=v_RM*units('m/s'))
SRH_1km = srh[2].magnitude
srh = mc.storm_relative_helicity(z*1000*units.m, u*units('m/s'), v*units('m/s'), depth=3000*units.m, storm_u=u_RM*units('m/s'), storm_v=v_RM*units('m/s'))
SRH_3km = srh[2].magnitude
srh = mc.storm_relative_helicity(z*1000*units.m, u*units('m/s'), v*units('m/s'), depth=(etop-ebot)*units.m, bottom=ebot*units.m, storm_u=u_RM*units('m/s'), storm_v=v_RM*units('m/s'))
ESRH = srh[2].magnitude

# Other convective parameters
crit = mc.critical_angle(prs*units.Pa, u*units('m/s'), v*units('m/s'), z*1000*units.m, u_RM*units('m/s'), v_RM*units('m/s'))
crit_angle = crit.magnitude
scp = mc.supercell_composite(mucape*units('J/kg'), ESRH*units('m^2/s^2'), EBS*units('m/s'))
SCP = scp[0].magnitude
stp = mc.significant_tornado(sbcape*units('J/kg'), z_lcl*units.km, SRH_500*units('m^2/s^2'), shr_06*units('m/s'))
STP = stp[0].magnitude




















