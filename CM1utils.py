
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 11:04:04 2022

@author: morgan.schneider

Functions for reading and plotting CM1 outputs.
"""

####################
### Load modules ###
####################

import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc
import cmocean
import pyart #need an earlier version of xarray -> 0.20.2 or earlier
import pickle
# import metpy.calc as mc
# from metpy.plots import SkewT, Hodograph
# from metpy.units import units
# from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# from mpl_toolkits.axes_grid1 import make_axes_locatable
# from glob import glob
import wrf
from scipy.ndimage import gaussian_filter

#%%
######################
### Define classes ###
######################

# Mimics the Matlab struct array
# Literally just a dict that uses dot indexing instead of bracket indexing
# so not actually useful but it took me far too long to discover dicts
class struct():
    def __init__(self,**kwargs):
        self.Set(**kwargs)
    def Set(self,**kwargs):
        self.__dict__.update(**kwargs)
    def SetAttr(self,var,val):
        self.__dict__[var] = val
    def keys(self):
        print(self.__dict__.keys())


#%%
########################
### Define functions ###
########################


# Define colormaps for common cm1 variables
cmaps = {
    'u':       {'cm': cmocean.cm.balance, 'label': "u (m s$^{-1}$)"},
    'v':       {'cm': cmocean.cm.balance, 'label': "v (m s$^{-1}$)"},
    'w':       {'cm': cmocean.cm.balance, 'label': "w (m s$^{-1}$)"},
    'wspd':    {'cm': 'pyart_HomeyerRainbow', 'label': "Wind speed (m s$^{-1}$)"},
    'th':      {'cm': 'pyart_HomeyerRainbow', 'label': "\u03B8 (K)"},
    'thpert':  {'cm': cmocean.cm.curl, 'label': "\u03B8' (K)"},
    'thr':     {'cm': 'pyart_HomeyerRainbow', 'label': "\u03B8\u1D68 (K)"},
    'thrpert': {'cm': cmocean.cm.curl, 'label': "\u03B8'\u1D68 (K)"},
    'qv':      {'cm': 'viridis_r', 'label': "$q_v$ (g kg$^{-1}$)"},
    'qvpert':  {'cm': 'BrBG', 'label': "$q_v$' (g kg$^{-1}$)"},
    'qx':      {'cm': 'viridis_r', 'label': "$q$ (g kg$^{-1}$)"},
    'prs':     {'cm': 'pyart_HomeyerRainbow', 'label': "p (hPa)"},
    'prspert': {'cm': cmocean.cm.balance, 'label': "p' (hPa)"},
    'pi':      {'cm': 'pyart_HomeyerRainbow', 'label': "\u03C0 (nondimensional)"},
    'pipert':  {'cm': cmocean.cm.balance, 'label': "\u03C0' (nondimensional)"},
    'rho':     {'cm': 'pyart_HomeyerRainbow', 'label': "\u03C1 (kg m$^{-3}$)"},
    'xvort':   {'cm': cmocean.cm.balance, 'label': "\u03BE (s$^{-1}$)"},
    'yvort':   {'cm': cmocean.cm.balance, 'label': "\u03B7 (s$^{-1}$)"},
    'zvort':   {'cm': cmocean.cm.balance, 'label': "\u03B6 (s$^{-1}$)"},
    'OW':      {'cm': cmocean.cm.balance, 'label': "OW (s$^{-2}$)"},
    'divh':    {'cm': cmocean.cm.balance, 'label': "\u25BD$_H$u (s$^{-1}$)"},
    'dbz':     {'cm': 'pyart_NWSRef', 'label': "$Z_H$ (dBZ)"},
    'pgfz':    {'cm': cmocean.cm.balance, 'label': "VPPGA (m s$^{-2}$)"},
    'cape':    {'cm': 'pyart_HomeyerRainbow', 'label': "CAPE (J kg$^{-1}$)"},
    'srh':     {'cm': 'pyart_HomeyerRainbow', 'label': "SRH (m$^{2}$ s$^{-2}$)"},
    'z':       {'cm': 'pyart_HomeyerRainbow', 'label': "Height (m)"},
    'scp':     {'cm': 'pyart_HomeyerRainbow', 'label': "SCP"},
    'stp':     {'cm': 'pyart_HomeyerRainbow', 'label': "STP"},
    'uh':      {'cm': 'pyart_HomeyerRainbow', 'label': "UH (m$^{2}$ s$^{-2}$)"}
}


# # Read CM1 output into struct object
# def read_cm1out(fname, dsvars=None):
#     # fname : full path to data file
#     # dsvars: list of the names of desired variables to load
    
#     # Open output file
#     ds = nc.Dataset(fname)
#     if dsvars is not None:
#         dsvars = np.array(dsvars)
#     else:
#         dsvars = np.array(list(ds.variables.keys()))
#     # Read data into a struct object (defined above)
#     df = struct()
#     #df = {}
#     for i in range(len(dsvars)):
#         df.SetAttr(dsvars[i], ds.variables[dsvars[i]][:].data)
#         #df.update({dsvars[i]: ds.variables[dsvars[i]][:].data})
#     ds.close()
    
#     return df


# Wrapper function for pcolormesh
def plot_cfill(x, y, data, field, ax, datalims=None, xlims=None, ylims=None,
               cmap=None, cbar=True, **kwargs):
    if cmap is None:
        cm, cb_label = cmaps[field]['cm'], cmaps[field]['label']
    else:
        cm, cb_label = cmap, cmaps[field]['label']
    
    if datalims is None:
        datamin = None
        datamax = None
    else:
        datamin = datalims[0]
        datamax = datalims[1]
    
    # Create the plot
    c = ax.pcolormesh(x, y, data, vmin=datamin, vmax=datamax, cmap=cm, **kwargs)

    # Format the colorbar
    # c.cmap.set_bad('grey', 1.0)
    if cbar:
        cb = plt.colorbar(c, ax=ax, extend='min')
        cb.set_label(cb_label)
    
    if xlims is not None:
        ax.set_xlim(xlims[0], xlims[1])
    if ylims is not None:
        ax.set_ylim(ylims[0], ylims[1])
    
    return c


# Get vertical cross-sections of 3D fields along any diagonal line, plus a 1D vector for the diagonal horizontal coordinate
def vert_cross_section(field, x, y, start=[0,0], end=[-1,-1]):
    # start: xy coordinates of cross-section start point (array or tuple)
    # end: xy coordinates of cross-section end point (array or tuple)
    ix1 = np.where(x >= start[0])[0][0]
    iy1 = np.where(y >= start[1])[0][0]
    ix2 = np.where(x >= end[0])[0][0]
    iy2 = np.where(y >= end[1])[0][0]
    
    xy = wrf.xy(field, start_point=(ix1,iy1), end_point=(ix2,iy2))
    x_cross = wrf.interp2dxy(np.tile(x, (field.shape[0],field.shape[1],1)), xy)
    field_cross = wrf.interp2dxy(field, xy)
    
    return field_cross.data, x_cross[0].data
    

# Get magnitude of wind vectors projected onto an angle alpha (meant for getting horizontal winds along diagonal cross sections)
def proj_winds(u, v, proj_angle):
    if proj_angle > 2*np.pi:
        proj_angle = proj_angle*np.pi/180 # convert to rads
    
    U_proj = u*np.sin(proj_angle) + v*np.cos(proj_angle)
    V_proj = u*np.cos(proj_angle) - v*np.sin(proj_angle)
    nu = U_proj * np.sin(proj_angle)
    nv = U_proj * np.cos(proj_angle)
    
    return U_proj,nu,nv


def save_to_pickle(data, pkl_fname, new_pkl=False):
    if new_pkl is True:
        dbfile = open(pkl_fname, 'wb')
        pickle.dump(data, dbfile)
        dbfile.close()
    else:
        dbfile = open(pkl_fname, 'rb')
        save_data = pickle.load(dbfile)
        dbfile.close()
        
        save_data.update(data)
        dbfile = open(pkl_fname, 'wb')
        pickle.dump(save_data, dbfile)
        dbfile.close()
    


# Calculate individual CM1 reflectivity contributions (NSSL 2-moment microphysics only)
def calc_dbz_contributions(qx, cx, rho, field, plot=False, **kwargs):
    # qx = hydrometeor mixing ratio in kg/kg
    # cx = hydrometeor number concentration in #/kg
    # rho = air density matrix in kg/m^3
    # field = 'rain', 'ice', 'snow', 'graupel', or 'hail'
    
    rho_l = 1000 # liquid water density
    
    if field == 'rain':
        z = np.where(cx>0, 20*1e18*(6/(np.pi*rho_l))**2 * (rho*qx)**2 / cx, 0)
        levs = [0.1, 1]
    elif field == 'ice':
        z = 750*(np.where(1e3*qx*rho > 1, 1, 1e3*qx*rho))**1.98
        levs = [0.01, 0.02]
    elif field == 'snow':
        z = np.where(cx>0, 1e18*323.3226*(0.106214**2)*0.887762*0.189*(rho*qx)**2 / (4.590844*(917**2)*(0.2**(4/3))*cx), 0)
        levs = [0.25, 1]
    elif field == 'graupel':
        z = np.where(cx>0, 0.224*20*1e18*(6/(np.pi*rho_l))**2 * (rho*qx)**2 / cx, 0)
        levs = [0.1, 0.2]
    elif field == 'hail':
        a = 0.5
        g1 = (6+a)*(5+a)*(4+a) / ((3+a)*(2+a)*(1+a))
        z = np.where(cx>0, 0.224*g1*1e18*(6/(np.pi*rho_l))**2 * (rho*qx)**2 / cx, 0)
        levs = [0.001, 0.003]
    
    if plot:
        dbz = np.where(z>0, 10*np.log10(z), 0)
        xf = kwargs.get('xf')
        xh = kwargs.get('xh')
        zf = kwargs.get('zf')
        zh = kwargs.get('zh')
        
        fig, (ax1,ax2) = plt.subplots(1,2,figsize=(14,6))
        
        c1 = plot_cfill(xf, zf, dbz[0,:,1000,:], 'dbz', ax1, datalims=[0,80])
        ax1.contour(xh, zh, 1000*qx[0,:,1000,:], levels=levs, colors='k', linewidths=1)
        ax1.set_xlabel('x (km)')
        ax1.set_ylabel('z (km)')
        ax1.set_title(f"Z_{field}", fontsize=14)
        ax1.set_xlim(-120,-20)
        ax1.set_ylim(0,10)

        c2 = plot_cfill(xf, zf, 1000*qx[0,:,1000,:], 'qr', ax2)
        ax2.contour(xh, zh, 1000*qx[0,:,1000,:], levels=levs, colors='k', linewidths=1)
        ax2.set_xlabel('x (km)')
        ax2.set_ylabel('z (km)')
        ax2.set_title(f"q_{field}", fontsize=14)
        ax2.set_xlim(-120,-20)
        ax2.set_ylim(0,10)
        plt.show()
        
    return z
    

# Recalculate brightband-corrected CM1 reflectivity (NSSL 2-moment microphysics only)
def calc_dbz(ds):
    # Load data
    rho = ds.variables['rho'][:].data # air density (kg/m3)
    rho_l = 1000 # liquid water density (kg/m3)
    
    # print("Reading mixing ratios...")
    qr = ds.variables['qr'][:].data # rain mixing ratio (kg/kg)
    qi = ds.variables['qi'][:].data # ice crystal mixing ratio
    qs = ds.variables['qs'][:].data # snow mixing ratio
    qg = ds.variables['qg'][:].data # graupel mixing ratio
    qhl = ds.variables['qhl'][:].data # hail mixing ratio
    
    # print("Reading number concentrations...")
    crw = ds.variables['crw'][:].data # rain number concentration (#/kg)
    #cci = ds.variables['cci'][:].data # ice crystal number concentration
    csw = ds.variables['csw'][:].data # snow number concentration
    chw = ds.variables['chw'][:].data # graupel number concentration
    chl = ds.variables['chl'][:].data # hail number concentration
    
    # print("Calculating reflectivity contributions...")
    with np.errstate(divide='ignore', invalid='ignore'):
        z_rain = np.where(crw>0, 20*1e18*(6/(np.pi*rho_l))**2 * (rho*qr)**2 / crw, 0)
        del crw,qr
        z_ice = np.where(qi>1e-4, 750*(np.where(1e3*qi*rho>1, 1, 1e3*qi*rho))**1.98, 0)
        del qi
        z_snow = np.where(csw>0, 1e18*323.3226*(0.106214**2)*0.887762*0.189*(rho*qs)**2 / (4.590844*(917**2)*(0.2**(4/3))*csw), 0)
        del csw,qs
        z_graupel = np.where(chw>0, 0.224*20*1e18*(6/(np.pi*rho_l))**2 * (rho*qg)**2 / chw, 0)
        del chw,qg
        a = 0.5
        g1 = (6+a)*(5+a)*(4+a) / ((3+a)*(2+a)*(1+a))
        z_hail = np.where(chl>0, 0.224*g1*1e18*(6/(np.pi*rho_l))**2 * (rho*qhl)**2 / chl, 0)
        del chl,qhl
        
        z = z_rain + z_ice + z_snow + z_graupel + z_hail
        dbz = np.where(z>0, 10*np.log10(z), 0)
        
        del z_rain,z_ice,z_snow,z_graupel,z_hail,z,rho
        
    return dbz


    
    






