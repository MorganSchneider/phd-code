# -*- coding: utf-8 -*-
"""
Created on Sat Sep 24 13:09:10 2022

@author: morgan.schneider
"""

# Load modules #

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from glob import glob
from metpy.calc import wind_components, bunkers_storm_motion
from metpy.calc import cape_cin, most_unstable_cape_cin, parcel_profile
from metpy.units import units
from metpy.plots import SkewT, Hodograph
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


#%%
save_flag = False
plt.close('all')
# Lv = 2.5e6; #J/kg
# Rv = 461.5; #J/kg/K
# e0 = 6.11; #hPa
# eps = 0.622;
# T0 = 273.; #K
k = 0.286;

fp = '/Users/morgan.schneider/Documents/PERiLS_LIDAR_Soundings/IOP4/'
#fp = '/Users/morgan.schneider/Documents/TORUS_Far_Field_Soundings/'

files = glob(fp+'*.csv')
#files = [files[3]] #t = 1800
for file in files:
    df = pd.read_csv(file, header=2)
    #heads = list(df.columns)
    #tmp = [heads[i] for i in range(33,38)]
    #df = df.drop(columns=tmp)
    
    T = list(df['Filtered Temperature (K)'])
    Td = list(df['Filtered Dewpoint (K)'])
    pres = list(df['Filtered Pressure (mb)'])
    wd = list(df['Filtered Wind Dir'])
    ws = list(df['Filtered Wind Spd (m/s)'])
    z = list(df['Filtered Altitude (m)'])
    qv = list(df['Calculated Mixing Ratio'])
    
    wnd = wind_components(ws*units('m/s'), wd*units.deg)
    u = wnd[0].magnitude[:]
    v = wnd[1].magnitude[:]
    
    # e = [e0*np.exp((Lv/Rv)*(1/T0 - 1/Td[i])) for i in range(len(Td))]
    # qv2 = [1000*eps*e[i]/pres[i] for i in range(len(pres))]
    
    theta = [T[i]*(1000/pres[i])**k for i in range(len(T))]
    
    ### header:    sfc pres (mb)   sfc theta (K)   sfc qv (g/kg)
    ### body:   z (m)    theta (K)    qv (g/kg)     u (m/s)     v(m/s)
    
    if save_flag:
        #fsave = file[60:-4]+'.txt' # for TORUS soundings
        fsave = file[62:-4]+'.txt'
        hd = np.zeros((1,3))
        hd[0][0] = pres[0]
        hd[0][1] = theta[0]
        hd[0][2] = qv[0]
        np.savetxt(fp+fsave, hd, fmt='%f')
        
        dat = {'z': z, 'theta': theta, 'qv': qv, 'u': u, 'v': v}
        dfs = pd.DataFrame(data=dat, dtype=float)
        
        with open(fp+fsave, 'a') as ff:
            ff.write(dfs.to_string(header=False, index=False))
    
    if False:
        print(f"---{file[-10:-4]}---")
        print(f"Profile top: {max(z):.0f} m")
        if max(z) >= 6123:
            bwnd = bunkers_storm_motion(pres*units.hPa, u*units('m/s'), v*units('m/s'), z*units.m)
            u_bunkersR = bwnd[0].magnitude[0]
            v_bunkersR = bwnd[0].magnitude[1]
            u_bunkersL = bwnd[1].magnitude[0]
            v_bunkersL = bwnd[1].magnitude[1]
            u_sfc6km = bwnd[2].magnitude[0]
            v_sfc6km = bwnd[2].magnitude[1]
            print(f"0-6 shear vector: {u_sfc6km:.1f}, {v_sfc6km:.1f}")
            print(f"Bunkers RM vector: {u_bunkersR:.1f}, {v_bunkersR:.1f}")
        else:
            print("0-6 shear vector: --")
            print("Bunkers RM vector: --")
        
        T_parcel = parcel_profile(pres*units.hPa, T[0]*units.K, Td[0]*units.K)
        thermo = cape_cin(pres*units.hPa, T*units.K, Td*units.K, T_parcel)
        thermo_mu = most_unstable_cape_cin(pres*units.hPa, T*units.K, Td*units.K)
        cape = thermo[0].magnitude
        cin = thermo[1].magnitude
        mucape = thermo_mu[0].magnitude
        mucin = thermo_mu[1].magnitude
        print(f"CAPE/CIN: {cape:.0f}/{cin:.0f}")
        print(f"MUCAPE/MUCIN: {mucape:.0f}/{mucin:.0f}")
    

    plt.rcParams['figure.figsize']=(9,9)
    
    T_parcel = parcel_profile(pres*units.hPa, T[0]*units.K, Td[0]*units.K)

    skew = SkewT()
    skew.plot(np.array(pres), np.array(Td)-273.15, '-g', linewidth=2)
    skew.plot(np.array(pres), np.array(T)-273.15, '-r', linewidth=2)
    skew.plot(np.array(pres), np.array(T_parcel.magnitude[:])-273.15, '-k', linewidth=2)
    # Add the relevant special lines
    skew.plot_dry_adiabats()
    skew.plot_moist_adiabats()
    skew.plot_mixing_lines()
    skew.ax.set_ylim(1000, 100)
    skew.ax.set_xlim(-40, 30)
    plt.title(f"{file[-10:-4]} UTC input sounding")
    # Create a hodograph
    ax_hod = inset_axes(skew.ax, '40%', '40%', loc=1)
    h = Hodograph(ax_hod, component_range=40.)
    h.add_grid(increment=20)
    h.plot(np.array(u), np.array(v), color='k', linewidth=1.5)

    plt.show()
    #plt.savefig(fp+'sounding'+file[-20:-4]+'.png', dpi=300)



