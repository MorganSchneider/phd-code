# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 13:22:01 2023

@author: morgan.schneider
"""

import numpy as np
import matplotlib.pyplot as plt

#%%

min_el = 1
max_el = 15
ntilts = 14


els15a = np.linspace(1, 15, 12)
els15b = np.linspace(1, 15, 13)
els15c = np.linspace(1, 15, 14)
els30a = np.linspace(1, 30, 12)
els30b = np.linspace(1, 30, 13)
els30c = np.linspace(1, 30, 14)
#els_sbu1 = np.linspace(1,15,30)
#els_sbu2 = np.linspace(1,30,30)
#els_srs = np.array([0.8, 1.8, 3.0, 4.2, 5.4, 6.7, 7.8, 9.0, 10.4])
#els_dow = np.array([0.5, 1.2, 1.9, 2.8, 3.8, 5.0, 6.5, 8.0, 9.5])

rr = np.linspace(0,20,41)
els = els15c
z5 = 5 * np.sin(els/180*np.pi)
z10 = 10 * np.sin(els/180*np.pi)
z15 = 15 * np.sin(els/180*np.pi)
z20 = 20 * np.sin(els/180*np.pi)

zz = [rr * np.sin(els[i]/180*np.pi) for i in range(len(els))]


fig = plt.figure(figsize=(14,8))
plt.axes().set_aspect(1.0)
for i in range(len(els)):
    plt.plot(rr, zz[i], '-k', linewidth=1)
    plt.text(rr[-5], zz[i][-5], f"{els[i]:.1f}$^o$")
    plt.text(rr[-1], zz[i][-1], f"{zz[i][-1]:.2f}")
    plt.text(rr[20], zz[i][20], f"{zz[i][20]:.2f}")
plt.xlabel('range (km)')
plt.ylabel('height (km)')
plt.title(f"Radar beam heights ({max_el} deg, {ntilts} tilts)")
#plt.title(f"Radar beam heights (SKYLER shallow)")
plt.grid(b=True)
plt.xlim(0,20)
plt.ylim(0,np.ceil(zz[-1][-1]))
plt.show()

#%%






