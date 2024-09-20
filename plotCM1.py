# -*- coding: utf-8 -*-
"""
Created on Wed May 31 11:26:54 2023
@author: morgan.schneider

Plotting CM1 outputs
"""

####################
### Load modules ###
####################

from CM1utils import *


#%% One big overview plot

xlims = [[-70,10], [-62,18], [-54,26], [-46,34], [-38,42]]
ylims = [[-140,-60], [-122,-42], [-104,-24], [-86,-6], [-68,12]]

w_lims = [-15,15]
thr_lims = [-14,0]
dbz_lims = [0,70]

fnums = [13, 28, 43, 58, 73]

fig1,axs1 = plt.subplots(3, 5, figsize=(15,8), subplot_kw=dict(box_aspect=1), layout='constrained')
fig2,axs2 = plt.subplots(3, 5, figsize=(15,8), subplot_kw=dict(box_aspect=1), layout='constrained')
fig3,axs3 = plt.subplots(3, 5, figsize=(15,8), subplot_kw=dict(box_aspect=1), layout='constrained')

# MERGER data
for i in np.arange(0,5):
    if i == 0:
        fp = '/Volumes/Promise_Pegasus_70TB/merger/merger-125m/base/'
    else:
        fp = '/Volumes/Promise_Pegasus_70TB/merger/merger-125m/'
    
    print(f"MERGER - cm1out_{fnums[i]:06d}")
    
    ds = nc.Dataset(fp+f"cm1out_{fnums[i]:06d}.nc")
    time = ds.variables['time'][:].data[0]/60
    xh = ds.variables['xh'][:].data
    xf = ds.variables['xf'][:].data
    yh = ds.variables['yh'][:].data
    yf = ds.variables['yf'][:].data
    z = ds.variables['z'][:].data
    zf = ds.variables['zf'][:].data
    
    xl = xlims[i]
    yl = ylims[i]
    ixw = np.where(xh >= xl[0])[0][0] # west bound
    ixe = np.where(xh >= xl[1])[0][0] # east bound
    iys = np.where(yh >= yl[0])[0][0] # south bound
    iyn = np.where(yh >= yl[1])[0][0] # north bound
    ix = slice(ixw,ixe+1)
    iy = slice(iys,iyn+1)

    iz1 = np.where(z >= 1)[0][0]
    iz = slice(0,iz1+1)
    
    qix = int(np.round(len(xh[(xh>=xl[0]) & (xh<=xl[1])])/60)*2)
    
    if i == 0:
        dbz = ds.variables['dbz2'][:].data[0,0,iy,ix]
        thr = ds.variables['th'][:].data[0,0,iy,ix] * (1 + 0.61*ds.variables['qv'][:].data[0,0,iy,ix] - 
                    (ds.variables['qc'][:].data[0,0,iy,ix] + ds.variables['qr'][:].data[0,0,iy,ix] + 
                     ds.variables['qi'][:].data[0,0,iy,ix] + ds.variables['qs'][:].data[0,0,iy,ix] + 
                     ds.variables['qg'][:].data[0,0,iy,ix] + ds.variables['qhl'][:].data[0,0,iy,ix]))
        thr0 = ds.variables['th0'][:].data[0,0,iy,ix] * (1 + 0.61*ds.variables['qv0'][:].data[0,0,iy,ix])
        thrpert = thr - thr0
        del thr,thr0
    else:
        dbz = ds.variables['dbz'][:].data[0,0,iy,ix]
    uinterp = ds.variables['uinterp'][:].data[0,iz,iy,ix]
    vinterp = ds.variables['vinterp'][:].data[0,iz,iy,ix]
    winterp = ds.variables['winterp'][:].data[0,iz,iy,ix]
    zvort = ds.variables['zvort'][:].data[0,iz,iy,ix]
    S_N = np.gradient(uinterp, xh[ix]*1000, axis=2) - np.gradient(vinterp, yh[iy]*1000, axis=1)
    S_S = np.gradient(vinterp, xh[ix]*1000, axis=2) + np.gradient(uinterp, yh[iy]*1000, axis=1)
    OW = S_N**2 + S_S**2 - zvort**2
    del S_N,S_S,zvort
    ds.close()
    
    if i > 0:
        ds = nc.Dataset(fp+f"pp/dyn_{fnums[i]:06d}.nc")
        thrpert = ds.variables['thrpert'][:].data[0,iy,ix]
        ds.close()
    
    
    c1 = plot_cfill(xh[ix], yh[iy], np.ma.masked_array(dbz, dbz<1), 'dbz', axs1[0,i], datalims=dbz_lims, xlims=xl, ylims=yl, cbar=False)
    axs1[0,i].contour(xh[ix], yh[iy], np.max(winterp,axis=0), levels=[5], colors='k', linewidths=1)
    # axs1[0,i].quiver(xh[ix][::qix], yh[iy][::qix], uinterp[0,::qix,::qix]+6, vinterp[0,::qix,::qix], color='k', scale=600, width=0.002, pivot='tail')
    axs1[0,i].set_title(f"t={time:.0f} min", fontsize=14)
    
    c2 = plot_cfill(xh[ix], yh[iy], thrpert, 'thrpert', axs2[0,i], datalims=thr_lims, cmap='YlGnBu_r', xlims=xl, ylims=yl, cbar=False)
    axs2[0,i].contour(xh[ix], yh[iy], np.max(winterp,axis=0), levels=[5], colors='k', linewidths=1, linestyles='-')
    axs2[0,i].contour(xh[ix], yh[iy], np.ma.masked_array(np.min(OW,axis=0), np.max(winterp,axis=0)<5), levels=[-0.001], colors='hotpink', linewidths=1.25, linestyles='-')
    axs2[0,i].quiver(xh[ix][::qix], yh[iy][::qix], uinterp[0,::qix,::qix]+6, vinterp[0,::qix,::qix], color='k', scale=600, width=0.002, pivot='tail')
    axs2[0,i].set_title(f"t={time:.0f} min", fontsize=14)
    
    c3 = plot_cfill(xh[ix], yh[iy], winterp[iz1,:,:], 'w', axs3[0,i], datalims=w_lims, xlims=xl, ylims=yl, cbar=False)
    axs3[0,i].contour(xh[ix], yh[iy], thrpert, levels=[-5], colors='g', linewidths=1, linestyles='-')
    axs3[0,i].contour(xh[ix], yh[iy], np.ma.masked_array(np.min(OW,axis=0), np.max(winterp,axis=0)<5), levels=[-0.001], colors='k', linewidths=1.25, linestyles='-')
    axs3[0,i].quiver(xh[ix][::qix], yh[iy][::qix], uinterp[iz1,::qix,::qix]+6, vinterp[iz1,::qix,::qix], color='k', scale=600, width=0.002, pivot='tail')
    axs3[0,i].set_title(f"t={time:.0f} min", fontsize=14)



# QLCS data
for i in np.arange(0,5):
    if i == 0:
        fp = '/Volumes/Promise_Pegasus_70TB/merger/qlcs-125m/base/'
    else:
        fp = '/Volumes/Promise_Pegasus_70TB/merger/qlcs-125m/'
    
    print(f"QLCS - cm1out_{fnums[i]:06d}")
    
    ds = nc.Dataset(fp+f"cm1out_{fnums[i]:06d}.nc")
    time = ds.variables['time'][:].data[0]/60
    xh = ds.variables['xh'][:].data
    xf = ds.variables['xf'][:].data
    yh = ds.variables['yh'][:].data
    yf = ds.variables['yf'][:].data
    z = ds.variables['z'][:].data
    zf = ds.variables['zf'][:].data
    
    xl = xlims[i]
    yl = ylims[i]
    ixw = np.where(xh >= xl[0])[0][0] # west bound
    ixe = np.where(xh >= xl[1])[0][0] # east bound
    iys = np.where(yh >= yl[0])[0][0] # south bound
    iyn = np.where(yh >= yl[1])[0][0] # north bound
    ix = slice(ixw,ixe+1)
    iy = slice(iys,iyn+1)

    iz1 = np.where(z >= 1)[0][0]
    iz = slice(0,iz1+1)
    
    if i == 0:
        dbz = ds.variables['dbz2'][:].data[0,0,iy,ix]
        thr = ds.variables['th'][:].data[0,0,iy,ix] * (1 + 0.61*ds.variables['qv'][:].data[0,0,iy,ix] - 
                    (ds.variables['qc'][:].data[0,0,iy,ix] + ds.variables['qr'][:].data[0,0,iy,ix] + 
                     ds.variables['qi'][:].data[0,0,iy,ix] + ds.variables['qs'][:].data[0,0,iy,ix] + 
                     ds.variables['qg'][:].data[0,0,iy,ix] + ds.variables['qhl'][:].data[0,0,iy,ix]))
        thr0 = ds.variables['th0'][:].data[0,0,iy,ix] * (1 + 0.61*ds.variables['qv0'][:].data[0,0,iy,ix])
        thrpert = thr - thr0
        del thr,thr0
    else:
        dbz = ds.variables['dbz'][:].data[0,0,iy,ix]
    uinterp = ds.variables['uinterp'][:].data[0,iz,iy,ix]
    vinterp = ds.variables['vinterp'][:].data[0,iz,iy,ix]
    winterp = ds.variables['winterp'][:].data[0,iz,iy,ix]
    zvort = ds.variables['zvort'][:].data[0,iz,iy,ix]
    S_N = np.gradient(uinterp, xh[ix]*1000, axis=2) - np.gradient(vinterp, yh[iy]*1000, axis=1)
    S_S = np.gradient(vinterp, xh[ix]*1000, axis=2) + np.gradient(uinterp, yh[iy]*1000, axis=1)
    OW = S_N**2 + S_S**2 - zvort**2
    del S_N,S_S,zvort
    ds.close()
    
    if i > 0:
        ds = nc.Dataset(fp+f"pp/dyn_{fnums[i]:06d}.nc")
        thrpert = ds.variables['thrpert'][:].data[0,iy,ix]
        ds.close()
    
    
    plot_cfill(xh[ix], yh[iy], np.ma.masked_array(dbz, dbz<1), 'dbz', axs1[1,i], datalims=dbz_lims, xlims=xl, ylims=yl, cbar=False)
    axs1[1,i].contour(xh[ix], yh[iy], np.max(winterp,axis=0), levels=[5], colors='k', linewidths=1)
    # axs1[1,i].quiver(xh[ix][::qix], yh[iy][::qix], uinterp[0,::qix,::qix]+6, vinterp[0,::qix,::qix], color='k', scale=600, width=0.002, pivot='tail')
    # axs1[1,i].set_title(f"t={time:.0f} min", fontsize=14)
    
    plot_cfill(xh[ix], yh[iy], thrpert, 'thrpert', axs2[1,i], datalims=thr_lims, cmap='YlGnBu_r', xlims=xl, ylims=yl, cbar=False)
    axs2[1,i].contour(xh[ix], yh[iy], np.max(winterp,axis=0), levels=[5], colors='k', linewidths=1, linestyles='-')
    axs2[1,i].contour(xh[ix], yh[iy], np.ma.masked_array(np.min(OW,axis=0), np.max(winterp,axis=0)<5), levels=[-0.001], colors='hotpink', linewidths=1.25, linestyles='-')
    axs2[1,i].quiver(xh[ix][::qix], yh[iy][::qix], uinterp[0,::qix,::qix]+6, vinterp[0,::qix,::qix], color='k', scale=600, width=0.002, pivot='tail')
    # axs2[1,i].set_title(f"t={time:.0f} min", fontsize=14)
    
    plot_cfill(xh[ix], yh[iy], winterp[iz1,:,:], 'w', axs3[1,i], datalims=w_lims, xlims=xl, ylims=yl, cbar=False)
    axs3[1,i].contour(xh[ix], yh[iy], thrpert, levels=[-5], colors='g', linewidths=1, linestyles='-')
    axs3[1,i].contour(xh[ix], yh[iy], np.ma.masked_array(np.min(OW,axis=0), np.max(winterp,axis=0)<5), levels=[-0.001], colors='k', linewidths=1.25, linestyles='-')
    axs3[1,i].quiver(xh[ix][::qix], yh[iy][::qix], uinterp[iz1,::qix,::qix]+6, vinterp[iz1,::qix,::qix], color='k', scale=600, width=0.002, pivot='tail')
    # axs3[1,i].set_title(f"t={time:.0f} min", fontsize=14)



# SUPERCELL data
for i in np.arange(0,5):
    if i == 0:
        fp = '/Volumes/Promise_Pegasus_70TB/merger/supercell-125m/base/'
    else:
        fp = '/Volumes/Promise_Pegasus_70TB/merger/supercell-125m/'
    
    print(f"SUPERCELL - cm1out_{fnums[i]:06d}")
    thr_lims = [-10,0]
    
    ds = nc.Dataset(fp+f"cm1out_{fnums[i]:06d}.nc")
    time = ds.variables['time'][:].data[0]/60
    xh = ds.variables['xh'][:].data
    xf = ds.variables['xf'][:].data
    yh = ds.variables['yh'][:].data
    yf = ds.variables['yf'][:].data
    z = ds.variables['z'][:].data
    zf = ds.variables['zf'][:].data
    
    xl = xlims[i]
    yl = ylims[i]
    ixw = np.where(xh >= xl[0])[0][0] # west bound
    ixe = np.where(xh >= xl[1])[0][0] # east bound
    iys = np.where(yh >= yl[0])[0][0] # south bound
    iyn = np.where(yh >= yl[1])[0][0] # north bound
    ix = slice(ixw,ixe+1)
    iy = slice(iys,iyn+1)

    iz1 = np.where(z >= 1)[0][0]
    iz = slice(0,iz1+1)
    
    if i == 0:
        dbz = ds.variables['dbz2'][:].data[0,0,iy,ix]
        thr = ds.variables['th'][:].data[0,0,iy,ix] * (1 + 0.61*ds.variables['qv'][:].data[0,0,iy,ix] - 
                    (ds.variables['qc'][:].data[0,0,iy,ix] + ds.variables['qr'][:].data[0,0,iy,ix] + 
                     ds.variables['qi'][:].data[0,0,iy,ix] + ds.variables['qs'][:].data[0,0,iy,ix] + 
                     ds.variables['qg'][:].data[0,0,iy,ix] + ds.variables['qhl'][:].data[0,0,iy,ix]))
        thr0 = ds.variables['th0'][:].data[0,0,iy,ix] * (1 + 0.61*ds.variables['qv0'][:].data[0,0,iy,ix])
        thrpert = thr - thr0
        del thr,thr0
    else:
        dbz = ds.variables['dbz'][:].data[0,0,iy,ix]
    uinterp = ds.variables['uinterp'][:].data[0,iz,iy,ix]
    vinterp = ds.variables['vinterp'][:].data[0,iz,iy,ix]
    winterp = ds.variables['winterp'][:].data[0,iz,iy,ix]
    zvort = ds.variables['zvort'][:].data[0,iz,iy,ix]
    S_N = np.gradient(uinterp, xh[ix]*1000, axis=2) - np.gradient(vinterp, yh[iy]*1000, axis=1)
    S_S = np.gradient(vinterp, xh[ix]*1000, axis=2) + np.gradient(uinterp, yh[iy]*1000, axis=1)
    OW = S_N**2 + S_S**2 - zvort**2
    del S_N,S_S,zvort
    ds.close()
    
    if i > 0:
        ds = nc.Dataset(fp+f"pp/dyn_{fnums[i]:06d}.nc")
        thrpert = ds.variables['thrpert'][:].data[0,iy,ix]
        ds.close()
    
    
    plot_cfill(xh[ix], yh[iy], np.ma.masked_array(dbz, dbz<1), 'dbz', axs1[2,i], datalims=dbz_lims, xlims=xl, ylims=yl, cbar=False)
    axs1[2,i].contour(xh[ix], yh[iy], np.max(winterp,axis=0), levels=[5], colors='k', linewidths=1)
    # axs1[2,i].quiver(xh[ix][::qix], yh[iy][::qix], uinterp[0,::qix,::qix]+6, vinterp[0,::qix,::qix], color='k', scale=600, width=0.002, pivot='tail')
    # axs1[2,i].set_title(f"t={time:.0f} min", fontsize=14)
    
    if i == 4:
        cb_flag = True
    else:
        cb_flag = False
    plot_cfill(xh[ix], yh[iy], thrpert, 'thrpert', axs2[2,i], datalims=thr_lims, cmap='YlGnBu_r', xlims=xl, ylims=yl, cbar=cb_flag)
    axs2[2,i].contour(xh[ix], yh[iy], np.max(winterp,axis=0), levels=[5], colors='k', linewidths=1, linestyles='-')
    axs2[2,i].contour(xh[ix], yh[iy], np.ma.masked_array(np.min(OW,axis=0), np.max(winterp,axis=0)<5), levels=[-0.001], colors='hotpink', linewidths=1.25, linestyles='-')
    axs2[2,i].quiver(xh[ix][::qix], yh[iy][::qix], uinterp[0,::qix,::qix]+6, vinterp[0,::qix,::qix], color='k', scale=600, width=0.002, pivot='tail')
    # axs2[2,i].set_title(f"t={time:.0f} min", fontsize=14)
    
    plot_cfill(xh[ix], yh[iy], winterp[iz1,:,:], 'w', axs3[2,i], datalims=w_lims, xlims=xl, ylims=yl, cbar=False)
    axs3[2,i].contour(xh[ix], yh[iy], thrpert, levels=[-5], colors='g', linewidths=1, linestyles='-')
    axs3[2,i].contour(xh[ix], yh[iy], np.ma.masked_array(np.min(OW,axis=0), np.max(winterp,axis=0)<5), levels=[-0.001], colors='k', linewidths=1.25, linestyles='-')
    axs3[2,i].quiver(xh[ix][::qix], yh[iy][::qix], uinterp[iz1,::qix,::qix]+6, vinterp[iz1,::qix,::qix], color='k', scale=600, width=0.002, pivot='tail')
    # axs3[2,i].set_title(f"t={time:.0f} min", fontsize=14)



cb1 = plt.colorbar(c1, ax=axs1[:,4], extend='min')
cb1.set_label('DBZ', fontsize=18)
cb2 = plt.colorbar(c2, ax=axs2[0:2,4], extend='min')
cb2.set_label("\u03B8'\u1D68 (K)", fontsize=18)
cb3 = plt.colorbar(c3, ax=axs3[:,4], extend='min')
cb3.set_label("w (m s$^{-1}$)", fontsize=18)


#%%

fig1.savefig('/Users/morgan.schneider/Documents/merger/overview_dbz.png', dpi=300)
fig2.savefig('/Users/morgan.schneider/Documents/merger/overview_thrpert.png', dpi=300)
fig3.savefig('/Users/morgan.schneider/Documents/merger/overview_w.png', dpi=300)



#%% Make overview plots (dBZ, w, thrpert, zvort)
######################
ip = '/Users/morgan.schneider/Documents/merger/merger-125m/'

wlims = [-15,15]
thr_lims = [-15,15]
pp_lims = [-20,20]
vort_lims = [-0.1,0.1]
dbz_lims = [1,70]



######################


xlims = [-180,180]
ylims = [-180,180]
qix = 35 #10

levs = [0, 0.1, 0.5, 1]
figsave = False

# for lev in levs:
#     iz = np.where(z >= lev)[0][0]
    
#     fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(12,9))
    
#     plot_cfill(xf, yf, dbz[iz,:,:], 'dbz', ax1, datalims=dbz_lims, xlims=xlims, ylims=ylims)
#     #ax1.quiver(xh[::qix], yh[::qix], uinterp[0,::qix,::qix], vinterp,0,::qix,::qix], scale=275, width=0.004, pivot='middle')
#     ax1.set_ylabel('y distance (km)')
#     ax1.set_title(f"$Z_H$ at {1000*z[iz]:.0f} m", fontsize=14)
    
#     plot_cfill(xf, yf, winterp[iz,:,:], 'w', ax2, datalims=wlims, xlims=xlims, ylims=ylims)
#     #ax2.quiver(xh[::qix], yh[::qix], uinterp[0,::qix,::qix], vinterp[0,::qix,::qix], scale=275, width=0.004, pivot='middle')
#     ax2.set_title(f"$w$ at {1000*z[iz]:.0f} m", fontsize=14)
    
#     plot_cfill(xf, yf, thrpert[iz,:,:], 'thrpert', ax3, datalims=thr_lims, xlims=xlims, ylims=ylims)
#     #ax3.quiver(xh[::qix], yh[::qix], uinterp[0,::qix,::qix], vinterp[0,::qix,::qix], scale=275, width=0.004, pivot='middle')
#     ax3.set_title(f"\u03B8\u1D68' at {1000*z[iz]:.0f} m", fontsize=14)
#     ax3.set_xlabel('x distance (km)')
#     ax3.set_ylabel('y distance (km)')
    
#     plot_cfill(xf, yf, zvort[iz,:,:], 'zvort', ax4, datalims=vort_lims, xlims=xlims, ylims=ylims)
#     #ax4.quiver(xh[::qix], yh[::qix], uinterp[0,::qix,::qix], vinterp[0,::qix,::qix], scale=275, width=0.004, pivot='middle')
#     ax4.set_xlabel('x distance (km)')
#     ax4.set_title(f"\u03B6 at {1000*z[iz]:.0f} m", fontsize=14)
    
#     plt.suptitle(f"$t={time:.0f}$ min", fontsize=14)
    
#     if figsave:
#         plt.savefig(ip+f"imgs_dbz/dbz_{1000*z[iz]:.0f}m_{time:03.0f}min.png", dpi=400)


xlims = [x1, x2]
ylims = [y1, y2]
iz1 = np.where(z>=1)[0][0]
qix = 6

# shade dbz, thrpert, or zvort?
fig,ax = plt.subplots(1,1,figsize=(8,6))

plot_cfill(xf, yf, thrpert[0,:,:], 'thrpert', ax, datalims=thr_lims, xlims=xlims, ylims=ylims)
ax.contour(xh, yh, dbz[0,:,:], levels=[10], colors='darkgray', linewidths=1, linestyles='--')
ax.contour(xh, yh, zvort[0,:,:], levels=[0.02,0.04,0.06,0.08,0.1], colors=['gold','darkorange','red','darkorchid','blue'], linewidths=1)
ax.contour(xh, yh, np.max(w[0:iz1,:,:], axis=0), levels=[10,20,30], colors='k', linewidths=1)
ax.contour(xh, yh, np.min(w[0:iz1,:,:], axis=0), levels=[-30,-20,-10], colors='dodgerblue', linewidths=1, linestyles='--')
ax.quiver(xh[::qix], yh[::qix], uinterp[0,::qix,::qix], vinterp[0,::qix,::qix], scale=275, width=0.004, pivot='middle')
ax.set_xlabel('x (km)')
ax.set_ylabel('y (km)')
ax.set_title(f"$t={time:.0f}$ min", fontsize=14)

if figsave:
    plt.savefig(ip+f"imgs_dbz/sfc_{time:03.0f}min.png", dpi=400)




#%% Load data for zoomed overview plots

from CM1utils import *

fp = '/Volumes/Promise_Pegasus_70TB/merger/supercell-125m/'
ip = '/Users/morgan.schneider/Documents/merger/supercell-125m/'

fn = 73
# 181 min -> 14 | 195 min -> 28 | 210 min -> 43 | 225 min -> 58 | 240 min -> 73

xlims = [-40,20] #[-55,25]
ylims = [-50,10] #[-100,-20]

# Read output file
ds = nc.Dataset(fp+f"cm1out_{fn:06d}.nc")
time = ds.variables['time'][:].data[0]/60
xh = ds.variables['xh'][:].data
xf = ds.variables['xf'][:].data
yh = ds.variables['yh'][:].data
yf = ds.variables['yf'][:].data
z = ds.variables['z'][:].data
zf = ds.variables['zf'][:].data

ixw = np.where(xh >= xlims[0])[0][0] # west bound
ixe = np.where(xh >= xlims[1])[0][0] # east bound
iys = np.where(yh >= ylims[0])[0][0] # south bound
iyn = np.where(yh >= ylims[1])[0][0] # north bound
ix = slice(ixw,ixe)
iy = slice(iys,iyn)

iz1 = np.where(z >= 1)[0][0]
# iz15 = np.where(z >= 1.5)[0][0]
# iz2 = np.where(z >= 2)[0][0]

if 'dbz2' in ds.variables:
    dbz_sfc = ds.variables['dbz2'][:].data[0,0,iy,ix]
    # dbz_1km = ds.variables['dbz2'][:].data[0,iz1,iy,ix]
    thr = ds.variables['th'][:].data[0,0,iy,ix] * (1 + 0.61*ds.variables['qv'][:].data[0,0,iy,ix] - 
                (ds.variables['qc'][:].data[0,0,iy,ix] + ds.variables['qr'][:].data[0,0,iy,ix] + 
                 ds.variables['qi'][:].data[0,0,iy,ix] + ds.variables['qs'][:].data[0,0,iy,ix] + 
                 ds.variables['qg'][:].data[0,0,iy,ix] + ds.variables['qhl'][:].data[0,0,iy,ix]))
    thr0 = ds.variables['th0'][:].data[0,0,iy,ix] * (1 + 0.61*ds.variables['qv0'][:].data[0,0,iy,ix])
    thrpert = thr - thr0
    del thr,thr0
else:
    dbz_sfc = ds.variables['dbz'][:].data[0,0,iy,ix]
    # dbz_1km = ds.variables['dbz'][:].data[0,iz1,iy,ix]
uinterp = ds.variables['uinterp'][:].data[0,slice(0,iz1+1),iy,ix]
vinterp = ds.variables['vinterp'][:].data[0,slice(0,iz1+1),iy,ix]
winterp = ds.variables['winterp'][:].data[0,slice(0,iz1+1),iy,ix]
# yvort = ds.variables['yvort'][:].data[0,slice(0,iz1+1),iy,ix]
zvort = ds.variables['zvort'][:].data[0,slice(0,iz1+1),iy,ix]
S_N = np.gradient(uinterp, xh[ix]*1000, axis=2) - np.gradient(vinterp, yh[iy]*1000, axis=1)
S_S = np.gradient(vinterp, xh[ix]*1000, axis=2) + np.gradient(uinterp, yh[iy]*1000, axis=1)
OW = S_N**2 + S_S**2 - zvort**2
del S_N,S_S,zvort,uinterp,vinterp

ds.close()


# Read dyn
ds = nc.Dataset(fp+f"pp/dyn_{fn:06d}.nc")
# p_dn = ds.variables['p_dn'][:].data[slice(0,iz1+1),iy,ix]
# p_b = ds.variables['p_b'][:].data[slice(0,iz1+1),iy,ix]
# p_dl = ds.variables['p_dl'][:].data[slice(0,iz1+1),iy,ix]
thrpert = ds.variables['thrpert'][:].data[0,iy,ix]
ds.close()


# iz200 = np.where(z >= 0.2)[0][0]

# vppga_dn = -1/1.1 * (p_dn[iz1,:,:] - p_dn[0,:,:]) / (1000 * (z[iz1]-z[0]))
# vppga_b = -1/1.1 * (p_b[iz1,:,:] - p_b[0,:,:]) / (1000 * (z[iz1]-z[0]))
# vppga_dl = -1/1.1 * (p_dl[iz1,:,:] - p_dl[0,:,:]) / (1000 * (z[iz1]-z[0]))

# nonlinear dynamic can be funky at sfc because of residual and buoyancy has a fucked up bottom bc
# use 150 or 200 m as lower bound instead
# vppga_dn_200 = -1/1.1 * (p_dn[iz1,:,:] - p_dn[iz200,:,:]) / (1000 * (z[iz1]-z[iz200]))
# vppga_b_200 = -1/1.1 * (p_b[iz1,:,:] - p_b[iz200,:,:]) / (1000 * (z[iz1]-z[iz200]))
# vppga_dl_200 = -1/1.1 * (p_dl[iz1,:,:] - p_dl[iz200,:,:]) / (1000 * (z[iz1]-z[iz200])




#%% Make zoomed overview plots (dBZ, w, thrpert, zvort)

w_lims = [-15,15]
thr_lims = [-14,0]
pp_lims = [-2,2]
vort_lims = [-0.1,0.1]
dbz_lims = [0,70]


xl = [-16.5,-6.5]
yl = [-9.5,0.5]
# xl = xlims
# yl = ylims

figsave = False



# qix = 8
qix = int(np.round(len(xh[(xh>=xl[0]) & (xh<=xl[1])])/60)*2)


# dbz, w and OW contoured
if True:
    # TLV criteria - dbz, 1km w and sfc OW contoured
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    
    plot_cfill(xh[ix], yh[iy], np.ma.masked_array(dbz_sfc, dbz_sfc<1), 'dbz', ax, datalims=dbz_lims)
    # ax.contour(xh[ix], yh[iy], thrpert, levels=[-5], colors='w', linewidths=1, linestyles='-')
    ax.contour(xh[ix], yh[iy], np.max(winterp[0:iz1,:,:],axis=0), levels=[5], colors='k', linewidths=1)
    # ax.contour(xh[ix], yh[iy], zvort[0,:,:], levels=[0.01,0.05,0.1], colors='k', linewidths=1)
    # ax.contour(xh[ix], yh[iy], np.min(OW[0:iz1,:,:],axis=0), levels=[-0.001], colors='silver', linewidths=1, linestyles='-')
    # ax.quiver(xh[ix][::qix], yh[iy][::qix], uinterp[0,::qix,::qix]+6, vinterp[0,::qix,::qix], color='k', scale=600, width=0.002, pivot='tail')
    ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')
    ax.set_title(f"$Z_{{sfc}}$, $w$ (black) \n t={time:.0f} min", fontsize=14)
    ax.set_xlim(xl)
    ax.set_ylim(yl)
    if figsave:
        plt.savefig(ip+f"imgs_overview/dbz_sfc_{time:03.0f}min.png", dpi=300)
    
    # MV criteria - 1 km w, 1 km OW
    # fig,ax = plt.subplots(1,1,figsize=(8,6))
    
    # plot_cfill(xh[ix], yh[iy], np.ma.masked_array(dbz_1km, dbz_1km<1), 'dbz', ax, datalims=dbz_lims)
    # ax.contour(xh[ix], yh[iy], winterp[iz1,:,:], levels=[-15,-10,10,15], colors='silver', linewidths=1)
    # # ax.contour(xh[ix], yh[iy], zvort[iz1,:,:], levels=[0.01,0.05,0.1], colors='k', linewidths=1)
    # ax.contour(xh[ix], yh[iy], OW[iz1,:,:], levels=[-0.001], colors='k', linewidths=1, linestyles='-')
    # # ax.quiver(xh[ix][::qix], yh[iy][::qix], uinterp[0,::qix,::qix], vinterp[0,::qix,::qix], scale=450, width=0.003, pivot='middle')
    # ax.set_xlabel('x (km)')
    # ax.set_ylabel('y (km)')
    # ax.set_xlim(xl)
    # ax.set_ylim(yl)
    # ax.set_title(f"$Z_{{1km}}$ (shaded), OW$_{{1km}}$ (black), $w_{{1km}}$ (gray)\n t={time:.0f} min", fontsize=14)
    
    # if figsave:
    #     plt.savefig(ip+f"1km_dbz_{time:03.0f}min.png", dpi=300)


# w, OW and thrpert contoured
# --> shade zvort to check if MVs are cyclonic or couplets! for formation mechanisms
if True:
    # wspd = ((uinterp[iz1,:,:]+6)**2 + vinterp[iz1,:,:]**2)**0.5
    
    # MV criteria - 1 km w, 1 km OW and 40 dBZ contoured
    # fig,ax = plt.subplots(1,1,figsize=(8,6))
    # plot_cfill(xh[ix], yh[iy], winterp[iz1,:,:], 'w', ax, datalims=w_lims)
    # ax.contour(xh[ix], yh[iy], thrpert, levels=[-5], colors='g', linewidths=1, linestyles='-')
    # ax.contour(xh[ix], yh[iy], np.ma.masked_array(np.min(OW[0:iz1,:,:],axis=0), winterp[iz1,:,:]<5), levels=[-0.001], colors='k', linewidths=1.25, linestyles='-')
    # # ax.quiver(xh[ix][::qix], yh[iy][::qix], uinterp[iz1,::qix,::qix]+6, vinterp[iz1,::qix,::qix], color='k', scale=600, width=0.002, pivot='tail')
    # ax.set_xlabel('x (km)')
    # ax.set_ylabel('y (km)')
    # ax.set_xlim(xl)
    # ax.set_ylim(yl)
    # ax.set_title(f"$w_{{1km}}$, OW (black), \u03B8'\u1D68 (green) \n t={time:.0f} min", fontsize=14)
    # if figsave:
    #     plt.savefig(ip+f"imgs_overview/w_1km_{time:03.0f}min.png", dpi=300)
    
    
    # wspd = ((uinterp[0,:,:])**2 + vinterp[0,:,:]**2)**0.5
    
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    plot_cfill(xh[ix], yh[iy], thrpert, 'thrpert', ax, datalims=thr_lims, cmap='YlGnBu_r')
    ax.contour(xh[ix], yh[iy], np.max(winterp[0:iz1,:,:],axis=0), levels=[5], colors='k', linewidths=1, linestyles='-')
    ax.contour(xh[ix], yh[iy], np.ma.masked_array(np.min(OW[0:iz1,:,:],axis=0), winterp[iz1,:,:]<5), levels=[-0.001], colors='hotpink', linewidths=1.25, linestyles='-')
    ax.quiver(xh[ix][::qix], yh[iy][::qix], uinterp[0,::qix,::qix]+6, vinterp[0,::qix,::qix], color='k', scale=600, width=0.002, pivot='tail')
    ax.quiver(xh[ix][::qix], yh[iy][::qix], uinterp[iz1,::qix,::qix]+6, vinterp[iz1,::qix,::qix], color='gray', scale=600, width=0.002, pivot='tail')
    ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')
    ax.set_xlim(xl)
    ax.set_ylim(yl)
    ax.set_title(f"\u03B8'\u1D68, $w$ (black), OW (pink) \n t={time:.0f} min", fontsize=14)
    if figsave:
        plt.savefig(ip+f"imgs_overview/thrpert_sfc_{time:03.0f}min.png", dpi=300)

# 200m p'dn, 1km p'dn, 200-1000m vppga dn
if False:
    # Sfc and 0-1 km mean p' dn
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    plot_cfill(xh[ix], yh[iy], p_dn[iz200,:,:]/100, 'prspert', ax, datalims=pp_lims)
    # ax.contour(xh[ix], yh[iy], thrpert, levels=[-8], colors='slategray', linewidths=1, linestyles='-')
    ax.contour(xh[ix], yh[iy], winterp[iz1,:,:], levels=[-15,-10,10,15], colors='yellow', linewidths=1)
    ax.contour(xh[ix], yh[iy], OW[0,:,:], levels=[-0.04,-0.0225,-0.01,-0.0025], colors='k', linewidths=1, linestyles='-')
    # ax.quiver(xh[ix][::qix], yh[iy][::qix], uinterp[0,::qix,::qix], vinterp[0,::qix,::qix], scale=275, width=0.003, pivot='middle')
    ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')
    ax.set_xlim(xl)
    ax.set_ylim(yl)
    ax.set_title(f"200 m $p_{{dn}}'$, OW$_{{sfc}}$ (black), $w_{{1km}}$ (yellow)\n t={time:.0f} min", fontsize=14)
    if figsave:
        plt.savefig(ip+f"imgs_tmp/sfc_pdn_{time:03.0f}min.png", dpi=300)
    
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    plot_cfill(xh[ix], yh[iy], p_dn[iz1,:,:]/100, 'prspert', ax, datalims=pp_lims)
    # ax.contour(xh[ix], yh[iy], thrpert, levels=[-8], colors='slategray', linewidths=1, linestyles='-')
    ax.contour(xh[ix], yh[iy], OW[iz1,:,:], levels=[-0.001], colors='k', linewidths=1, linestyles='-')
    ax.contour(xh[ix], yh[iy], winterp[iz1,:,:], levels=[-15,-10,10,15], colors='yellow', linewidths=1)
    # ax.contour(xh[ix], yh[iy], zvort[iz1,:,:], levels=[0.05,0.1], colors='red', linewidths=1)
    # ax.contour(xh[ix], yh[iy], yvort[iz1,:,:], levels=[-0.1,-0.05,0.05,0.1], colors='yellow', linewidths=1)
    # ax.quiver(xh[ix][::qix], yh[iy][::qix], uinterp[iz1,::qix,::qix], vinterp[iz1,::qix,::qix], scale=450, width=0.003, pivot='middle')
    ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')
    ax.set_xlim(xl)
    ax.set_ylim(yl)
    ax.set_title(f"1 km $p_{{dn}}'$, OW$_{{1km}}$ (black), $w_{{1km}}$ (yellow)\n t={time:.0f} min", fontsize=14)
    if figsave:
        plt.savefig(ip+f"imgs_tmp/1km_pdn_{time:03.0f}min.png", dpi=300)
    
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    plot_cfill(xh[ix], yh[iy], vppga_dn_200, 'pgfz', ax, datalims=[-0.2,0.2])
    # ax.contour(xh[ix], yh[iy], thrpert, levels=[-8], colors='slategray', linewidths=1, linestyles='-')
    ax.contour(xh[ix], yh[iy], OW[iz1,:,:], levels=[-0.001], colors='k', linewidths=1, linestyles='-')
    ax.contour(xh[ix], yh[iy], winterp[iz1,:,:], levels=[-15,-10,5,10,15], colors='yellow', linewidths=1)
    # ax.quiver(xh[ix][::qix], yh[iy][::qix], uinterp[iz1,::qix,::qix], vinterp[iz1,::qix,::qix], scale=450, width=0.003, pivot='middle')
    ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')
    ax.set_xlim(xl)
    ax.set_ylim(yl)
    ax.set_title(f"200-1000 m VPPGA$_{{dn}}$, OW$_{{1km}}$ (black), $w_{{1km}}$ (yellow)\n t={time:.0f} min", fontsize=14)
    if figsave:
        plt.savefig(ip+f"imgs_tmp/1km_vppga_dn_{time:03.0f}min.png", dpi=300)

# Sfc p', 1km p', 0-1km vppga
if False:
    fig,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(16,4))
    plot_cfill(xh[ix], yh[iy], p_dn[0,:,:]/100, 'prspert', ax1, datalims=pp_lims)
    # ax1.contour(xh[ix], yh[iy], dbz_sfc, levels=[40], colors='silver', linewidths=1, linestyles='-')
    ax1.contour(xh[ix], yh[iy], winterp[iz1,:,:], levels=[-15,-10,10,15], colors='yellow', linewidths=1)
    ax1.contour(xh[ix], yh[iy], OW[0,:,:], levels=[-0.04,-0.0225,-0.01,-0.0025], colors='k', linewidths=1, linestyles='-')
    # ax1.quiver(xh[ix][::qix], yh[iy][::qix], uinterp[0,::qix,::qix], vinterp[0,::qix,::qix], scale=275, width=0.003, pivot='middle')
    ax1.set_xlabel('x (km)')
    ax1.set_ylabel('y (km)')
    ax1.set_xlim(xl)
    ax1.set_ylim(yl)
    ax1.set_title(f"Sfc $p_{{dn}}'$, {time:.0f} min", fontsize=14)
    
    plot_cfill(xh[ix], yh[iy], p_dl[0,:,:]/100, 'prspert', ax2, datalims=pp_lims)
    # ax2.contour(xh[ix], yh[iy], dbz_sfc, levels=[40], colors='silver', linewidths=1, linestyles='-')
    ax2.contour(xh[ix], yh[iy], winterp[iz1,:,:], levels=[-15,-10,10,15], colors='yellow', linewidths=1)
    ax2.contour(xh[ix], yh[iy], OW[0,:,:], levels=[-0.04,-0.0225,-0.01,-0.0025], colors='k', linewidths=1, linestyles='-')
    # ax2.quiver(xh[ix][::qix], yh[iy][::qix], uinterp[0,::qix,::qix], vinterp[0,::qix,::qix], scale=275, width=0.003, pivot='middle')
    ax2.set_xlabel('x (km)')
    ax2.set_ylabel('y (km)')
    ax2.set_xlim(xl)
    ax2.set_ylim(yl)
    ax2.set_title(f"Sfc $p_{{dl}}'$, {time:.0f} min", fontsize=14)
    
    plot_cfill(xh[ix], yh[iy], p_b[0,:,:]/100, 'prspert', ax3, datalims=pp_lims)
    # ax3.contour(xh[ix], yh[iy], dbz_sfc, levels=[40], colors='silver', linewidths=1, linestyles='-')
    ax3.contour(xh[ix], yh[iy], winterp[iz1,:,:], levels=[-15,-10,10,15], colors='yellow', linewidths=1)
    ax3.contour(xh[ix], yh[iy], OW[0,:,:], levels=[-0.04,-0.0225,-0.01,-0.0025], colors='k', linewidths=1, linestyles='-')
    # ax3.quiver(xh[ix][::qix], yh[iy][::qix], uinterp[0,::qix,::qix], vinterp[0,::qix,::qix], scale=275, width=0.003, pivot='middle')
    ax3.set_xlabel('x (km)')
    ax3.set_ylabel('y (km)')
    ax3.set_xlim(xl)
    ax3.set_ylim(yl)
    ax3.set_title(f"Sfc $p_{{b}}'$, {time:.0f} min", fontsize=14)
    
    if figsave:
        plt.savefig(ip+f"imgs_tmp/sfc_pp_all_{time:03.0f}min.png", dpi=300)
    
    
    
    fig,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(16,4))
    plot_cfill(xh[ix], yh[iy], p_dn[iz1,:,:]/100, 'prspert', ax1, datalims=pp_lims)
    # ax1.contour(xh[ix], yh[iy], dbz_sfc, levels=[40], colors='silver', linewidths=1, linestyles='-')
    ax1.contour(xh[ix], yh[iy], winterp[iz1,:,:], levels=[-15,-10,10,15], colors='yellow', linewidths=1)
    ax1.contour(xh[ix], yh[iy], OW[iz1,:,:], levels=[-0.04,-0.0225,-0.01,-0.0025], colors='k', linewidths=1, linestyles='-')
    # ax1.quiver(xh[ix][::qix], yh[iy][::qix], uinterp[iz1,::qix,::qix], vinterp[iz1,::qix,::qix], scale=275, width=0.003, pivot='middle')
    ax1.set_xlabel('x (km)')
    ax1.set_ylabel('y (km)')
    ax1.set_xlim(xl)
    ax1.set_ylim(yl)
    ax1.set_title(f"1 km $p_{{dn}}'$, {time:.0f} min", fontsize=14)
    
    plot_cfill(xh[ix], yh[iy], p_dl[iz1,:,:]/100, 'prspert', ax2, datalims=pp_lims)
    # ax2.contour(xh[ix], yh[iy], dbz_sfc, levels=[40], colors='silver', linewidths=1, linestyles='-')
    ax2.contour(xh[ix], yh[iy], winterp[iz1,:,:], levels=[-15,-10,10,15], colors='yellow', linewidths=1)
    ax2.contour(xh[ix], yh[iy], OW[iz1,:,:], levels=[-0.04,-0.0225,-0.01,-0.0025], colors='k', linewidths=1, linestyles='-')
    # ax2.quiver(xh[ix][::qix], yh[iy][::qix], uinterp[iz1,::qix,::qix], vinterp[iz1,::qix,::qix], scale=275, width=0.003, pivot='middle')
    ax2.set_xlabel('x (km)')
    ax2.set_ylabel('y (km)')
    ax2.set_xlim(xl)
    ax2.set_ylim(yl)
    ax2.set_title(f"1 km $p_{{dl}}'$, {time:.0f} min", fontsize=14)
    
    plot_cfill(xh[ix], yh[iy], p_b[iz1,:,:]/100, 'prspert', ax3, datalims=pp_lims)
    # ax3.contour(xh[ix], yh[iy], dbz_sfc, levels=[40], colors='silver', linewidths=1, linestyles='-')
    ax3.contour(xh[ix], yh[iy], winterp[iz1,:,:], levels=[-15,-10,10,15], colors='yellow', linewidths=1)
    ax3.contour(xh[ix], yh[iy], OW[iz1,:,:], levels=[-0.04,-0.0225,-0.01,-0.0025], colors='k', linewidths=1, linestyles='-')
    # ax3.quiver(xh[ix][::qix], yh[iy][::qix], uinterp[iz1,::qix,::qix], vinterp[iz1,::qix,::qix], scale=275, width=0.003, pivot='middle')
    ax3.set_xlabel('x (km)')
    ax3.set_ylabel('y (km)')
    ax3.set_xlim(xl)
    ax3.set_ylim(yl)
    ax3.set_title(f"1 km $p_{{b}}'$, {time:.0f} min", fontsize=14)
    
    if figsave:
        plt.savefig(ip+f"imgs_tmp/1km_pp_all_{time:03.0f}min.png", dpi=300)
    
    
    
    fig,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(16,4))
    plot_cfill(xh[ix], yh[iy], vppga_dn, 'pgfz', ax1, datalims=[-0.2,0.2])
    # ax1.contour(xh[ix], yh[iy], dbz_sfc, levels=[40], colors='silver', linewidths=1, linestyles='-')
    ax1.contour(xh[ix], yh[iy], winterp[iz1,:,:], levels=[-15,-10,10,15], colors='yellow', linewidths=1)
    ax1.contour(xh[ix], yh[iy], OW[iz1,:,:], levels=[-0.04,-0.0225,-0.01,-0.0025], colors='k', linewidths=1, linestyles='-')
    # ax1.quiver(xh[ix][::qix], yh[iy][::qix], uinterp[iz1,::qix,::qix], vinterp[iz1,::qix,::qix], scale=275, width=0.003, pivot='middle')
    ax1.set_xlabel('x (km)')
    ax1.set_ylabel('y (km)')
    ax1.set_xlim(xl)
    ax1.set_ylim(yl)
    ax1.set_title(f"0-1 km VPPGA$_{{dn}}$, {time:.0f} min", fontsize=14)
    
    plot_cfill(xh[ix], yh[iy], vppga_dl, 'pgfz', ax2, datalims=[-0.2,0.2])
    # ax2.contour(xh[ix], yh[iy], dbz_sfc, levels=[40], colors='silver', linewidths=1, linestyles='-')
    ax2.contour(xh[ix], yh[iy], winterp[iz1,:,:], levels=[-15,-10,10,15], colors='yellow', linewidths=1)
    ax2.contour(xh[ix], yh[iy], OW[iz1,:,:], levels=[-0.04,-0.0225,-0.01,-0.0025], colors='k', linewidths=1, linestyles='-')
    # ax2.quiver(xh[ix][::qix], yh[iy][::qix], uinterp[iz1,::qix,::qix], vinterp[iz1,::qix,::qix], scale=275, width=0.003, pivot='middle')
    ax2.set_xlabel('x (km)')
    ax2.set_ylabel('y (km)')
    ax2.set_xlim(xl)
    ax2.set_ylim(yl)
    ax2.set_title(f"0-1 km VPPGA$_{{dl}}$, {time:.0f} min", fontsize=14)
    
    plot_cfill(xh[ix], yh[iy], vppga_b, 'pgfz', ax3, datalims=[-0.2,0.2])
    # ax3.contour(xh[ix], yh[iy], dbz_sfc, levels=[40], colors='silver', linewidths=1, linestyles='-')
    ax3.contour(xh[ix], yh[iy], winterp[iz1,:,:], levels=[-15,-10,10,15], colors='yellow', linewidths=1)
    ax3.contour(xh[ix], yh[iy], OW[iz1,:,:], levels=[-0.04,-0.0225,-0.01,-0.0025], colors='k', linewidths=1, linestyles='-')
    # ax3.quiver(xh[ix][::qix], yh[iy][::qix], uinterp[iz1,::qix,::qix], vinterp[iz1,::qix,::qix], scale=275, width=0.003, pivot='middle')
    ax3.set_xlabel('x (km)')
    ax3.set_ylabel('y (km)')
    ax3.set_xlim(xl)
    ax3.set_ylim(yl)
    ax3.set_title(f"0-1 km VPPGA$_{{b}}$, {time:.0f} min", fontsize=14)
    
    if figsave:
        plt.savefig(ip+f"imgs_tmp/1km_vppga_all_{time:03.0f}min.png", dpi=300)
    

if False:
    iy1 = np.where(yh[iy]>=-112.8)[0][0]
    iy1 = np.where(yh[iy]>=-109.5)[0][0]
    
    fig,ax = plt.subplots(1,1,figsize=(10,4))
    plot_cfill(xh[ix], z[0:iz1+1], winterp[:,iy1,:], 'w', ax, datalims=w_lims, xlims=[-41,-31], ylims=[0,1])
    # ax.contour(xh[ix], z[0:iz1+1], thrpert[:,iy1,:], levels=[-5], colors='g', linewidths=1, linestyles='-')
    # ax.contour(xh, z, winterp[:,iy1,:], levels=[-10,-5,5,10], colors='w', linewidths=1)
    # ax.contour(xh, z, np.mean(zvort[:,iy1-4:iy1+4,:], axis=1), levels=[0.025,0.05,0.075,0.1], colors='k', linewidths=1)
    ax.contour(xh[ix], z[0:iz1+1], np.min(OW[:,iy1-4:iy1+4,:], axis=1), levels=[-0.01,-0.0025,-0.001], colors='k', linewidths=1, linestyles='-')
    ax.quiver(xh[ix][::3], z[0:iz1+1][::1], uinterp[::1,iy1,::3]+6, winterp[::1,iy1,::3], scale=600, width=0.002, pivot='middle')
    ax.set_xlabel('x (km)')
    ax.set_ylabel('z (km)')
    # ax.set_title(f"$Z_H$ (shaded), OW (black), $w$ (white), \u03B8\u1D68'=-8 K (blue)\n t={time:.0f} min", fontsize=14)
    


#%% Zoomed-in plan view with shared-axis vertical cross sections

from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1 import make_axes_locatable

fp = '/Volumes/Promise_Pegasus_70TB/merger/merger-125m/'
ip = '/Users/morgan.schneider/Documents/merger/merger-125m/'

fn = 43
# 181 min -> 14 | 195 min -> 28 | 210 min -> 43 | 225 min -> 58 | 240 min -> 73

# Read output file
ds = nc.Dataset(fp+f"cm1out_{fn:06d}.nc")
time = ds.variables['time'][:].data[0]/60
xh = ds.variables['xh'][:].data
xf = ds.variables['xf'][:].data
yh = ds.variables['yh'][:].data
yf = ds.variables['yf'][:].data
z = ds.variables['z'][:].data
zf = ds.variables['zf'][:].data

xlims = [-50,10] #[-55,25]
ylims = [-90,-30] #[-100,-20]

ixw = np.where(xh >= xlims[0])[0][0] # west bound
ixe = np.where(xh >= xlims[1])[0][0] # east bound
iys = np.where(yh >= ylims[0])[0][0] # south bound
iyn = np.where(yh >= ylims[1])[0][0] # north bound
ix = slice(ixw,ixe+1)
iy = slice(iys,iyn+1)

iz1 = np.where(z >= 1)[0][0]
iz15 = np.where(z >= 1.5)[0][0]
iz2 = np.where(z >= 2)[0][0]
iz5 = np.where(z >= 5)[0][0]
iz200 = np.where(z >= 0.2)[0][0]
iz = slice(0,iz5+1)

if 'dbz2' in ds.variables:
    dbz = ds.variables['dbz2'][:].data[0,iz,iy,ix]
    thr = ds.variables['th'][:].data[0,iz,iy,ix] * (1 + 0.61*ds.variables['qv'][:].data[0,iz,iy,ix] - 
                (ds.variables['qc'][:].data[0,iz,iy,ix] + ds.variables['qr'][:].data[0,iz,iy,ix] + 
                 ds.variables['qi'][:].data[0,iz,iy,ix] + ds.variables['qs'][:].data[0,iz,iy,ix] + 
                 ds.variables['qg'][:].data[0,iz,iy,ix] + ds.variables['qhl'][:].data[0,iz,iy,ix]))
    thr0 = ds.variables['th0'][:].data[0,iz,iy,ix] * (1 + 0.61*ds.variables['qv0'][:].data[0,iz,iy,ix])
    thrpert = thr - thr0
    del thr,thr0
else:
    dbz = ds.variables['dbz'][:].data[0,iz,iy,ix]
# prspert = ds.variables['prspert'][:].data[0,iz,iy,ix]/100
uinterp = ds.variables['uinterp'][:].data[0,iz,iy,ix]
vinterp = ds.variables['vinterp'][:].data[0,iz,iy,ix]
winterp = ds.variables['winterp'][:].data[0,iz,iy,ix]
# zvort = ds.variables['zvort'][:].data[0,iz,iy,ix]
# S_N = np.gradient(uinterp, xh[ix]*1000, axis=2) - np.gradient(vinterp, yh[iy]*1000, axis=1)
# S_S = np.gradient(vinterp, xh[ix]*1000, axis=2) + np.gradient(uinterp, yh[iy]*1000, axis=1)
# OW = S_N**2 + S_S**2 - zvort**2
# del S_N,S_S,zvort
ds.close()

# Read dyn
ds = nc.Dataset(fp+f"pp/dyn_{fn:06d}.nc")
thrpert = ds.variables['thrpert'][:].data[iz,iy,ix]
ds.close()

img_str = 's1'


wlims = [-15,15]
thr_lims = [-15,15]
pp_lims = [-20,20]
vort_lims = [-0.1,0.1]
dbz_lims = [0,70]

# xlims = [-50,-30] # [-60,0]
# ylims = [-130,-110] # [-150,-90]

qix = int(np.round(len(xh[(xh>=xlims[0]) & (xh<=xlims[1])])/60)*2) # scale wind vector spacing to number of grid pts
qiz = 3


figsave = False



if False:
    iy1 = np.where(xh[ix] >= -17.8)[0][0]
    ix1 = np.where(yh[iy] >= -89.7)[0][0]

    xlims = [xh[ix][ix1-80], xh[ix][ix1+81]]
    ylims = [yh[iy][iy1-80], yh[iy][iy1+81]]
    
    fig,ax = plt.subplots(figsize=(11,8))
    
    divider = make_axes_locatable(ax)
    top_ax = divider.append_axes('top', 1.85, pad=0.2, sharex=ax)
    side_ax = divider.append_axes('left', 2, pad=0.2, sharey=ax)
    
    top_ax.xaxis.set_tick_params(labelbottom=False)
    ax.yaxis.set_tick_params(labelleft=False)
    
    ax.set_xlabel('x (km)', fontsize=14)
    side_ax.set_ylabel('y (km)', fontsize=14)
    top_ax.set_ylabel('z (km)', fontsize=14)
    side_ax.set_xlabel('z (km)', fontsize=14)
    
    
    plot_cfill(xh[ix], yh[iy], np.ma.masked_array(dbz[0,:,:], dbz[0,:,:]<1), 'dbz', ax, datalims=dbz_lims, xlims=xlims, ylims=ylims)
    # ax.contour(xh[ix], yh[iy], thrpert[0,:,:], levels=[-8], colors='blue', linewidths=1, linestyles='-')
    ax.contour(xh[ix], yh[iy], winterp[iz1,:,:], levels=[-10,-5,5,10], colors='silver', linewidths=1)
    # ax.contour(xh[ix], yh[iy], OW[iz1,:,:], levels=[-0.0001], colors='k', linewidths=1, linestyles='-')
    ax.contour(xh[ix], yh[iy], OW[iz1,:,:], levels=[-0.001], colors='k', linewidths=1, linestyles='-')
    ax.quiver(xh[ix][::qix], yh[iy][::qix], uinterp[0,::qix,::qix], vinterp[0,::qix,::qix], scale=600, width=0.002, pivot='tail')
    hline = ax.axhline(yh[iy], color='k', linewidth=0.75, linestyle='--')
    vline = ax.axvline(xh[ix], color='k', linewidth=0.75, linestyle='--')
    
    plot_cfill(xh[ix], z[iz], np.ma.masked_array(dbz[:,iy1,:], dbz[:,iy1,:]<1), 'dbz', top_ax, datalims=dbz_lims, xlims=xlims, ylims=[0,4])
    top_ax.contour(xh[ix], z[iz], np.mean(thrpert[:,iy1-4:iy1+4,:], axis=1), levels=[-8], colors='blue', linewidths=1, linestyles='-')
    top_ax.contour(xh[ix], z[iz], np.mean(winterp[:,iy1-4:iy1+4,:], axis=1), levels=[-15,-10,-5,5,10,15], colors='lightgrey', linewidths=1)
    top_ax.contour(xh[ix], z[iz], np.min(OW[:,iy1-4:iy1+4,:], axis=1), levels=[-0.01,-0.0025,-0.001], colors='k', linewidths=1, linestyles='-')
    # top_ax.quiver(xh[ix][::qix], z[iz][::qiz], uinterp[::qiz,iy1,::qix], winterp[::qiz,iy1,::qix], scale=400, width=0.002, pivot='tail')
    
    plot_cfill(z[iz], yh[iy], np.ma.masked_array(dbz[:,:,ix1], dbz[:,:,ix1]<1).transpose(), 'dbz', side_ax, datalims=dbz_lims, xlims=[0,4], ylims=ylims)
    side_ax.contour(z[iz], yh[iy], np.mean(thrpert[:,:,ix1-4:ix1+4], axis=2).transpose(), levels=[-8], colors='blue', linewidths=1, linestyles='-')
    side_ax.contour(z[iz], yh[iy], np.mean(winterp[:,:,ix1-4:ix1+4], axis=2).transpose(), levels=[-15,-10,-5,5,10,15], colors='lightgrey', linewidths=1)
    side_ax.contour(z[iz], yh[iy], np.min(OW[:,:,ix1-4:ix1+4], axis=2).transpose(), levels=[-0.01,-0.0025,-0.001], colors='k', linewidths=1, linestyles='-')
    # side_ax.quiver(z[iz][::qiz], yh[iy][::qix], -1*winterp[::qiz,::qix,ix1].transpose(), uinterp[::qiz,::qix,ix1].transpose(), scale=300, width=0.005, pivot='tail')
    # side_ax.set_xticks(np.arange(0,4))
    side_ax.invert_xaxis()
    
    plt.suptitle(f"$Z_{{sfc}}$ (shaded), OW$_{{sfc}}$ (black), $w_{{1km}}$ (gray)\n t={time:.0f} min", fontsize=12) 
    
    if figsave:
        plt.savefig(ip+f"imgs_overview/dbz_{time:03.0f}min_{img_str}_cross-sects.png", dpi=400)

#%% Load data for diagonal cross sections

from CM1utils import *


fp = '/Volumes/Promise_Pegasus_70TB/merger/merger-125m/'
ip = '/Users/morgan.schneider/Documents/merger/merger-125m/'

fn = 43
# 181 min -> 14 | 195 min -> 28 | 210 min -> 43 | 225 min -> 58 | 240 min -> 73


# Read output file
ds = nc.Dataset(fp+f"cm1out_{fn:06d}.nc")
time = ds.variables['time'][:].data[0]/60
xh = ds.variables['xh'][:].data
xf = ds.variables['xf'][:].data
yh = ds.variables['yh'][:].data
yf = ds.variables['yf'][:].data
z = ds.variables['z'][:].data
zf = ds.variables['zf'][:].data

xlims = [-50,10] #[-55,25]
ylims = [-90,-30] #[-100,-20]

ixw = np.where(xh >= xlims[0])[0][0] # west bound
ixe = np.where(xh >= xlims[1])[0][0] # east bound
iys = np.where(yh >= ylims[0])[0][0] # south bound
iyn = np.where(yh >= ylims[1])[0][0] # north bound
ix = slice(ixw,ixe+1)
iy = slice(iys,iyn+1)

iz1 = np.where(z >= 1)[0][0]
iz15 = np.where(z >= 1.5)[0][0]
iz2 = np.where(z >= 2)[0][0]
iz5 = np.where(z >= 5)[0][0]
iz200 = np.where(z >= 0.2)[0][0]
iz = slice(0,iz5+1)

prs = ds.variables['prs'][:].data[0,iz,iy,ix]
uinterp = ds.variables['uinterp'][:].data[0,iz,iy,ix]
vinterp = ds.variables['vinterp'][:].data[0,iz,iy,ix]
winterp = ds.variables['winterp'][:].data[0,iz,iy,ix]
# zvort = ds.variables['zvort'][:].data[0,iz,iy,ix]
# S_N = np.gradient(uinterp, xh[ix]*1000, axis=2) - np.gradient(vinterp, yh[iy]*1000, axis=1)
# S_S = np.gradient(vinterp, xh[ix]*1000, axis=2) + np.gradient(uinterp, yh[iy]*1000, axis=1)
# OW = S_N**2 + S_S**2 - zvort**2
# del S_N,S_S,zvort
# xvort = ds.variables['xvort'][:].data[0,iz,iy,ix]
# yvort = ds.variables['yvort'][:].data[0,iz,iy,ix]

if 'dbz2' in ds.variables:
    dbz = ds.variables['dbz2'][:].data[0,iz,iy,ix]
    thr = ds.variables['th'][:].data[0,iz,iy,ix] * (1 + 0.61*ds.variables['qv'][:].data[0,iz,iy,ix] - 
                (ds.variables['qc'][:].data[0,iz,iy,ix] + ds.variables['qr'][:].data[0,iz,iy,ix] + 
                  ds.variables['qi'][:].data[0,iz,iy,ix] + ds.variables['qs'][:].data[0,iz,iy,ix] + 
                  ds.variables['qg'][:].data[0,iz,iy,ix] + ds.variables['qhl'][:].data[0,iz,iy,ix]))
    thr0 = ds.variables['th0'][:].data[0,iz,iy,ix] * (1 + 0.61*ds.variables['qv0'][:].data[0,iz,iy,ix])
    thrpert = thr - thr0
    # B = 9.8 * (thrpert/thr0)
    del thr,thr0
    prspert = (prs - ds.variables['prs0'][:].data[0,iz,iy,ix])/100 # hPa
    ds.close()
else:
    dbz = ds.variables['dbz'][:].data[0,iz,iy,ix]
    ds.close()
    
    ds = nc.Dataset(fp+f"base/cm1out_000001.nc")
    prspert = (prs - ds.variables['prs0'][:].data[0,iz,iy,ix])/100 # hPa
    ds.close()
    
    ds = nc.Dataset(fp+f"pp/dyn_{fn:06d}.nc")
    thrpert = ds.variables['thrpert'][:].data[iz,iy,ix]
    # thr0 = np.moveaxis(np.tile(ds.variables['thr0'][:].data, (thrpert1.shape[1], thrpert1.shape[2], 1)), -1, 0)
    # B = 9.8 * (thrpert/thr0)
    # del thr0
    ds.close()

img_str = 's1'


wlims = [-15,15]
thr_lims = [-15,15]
pp_lims = [-20,20]
vort_lims = [-0.1,0.1]
dbz_lims = [0,70]


qix = int(np.round(len(xh[(xh>=xlims[0]) & (xh<=xlims[1])])/60)*2) # scale wind vector spacing to number of grid pts
qiz = 3

#%% Calculate and plot the diagonal cross sections

# from CM1utils import *

# Cross sections!
x1 = -25 # cross section start point
y1 = -81
x2 = 0 # cross section end point
y2 = -58


proj_angle = np.arctan2(x2-x1, y2-y1) # angle of cross section surface
U_proj,nu,nv = proj_winds(uinterp, vinterp, proj_angle) # horizontal wind projected onto cross section surface

u_cross,x_cross = vert_cross_section(U_proj, xh[ix], yh[iy], start=[x1,y1], end=[x2,y2])
w_cross,x_cross = vert_cross_section(winterp, xh[ix], yh[iy], start=[x1,y1], end=[x2,y2])
thr_cross,x_cross = vert_cross_section(thrpert, xh[ix], yh[iy], start=[x1,y1], end=[x2,y2])
dbz_cross,x_cross = vert_cross_section(dbz, xh[ix], yh[iy], start=[x1,y1], end=[x2,y2])

# zvort_cross,x_cross = vert_cross_section(zvort, xh[ix], yh[iy], start=[x1,y1], end=[x2,y2])
# OW_cross,x_cross = vert_cross_section(OW, xh[ix], yh[iy], start=[x1,y1], end=[x2,y2])
# pp_cross,x_cross = vert_cross_section(prspert, xh[ix], yh[iy], start=[x1,y1], end=[x2,y2])
# hvort_proj,nvx,nvy = proj_winds(xvort, yvort, proj_angle) # horizontal vorticity along cross section
# hvort_cross,x_cross = vert_cross_section(hvort_proj, xh[ix], yh[iy], start=[x1,y1], end=[x2,y2])
# B_cross,x_cross = vert_cross_section(B, xh[ix], yh[iy], start=[x1,y1], end=[x2,y2])





fig,ax = plt.subplots(1,1,figsize=(8,6))
plot_cfill(xh[ix], yh[iy], np.ma.masked_array(dbz[0,:,:], dbz[0,:,:]<1), 'dbz', ax, datalims=dbz_lims, xlims=xlims, ylims=ylims)
# ax.contour(xh[ix], yh[iy], winterp[iz1,:,:], levels=[-10,-5,5,10], colors='k', linewidths=1)
ax.quiver(xh[ix][::qix], yh[iy][::qix], uinterp[iz1,::qix,::qix], vinterp[iz1,::qix,::qix], scale=600, width=0.002, pivot='tail')
ax.plot([x1,x2], [y1,y2], '-k', linewidth=1.5)
plt.show()

fig,ax = plt.subplots(1,1,figsize=(8,6))
plot_cfill(xh[ix], yh[iy], thrpert[0,:,:], 'thrpert', ax, datalims=[-14,0], cmap='YlGnBu_r', xlims=xlims, ylims=ylims)
ax.contour(xh[ix], yh[iy], np.max(winterp[0:iz1+1,:,:], axis=0), levels=[5], colors='k', linewidths=1)
ax.quiver(xh[ix][::qix], yh[iy][::qix], uinterp[0,::qix,::qix], vinterp[0,::qix,::qix], scale=600, width=0.002, pivot='tail')
ax.plot([x1,x2], [y1,y2], '-k', linewidth=1.5)
plt.show()

# fig,ax = plt.subplots(1,1,figsize=(8,6))
# plot_cfill(xh[ix], yh[iy], winterp[iz1,:,:], 'w', ax, datalims=w_lims, xlims=xlims, ylims=ylims)
# # ax.contour(xh[ix], yh[iy], winterp[iz1,:,:], levels=[-10,-5,5,10], colors='k', linewidths=1)
# ax.quiver(xh[ix][::qix], yh[iy][::qix], uinterp[iz1,::qix,::qix], vinterp[iz1,::qix,::qix], scale=600, width=0.002, pivot='tail')
# ax.plot([x1,x2], [y1,y2], '-k', linewidth=2)
# plt.show()

fig,ax = plt.subplots(1,1,figsize=(8,6))
plot_cfill(xh[ix], yh[iy], np.ma.masked_array(dbz[0,:,:], dbz[0,:,:]<1), 'dbz', ax, datalims=dbz_lims, xlims=xlims, ylims=ylims)
# ax.contour(xh[ix], yh[iy], winterp[iz1,:,:], levels=[-10,-5,5,10], colors='k', linewidths=1)
ax.quiver(xh[ix][::qix], yh[iy][::qix], nu[iz1,::qix,::qix], nv[iz1,::qix,::qix], scale=600, width=0.002, pivot='tail')
ax.plot([x1,x2], [y1,y2], '-k', linewidth=1.5)
plt.show()


xl = [x1+5,x2-5]

qih = 6
qiz = 3

fig,ax = plt.subplots(1,1,figsize=(8,4))
plot_cfill(x_cross, z[iz], w_cross, 'w', ax, datalims=w_lims, xlims=xl, ylims=[0,4])
ax.quiver(x_cross[::qih], z[iz][::qiz], u_cross[::qiz,::qih], w_cross[::qiz,::qih], scale=600, width=0.002, pivot='tail')
plt.show()

fig,ax = plt.subplots(1,1,figsize=(8,4))
plot_cfill(x_cross, z[iz], thr_cross, 'thrpert', ax, datalims=thr_lims, xlims=xl, ylims=[0,4])
ax.quiver(x_cross[::qih], z[iz][::qiz], u_cross[::qiz,::qih], w_cross[::qiz,::qih], scale=600, width=0.002, pivot='tail')
plt.show()

fig,ax = plt.subplots(1,1,figsize=(8,4))
plot_cfill(x_cross, z[iz], w_cross, 'w', ax, datalims=w_lims, xlims=xl, ylims=[0,4])
ax.quiver(x_cross[::qih], z[iz][::qiz], 0*u_cross[::qiz,::qih], w_cross[::qiz,::qih], scale=250, width=0.0025, pivot='tail')
plt.show()

qih = 5
# qiz = 3
fig,ax = plt.subplots(1,1,figsize=(8,4))
plot_cfill(x_cross, z[iz], thr_cross, 'thrpert', ax, datalims=thr_lims, xlims=xl, ylims=[0,4])
ax.quiver(x_cross[::qih], z[iz][::qiz], 0*u_cross[::qiz,::qih], np.ma.masked_array(w_cross, w_cross<0)[::qiz,::qih], color='r', scale=250, width=0.0025, pivot='tail')
ax.quiver(x_cross[::qih], z[iz][::qiz], 0*u_cross[::qiz,::qih], np.ma.masked_array(w_cross, w_cross>-0)[::qiz,::qih], color='b', scale=250, width=0.0025, pivot='tail')
plt.show()


#%% Separate cross sections

fp = '/Volumes/Promise_Pegasus_70TB/merger/merger-125m/base/'
ip = '/Users/morgan.schneider/Documents/merger/merger-125m/'

fn = 13
# 181 min -> 14 | 195 min -> 28 | 210 min -> 43 | 225 min -> 58 | 240 min -> 73

# Read output file
ds = nc.Dataset(fp+f"cm1out_{fn:06d}.nc")
time = ds.variables['time'][:].data[0]/60
xh = ds.variables['xh'][:].data
xf = ds.variables['xf'][:].data
yh = ds.variables['yh'][:].data
yf = ds.variables['yf'][:].data
z = ds.variables['z'][:].data
zf = ds.variables['zf'][:].data

xlims = [-70,10]
# ylims = [-140,-60]
iy1 = np.where(yh >= -105)[0][0]
iy = slice(iy1-40,iy1+41)
ylims = [yh[iy1-40], yh[iy1+40]]


ixw = np.where(xh >= xlims[0])[0][0]
ixe = np.where(xh >= xlims[1])[0][0]
ix = slice(ixw,ixe+1)
# iys = np.where(yh >= ylims[0])[0][0]
# iyn = np.where(yh >= ylims[1])[0][0]
# iy = slice(iys,iyn+1)

iz = np.where(z >= 5)[0][0]


if 'dbz2' in ds.variables:
    dbz = ds.variables['dbz2'][:].data[0,slice(0,iz+1),iy,ix]
    thr = ds.variables['th'][:].data[0,slice(0,iz+1),iy,ix] * (1 + 0.61*ds.variables['qv'][:].data[0,slice(0,iz+1),iy,ix] - 
                (ds.variables['qc'][:].data[0,slice(0,iz+1),iy,ix] + ds.variables['qr'][:].data[0,slice(0,iz+1),iy,ix] + 
                 ds.variables['qi'][:].data[0,slice(0,iz+1),iy,ix] + ds.variables['qs'][:].data[0,slice(0,iz+1),iy,ix] + 
                 ds.variables['qg'][:].data[0,slice(0,iz+1),iy,ix] + ds.variables['qhl'][:].data[0,slice(0,iz+1),iy,ix]))
    thr0 = ds.variables['th0'][:].data[0,slice(0,iz+1),iy,ix] * (1 + 0.61*ds.variables['qv0'][:].data[0,slice(0,iz+1),iy,ix])
    thrpert = thr - thr0
    del thr,thr0
else:
    dbz = ds.variables['dbz'][:].data[0,slice(0,iz+1),iy,ix]
# prspert = ds.variables['prspert'][:].data[0,slice(0,iz+1),iy,ix]/100
uinterp = ds.variables['uinterp'][:].data[0,slice(0,iz+1),iy,ix]
vinterp = ds.variables['vinterp'][:].data[0,slice(0,iz+1),iy,ix]
winterp = ds.variables['winterp'][:].data[0,slice(0,iz+1),iy,ix]
zvort = ds.variables['zvort'][:].data[0,slice(0,iz+1),iy,ix]
S_N = np.gradient(uinterp, xh*1000, axis=2) - np.gradient(vinterp, yh*1000, axis=1)
S_S = np.gradient(vinterp, xh*1000, axis-2) + np.gradient(uinterp, yh*1000, axis=1)
OW = S_N**2 + S_S**2 - zvort**2
del S_N,S_S,zvort,vinterp
ds.close()

# Read dyn
ds = nc.Dataset(fp+f"pp/dyn_{fn:06d}.nc")
thrpert = ds.variables['thrpert'][:].data[slice(0,iz+1),iy,ix]
ds.close()

img_str = 's1'


wlims = [-15,15]
thr_lims = [-15,15]
pp_lims = [-20,20]
vort_lims = [-0.1,0.1]
dbz_lims = [0,70]


qix = int(np.round(len(xh[(xh>=xlims[0]) & (xh<=xlims[1])])/60)*2) # scale wind vector spacing to number of grid pts
qiz = 3


figsave = False


# Separate x-z and y-z cross sections
# LOOK INTO WRF FUNCTION wrf.interp2dxy and wrf.vertcross for diagonal cross sections!
if False:
    
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    plot_cfill(xh, z, dbz[:,iy1,:], 'dbz', ax, datalims=dbz_lims, xlims=xlims, ylims=[0,5])
    # ax.contour(xh, z, np.mean(thrpert[:,iy1-4:iy1+4,:], axis=1), levels=[-5], colors='b', linewidths=1, linestyles='-')
    ax.contour(xh, z, winterp[:,iy1,:], levels=[-10,-5,5,10], colors='w', linewidths=1)
    # ax.contour(xh, z, np.mean(zvort[:,iy1-4:iy1+4,:], axis=1), levels=[0.025,0.05,0.075,0.1], colors='k', linewidths=1)
    ax.contour(xh, z, np.min(OW[:,iy1-4:iy1+4,:], axis=1), levels=[-0.01,-0.0025,-0.001], colors='k', linewidths=1, linestyles='-')
    ax.quiver(xh[::qix], z[::qiz], uinterp[::qiz,iy1,::qix], winterp[::qiz,iy1,::qix], scale=600, width=0.002, pivot='tail')
    ax.set_xlabel('x (km)')
    ax.set_ylabel('z (km)')
    ax.set_title(f"$Z_H$ (shaded), OW (black), $w$ (white), \u03B8\u1D68'=-8 K (blue)\n t={time:.0f} min", fontsize=14)
    
    
    # fig,ax = plt.subplots(1,1,figsize=(8,6))
    # plot_cfill(yf, zf, np.ma.masked_array(dbz[:,:,ix1], dbz[:,:,ix1]<1).transpose(), 'dbz', ax, datalims=dbz_lims, xlims=ylims, ylims=[0,5])
    # ax.contour(yh, z, np.mean(thrpert[:,:,ix1-4:ix1+4], axis=2), levels=[-5], colors='b', linewidths=1, linestyles='-')
    # ax.contour(yh, z, np.mean(winterp[:,:,ix1-4:ix1+4], axis=2), levels=[-15,-10,-5,5,10,15], colors='w', linewidths=1)
    # # ax.contour(yh, z, np.mean(zvort[:,:,ix1-4:ix1+4], axis=2), levels=[0.025,0.05,0.075,0.1], colors='k', linewidths=1)
    # ax.contour(z, yh, np.mean(OW[:,:,ix1-4:ix1+4], axis=2), levels=[-0.01,-0.0025,-0.001], colors='k', linewidths=1, linestyles='-')
    # ax.quiver(yh[::qix], z[::qiz], uinterp[::qiz,::qix,ix1], winterp[::qiz,::qix,ix1], scale=600, width=0.002, pivot='tail')
    # ax.set_xlabel('y (km)')
    # ax.set_ylabel('z (km)')
    # ax.set_title(f"$Z_H$ (shaded), OW (black), $w$ (white), \u03B8\u1D68'=-8 K (blue)\n t={time:.0f} min", fontsize=14)

if False:
    fig,ax = plt.subplots(1,1,figsize=(10,4))
    plot_cfill(xf, zf, winterp[:,iy1,:], 'w', ax, datalims=w_lims, xlims=xlims, ylims=[0,5])
    ax.contour(xh, z, thrpert[:,iy1,:], levels=[-5], colors='g', linewidths=1, linestyles='-')
    # ax.contour(xh, z, winterp[:,iy1,:], levels=[-10,-5,5,10], colors='w', linewidths=1)
    # ax.contour(xh, z, np.mean(zvort[:,iy1-4:iy1+4,:], axis=1), levels=[0.025,0.05,0.075,0.1], colors='k', linewidths=1)
    ax.contour(xh, z, np.min(OW[:,iy1-4:iy1+4,:], axis=1), levels=[-0.01,-0.0025,-0.001], colors='hotpink', linewidths=1, linestyles='-')
    ax.quiver(xh[::qix], z[::qiz], uinterp[::qiz,iy1,::qix], winterp[::qiz,iy1,::qix], scale=600, width=0.002, pivot='tail')
    ax.set_xlabel('x (km)')
    ax.set_ylabel('z (km)')
    # ax.set_title(f"$Z_H$ (shaded), OW (black), $w$ (white), \u03B8\u1D68'=-8 K (blue)\n t={time:.0f} min", fontsize=14)
    
    fig,ax = plt.subplots(1,1,figsize=(10,4))
    plot_cfill(xf, zf, thrpert[:,iy1,:], 'thrpert', ax, datalims=[-10,10], cmap='YlGnBu_r', xlims=xlims, ylims=[0,5])
    ax.contour(xh, z, winterp[:,iy1,:], levels=[-10,-5,5,10], colors='w', linewidths=1)
    ax.contour(xh, z, prspert[:,iy1,:], levels=[1,2,3,4,5], colors='r', linewidths=1)
    ax.contour(xh, z, prspert[:,iy1,:], levels=[-5,-4,-3,-2,-1], colors='k', linewidths=1, linestyles='-')
    # ax.contour(xh, z, np.mean(zvort[:,iy1-4:iy1+4,:], axis=1), levels=[0.025,0.05,0.075,0.1], colors='k', linewidths=1)
    ax.contour(xh, z, np.min(OW[:,iy1-4:iy1+4,:], axis=1), levels=[-0.01,-0.0025,-0.001], colors='hotpink', linewidths=1, linestyles='-')
    ax.quiver(xh[::qix], z[::qiz], uinterp[::qiz,iy1,::qix], winterp[::qiz,iy1,::qix], scale=600, width=0.002, pivot='tail')
    ax.set_xlabel('x (km)')
    ax.set_ylabel('z (km)')
    # ax.set_title(f"$Z_H$ (shaded), OW (black), $w$ (white), \u03B8\u1D68'=-8 K (blue)\n t={time:.0f} min", fontsize=14)
    

    


#%% MVs and TLVs lol good luck


# Lovell and Parker 2022 criteria --> at 500, 1000, and 3000 m
ow0_lp22 = -0.0025 # Null TLV - sfc OW
ow1_lp22 = -0.01 # Weak TLV - sfc OW
ow2_lp22 = -0.0225 # Strong TLV - sfc OW
ow3_lp22 = -0.04 # Intense TLV - sfc OW

# Sherburn and Parker 2019 criteria --> must be over 1.5 km^2 area for >=5 min
ow_mv_sp19 = -0.001 # MV vortex - 1.5 km OW
w_mv_sp19 = 5 # MV updraft - 1.5 km w


fp = '/Volumes/Promise_Pegasus_70TB/merger/merger-125m/'
ip = '/Users/morgan.schneider/Documents/merger/merger-125m/'

from skimage.feature import peak_local_max
from skimage.color import rgb2gray
from PIL import Image


ds = nc.Dataset(fp+'cm1out_000014.nc')
xh = ds.variables['xh'][:].data
xf = ds.variables['xf'][:].data
yh = ds.variables['yh'][:].data
yf = ds.variables['yf'][:].data
z = ds.variables['z'][:].data
zf = ds.variables['zf'][:].data
ds.close()

ix1 = np.where(xh >= -60)[0][0]
ix2 = np.where(xh >= 30)[0][0]
iy1 = np.where(yh >= -110)[0][0]
iy2 = np.where(yh >= -20)[0][0]
iz15 = np.where(z >= 1.5)[0][0]
ix = slice(ix1,ix2)
iy = slice(iy1,iy2)


# attempted sherburn and parker 2019 mv object identification
# i hate my life
for fn in np.arange(14,74):
    print(f"cm1out_{fn:06d}")
    ds = nc.Dataset(fp+f"cm1out_{fn:06d}.nc")
    time = ds.variables['time'][:].data[0]/60
    w = ds.variables['winterp'][:].data[0,0:iz1+1,iy,ix]
    zvort = ds.variables['zvort'][:].data[0,0:iz1+1,iy,ix]
    u = ds.variables['uinterp'][:].data[0,0:iz1+1,iy,ix]
    v = ds.variables['vinterp'][:].data[0,0:iz1+1,iy,ix]
    S_N = np.gradient(u, xh[ix]*1000, axis=2) - np.gradient(v, yh[iy]*1000, axis=1)
    S_S = np.gradient(v, xh[ix]*1000, axis=2) + np.gradient(u, yh[iy]*1000, axis=1)
    OW = S_N**2 + S_S**2 - zvort**2
    del u,v,S_N,S_S,zvort
    ds.close()
    
    w_max = np.max(w, axis=0)
    OW_min = np.min(OW, axis=0)
    
    # Lovell and Parker 2022 TLV criteria
    isubtlv = (OW[0,:,:] <= -0.0025) & (OW[0,:,:] > -0.01)
    itlv1 = (OW[0,:,:] <= -0.01) & (OW[0,:,:] > -0.0225)
    itlv2 = (OW[0,:,:] <= -0.0225) & (OW[0,:,:] > -0.04)
    itlv3 = (OW[0,:,:] <= -0.04)
    itlv = (OW[0,:,:] <= -0.01)
    # Sherburn and Parker 2019 mesovortex criteria
    imv = (OW_min <= -0.001) & (w_max >= 5)
    
    mv_flg = np.where(imv, 1, 0)
    im = plt.imshow(mv_flg, cmap='gray', origin='lower', vmin=0, vmax=1)
    img = Image.fromarray(np.uint8(im.get_cmap()(im.get_array())*255)).convert('RGB')
    mv_img = rgb2gray(img)
    xy_mv = peak_local_max(mv_img, min_distance=10, threshold_abs=0.05)
    
    xc_mv = np.zeros(shape=(len(xy_mv),), dtype=float)
    yc_mv = np.zeros(shape=(len(xy_mv),), dtype=float)
    for i in range(len(xy_mv)):
        xi = np.where(xh[ix] == xy_mv[i,1])[0][0]
        yi = np.where(yh[iy] == xy_mv[i,0])[0][0]
        idx = slice(xi-5,xi+6)
        jdy = slice(yi-5,yi+6)
        w_obj = w_max[jdy,idx]
        OW_obj = OW_min[jdy,idx]
        if (np.percentile(w_obj,90) >= 5) and (np.percentile(OW_obj,10) <= -0.001):
            xc_mv[i] = xy_mv[i,1]
            yc_mv[i] = xy_mv[i,0]
    
    # make an array of the mv indices?
    
    tlv_flg = np.where(itlv, 1, 0)
    im = plt.imshow(tlv_flg, cmap='gray', origin='lower', vmin=0, vmax=1)
    img = Image.fromarray(np.uint8(im.get_cmap()(im.get_array())*255)).convert('RGB')
    tlv_img = rgb2gray(img)
    xy_tlv = peak_local_max(tlv_img, min_distance=10, threshold_abs=0.05)
    
    xc_tlv = np.zeros(shape=(len(xy_tlv),), dtype=float)
    yc_tlv = np.zeros(shape=(len(xy_tlv),), dtype=float)
    for i in range(len(xy_tlv)):
        xi = np.where(xh[ix] == xy_tlv[i,1])[0][0]
        yi = np.where(yh[iy] == xy_tlv[i,0])[0][0]
        idx = slice(xi-5,xi+6)
        jdy = slice(yi-5,yi+6)
        OW_obj = OW[0,jdy,idx]
        w_obj = w_max[jdy,idx] # colocated with an updraft?
        if (np.percentile(OW_obj,10) <= -0.01) and (np.percentile(w_obj,90) >= 5):
            xc_tlv[i] = xy_tlv[i,1]
            yc_tlv[i] = xy_tlv[i,0]
    
    # make an array of the tlv indices?
    
    subtlv_flag = np.where(isubtlv, 1, 0)
    im = plt.imshow(subtlv_flg, cmap='gray', origin='lower', vmin=0, vmax=1)
    img = Image.fromarray(np.uint8(im.get_cmap()(im.get_array())*255)).convert('RGB')
    subtlv_img = rgb2gray(img)
    xy_subtlv = peak_local_max(subtlv_img, min_distance=10, threshold_abs=0.05)
    
    xc_subtlv = np.zeros(shape=(len(xy_subtlv),), dtype=float)
    yc_subtlv = np.zeros(shape=(len(xy_subtlv),), dtype=float)
    for i in range(len(xy_subtlv)):
        xi = np.where(xh[ix] == xy_subtlv[i,1])[0][0]
        yi = np.where(yh[iy] == xy_subtlv[i,0])[0][0]
        idx = slice(xi-5,xi+6)
        jdy = slice(yi-5,yi+6)
        OW_obj = OW[0,jdy,idx]
        w_obj = w_max[jdy,idx] # colocated with an updraft?
        if (np.percentile(OW_obj,10) <= -0.01) and (np.percentile(w_obj,90) >= 5):
            xc_subtlv[i] = xy_subtlv[i,1]
            yc_subtlv[i] = xy_subtlv[i,0]
    
    # make an array of the sub tlv indices?




#%% Base state profile and CI time series

from metpy.plots import SkewT,Hodograph
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

figsave = False

if False:
    if 'T0' not in locals():
        ds = nc.Dataset('/Volumes/Promise_Pegasus_70TB/merger/merger-125m/base/cm1out_000013.nc')
        th0 = ds.variables['th0'][:].data[0,:,-1,-1]
        qv0 = ds.variables['qv0'][:].data[0,:,-1,-1]
        prs0 = ds.variables['prs0'][:].data[0,:,-1,-1]
        u0 = ds.variables['u0'][:].data[0,:,-1,-1]
        v0 = ds.variables['v0'][:].data[0,:,-1,-1]
        ds.close()
        
        T0 = th0 * (prs0/100000)**0.286
        e0 = (qv0*prs0/100) / (0.622+qv0)
        Td0 = 243.5 / ((17.67 / (np.log(e0/6.112))) - 1) + 273.15
        T0_parcel = mc.parcel_profile((prs0/100)*units.hPa, T0[0]*units.K, Td0[0]*units.K)
    
    
    fig = plt.figure(figsize=(12,10))
    skew = SkewT(fig=fig)
    skew.plot(prs0/100, (T0-273.15), '-r', linewidth=2)
    skew.plot(prs0/100, (Td0-273.15), '-g', linewidth=2)
    skew.plot(prs0/100, np.array(T0_parcel.magnitude[:])-273.15, '-k', linewidth=2)
    skew.plot_dry_adiabats(linewidths=1)
    skew.plot_moist_adiabats(linewidths=1)
    # skew.plot_mixing_lines()
    skew.shade_cape(prs0/100, T0-273.15, T0_parcel.magnitude[:]-273.15)
    skew.ax.set_ylim(1000, 100)
    skew.ax.set_xlim(-40, 30)
    ax_hod = inset_axes(skew.ax, '42%', '42%', loc=1)
    H = Hodograph(ax_hod, component_range=40.)
    H.add_grid(increment=10)
    H.plot(u0, v0, color='k', linewidth=1.5)
    # skew.ax.set_title('Base state environment (1959 UTC 30 March 2022 PERiLS IOP2)')
    skew.ax.set_xlabel('Temperature (C)', fontsize=14)
    skew.ax.set_ylabel('Pressure (hPa)', fontsize=14)
    if figsave:
        plt.savefig('/Users/morgan.schneider/Documents/merger/basestate.png', dpi=400)

if False:
    tmp1 = np.linspace(0,3600,3601)
    tmp2 = np.zeros(shape=(3601,), dtype=float)
    tmp3 = np.zeros(shape=(3601,), dtype=float)
    S0 = 0.021
    
    tmp2[1:1800] = S0/1800 * tmp1[1:1800]
    tmp2[1800] = S0
    tmp2[1801] = S0
    tmp2[1802:3600] = -S0/1799 * tmp1[1802:3600] + 3600*S0/1799
    
    
    fig,ax = plt.subplots(1,1,figsize=(10,4))
    ax.plot(tmp1, tmp2, 'k', linewidth=1.5)
    ax.axhline(y=S0, color='k', linestyle='--', linewidth=1)
    ax.set_xlabel('Time (s)', fontsize=14)
    ax.set_ylabel("S (K s$^{-1}$)", fontsize=14)
    # ax.set_title('Heating tendency', fontsize=14)
    ax.set_xlim([0,3600])
    ax.set_xticks(np.linspace(0,3600,7))
    ax.set_ylim([0,0.025])
    ax.text(600, 0.0215, "$S_0 = 0.021$", fontsize=14)
    if figsave:
        plt.savefig('/Users/morgan.schneider/Documents/merger/CI_tendency.png', dpi=400)
    
    
    for i in range(len(tmp1)):
        if i <= 1800:
            tmp3[i] = 0.5*S0/1800 * tmp1[i]**2
        elif i == 1801:
            tmp3[i] = 901*S0
        elif i > 1801:
            tmp3[i] = 901*S0 - 0.5*S0/1799*(tmp1[i]**2 - 1801**2) + 3600*S0/1799*(tmp1[i]-1801)
    
    fig,ax = plt.subplots(1,1,figsize=(10,4))
    ax.plot(tmp1, tmp3, 'k', linewidth=1.5)
    ax.set_xlabel('Time (s)', fontsize=14)
    ax.set_ylabel("\u03B8' (K)", fontsize=14)
    # ax.set_title("\u03B8'", fontsize=14)
    ax.set_xlim([0,3600])
    ax.set_ylim([0,40])
    ax.set_xticks(np.linspace(0,3600,7))
    if figsave:
        plt.savefig('/Users/morgan.schneider/Documents/merger/CI_theta.png', dpi=400)
    
    
    
    fig,((ax1),(ax2)) = plt.subplots(2,1,figsize=(10,8))
    
    ax1.plot(tmp1, tmp2, 'k', linewidth=1.5)
    ax1.axhline(y=S0, color='k', linestyle='--', linewidth=1)
    # ax1.set_xlabel('Time (s)')
    ax1.set_ylabel("S (K s$^{-1}$)", fontsize=14)
    # ax1.set_title('Heating tendency magnitude', fontsize=14)
    ax1.set_xlim([0,3600])
    ax1.set_xticks(np.linspace(0,3600,7))
    ax1.set_ylim([0,0.025])
    ax1.text(600, 0.0215, "$S_0 = 0.021$", fontsize=14)
    
    ax2.plot(tmp1, tmp3, 'k', linewidth=1.5)
    ax2.set_xlabel('Time (s)', fontsize=14)
    ax2.set_ylabel("\u03B8' (K)", fontsize=14)
    # ax2.set_title("\u03B8' magnitude", fontsize=14)
    ax2.set_xlim([0,3600])
    ax2.set_ylim([0,40])
    ax2.set_xticks(np.linspace(0,3600,7))
    
    if figsave:
        plt.savefig('/Users/morgan.schneider/Documents/merger/CI.png', dpi=400)
    

#%% Pickle files! for zoomed boxes

update_box = False
new_box = False



x1_pp = []
x2_pp = [x1_pp[i] + 10 for i in range(len(x1_pp))]

y1_pp = []
y2_pp = [y1_pp[i] + 10 for i in range(len(y1_pp))]


if update_box:
    dbfile = open('/Users/morgan.schneider/Documents/merger/supercell-125m/boxes_s2.pkl', 'rb')
    box = pickle.load(dbfile)
    dbfile.close()
    
    new_vars = {'x1_pp':x1_pp, 'x2_pp':x2_pp, 'y1_pp':y1_pp, 'y2_pp':y2_pp}
    box.update(new_vars)
    dbfile = open('/Users/morgan.schneider/Documents/merger/supercell-125m/boxes_s2.pkl', 'wb')
    pickle.dump(box, dbfile)
    dbfile.close()

if new_box:
    dbfile = open('/Users/morgan.schneider/Documents/merger/supercell-125m/boxes_s2.pkl', 'wb')
    box = {'x1_pp':x1_pp, 'x2_pp':x2_pp, 'y1_pp':y1_pp, 'y2_pp':y2_pp}
    pickle.dump(box, dbfile)
    dbfile.close()

    

#%% Random junk that I should probably keep around

# Reprocess reflectivity if needed
if False:
    import os
    
    fp = '/Volumes/Promise_Pegasus_70TB/merger/merger-125m/'
    files = glob(fp+'cm1out_0000*.nc')
    for n in range(len(files)):
        ds = nc.Dataset(files[n], 'r+')
        print(f"File {files[n]} ...")
        
        if 'dbz2' not in ds.variables:
            print(f"... Recalculating dBZ ...")
            dbz_new = calc_dbz(ds)
            
            print("... Writing dBZ to file ...")
            dbz2 = ds.createVariable('dbz2', 'f4', ('time', 'nk', 'nj', 'ni'))
            dbz2.units = 'dBZ'
            dbz2.long_name = 'recalculated reflectivity'
            dbz2[:,:,:,:] = dbz_new[:,:,:,:]
            
            del dbz2,dbz_new
            ds.close()
            print("... Done!")
            print("---")


# Check dyn pp values
if False:
    files = glob('/Volumes/Promise_Pegasus_70TB/merger/supercell-125m/pp/*')
    for ff in files:
        ds = nc.Dataset(ff)
        d = np.min(ds.variables['pgfz_b'][:].data)
        d2 = np.max(ds.variables['pgfz_b'][:].data[0:-2,:,:])
        print(f"{ff[-5:-3]}: {d:.3f}, {d2:.3f}")
        ds.close()




