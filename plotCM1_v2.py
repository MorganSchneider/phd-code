# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 17:00:30 2023

@author: morgan.schneider

Plotting CM1 outputs
"""

####################
### Load modules ###
####################

from CM1utils import *


#%% Load data

figsave = True

fp = 'D:/cm1out/merger/supercell3-500m/'
#fp = '/Users/morgan.schneider/Documents/cm1/merger/merger1/'
#stats = read_cm1out(fp+'cm1out_stats.nc')
fn = fp+'cm1out_000049.nc'

dsvars = ['time','xh','yh','z','xf','yf','zf',
          'winterp','thpert','dbz','zvort',
          'qvpert','qc',]
#          'prspert','th','xvort','uinterp','vinterp']

data = read_cm1out(fn, dsvars)

thrpert = data.th * (1 + 0.61*data.qv - data.qc)
data.SetAttr('thrpert', thrpert)

#%%

ulims = [-40., 40.]
wlims = [-15., 15.]
thlims = [294., 304.]
pthlims = [-10., 10.]
pplims = [-5., 5.]
hvlims = [-0.1, 0.1]
svlims = [-0.02, 0.02]
dbzlims = [0., 80.]
qvlims = [0., 4.]
nclims = [0., 20000.]


#%%
lev = 3
iz = np.where(data.z >= lev)[0][0]
iy = np.where(data.yh >= -70)[0][0]
xlims = [-180, 180] #[-100, 0]
ylims = [-180, 180] #[-140, -40]
qix = 35 #10


fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(12,8))

c1 = plot_cfill(data.xf, data.yf, data.dbz[0,iz,:,:], 'dbz', ax1, datalims=dbzlims)
#ax1.quiver(data.xh[::qix], data.yh[::qix], data.uinterp[0,0,::qix,::qix], data.vinterp[0,0,::qix,::qix], scale=275, width=0.004, pivot='middle')
ax1.set_ylabel('y distance (km)')
ax1.set_title(f"$Z_H$ at {1000*data.z[iz]:.0f} m", fontsize=14)
#rect = patches.Rectangle((30,20), 20, 20, linewidth=2, edgecolor='k', facecolor='none')
#ax1.add_patch(rect)
ax1.set_xlim(xlims[0],xlims[1])
ax1.set_ylim(ylims[0],ylims[1])

#c2 = plot_cfill(data.xf, data.yf, np.mean(data.winterp[0,:,:,:],0), 'w', ax2, datalims=wlims)
c2 = plot_cfill(data.xf, data.yf, data.winterp[0,iz,:,:], 'w', ax2, datalims=wlims)
#ax2.quiver(data.xh[::qix], data.yh[::qix], data.uinterp[0,0,::qix,::qix], data.vinterp[0,0,::qix,::qix], scale=275, width=0.004, pivot='middle')
ax2.set_title(f"$w$ at {1000*data.z[iz]:.0f} m", fontsize=14)
ax2.set_xlim(xlims[0],xlims[1])
ax2.set_ylim(ylims[0],ylims[1])

c3 = plot_cfill(data.xf, data.yf, data.thpert[0,iz,:,:], 'thpert', ax3, datalims=pthlims)
#ax3.quiver(data.xh[::qix], data.yh[::qix], data.uinterp[0,0,::qix,::qix], data.vinterp[0,0,::qix,::qix], scale=275, width=0.004, pivot='middle')
ax3.set_title(f"\u03B8' at {1000*data.z[iz]:.0f} m", fontsize=14)
#c3 = plot_cfill(data.xf, data.yf, data.prspert[0,iz,:,:]/100, 'prspert', ax3, datalims=pplims)
#ax3.set_title(f"p' at {1000*data.z[iz]:.0f} m (hPa)", fontsize=14)
ax3.set_xlabel('x distance (km)')
ax3.set_ylabel('y distance (km)')
ax3.set_xlim(xlims[0],xlims[1])
ax3.set_ylim(ylims[0],ylims[1])

c4 = plot_cfill(data.xf, data.yf, data.zvort[0,iz,:,:], 'zvort', ax4, datalims=svlims)
#ax4.quiver(data.xh[::qix], data.yh[::qix], data.uinterp[0,0,::qix,::qix], data.vinterp[0,0,::qix,::qix], scale=275, width=0.004, pivot='middle')
ax4.set_xlabel('x distance (km)')
ax4.set_title(f"\u03B6 at {1000*data.z[iz]:.0f} m", fontsize=14)
ax4.set_xlim(xlims[0],xlims[1])
ax4.set_ylim(ylims[0],ylims[1])

plt.suptitle(f"$t={data.time[0]/60:.0f}$ min", fontsize=14)
plt.show()
if figsave:
    plt.savefig(fp+f"dataH_{data.time[0]/60:.0f}min_{round(1000*data.z[iz],-2):.0f}m.png")

zlims = [0, 15] #[0, 6]
xlims = [-180, 180] #[-80, -20]
qix = 30 #6
qiz = 15 #8

# fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(12,8))

# c1 = plot_cfill(data.xf, data.zf, np.mean(data.dbz[0,:,:,:],1), 'dbz', ax1, datalims=dbzlims)
# #c1 = plot_cfill(data.xf, data.zf, data.dbz[0,:,iy,:], 'dbz', ax1, datalims=dbzlims)
# #ax1.quiver(data.xh[::qix], data.z[::qiz], data.uinterp[0,::qiz,iy,::qix], data.winterp[0,::qiz,iy,::qix], scale=250, width=0.004, pivot='middle')
# ax1.set_ylabel('Height (km)')
# ax1.set_title('$Z_H$', fontsize=14)
# ax1.set_ylim(zlims[0],zlims[1])
# ax1.set_xlim(xlims[0],xlims[1])

# c2 = plot_cfill(data.xf, data.zf, np.mean(data.winterp[0,:,:,:],1), 'w', ax2, datalims=wlims)
# #c2 = plot_cfill(data.xf, data.zf, data.winterp[0,:,iy,:], 'w', ax2, datalims=wlims)
# #ax2.quiver(data.xh[::qix], data.z[::qiz], data.uinterp[0,::qiz,iy,::qix], data.winterp[0,::qiz,iy,::qix], scale=250, width=0.004, pivot='middle')
# ax2.set_title("$w$", fontsize=14)
# ax2.set_ylim(zlims[0],zlims[1])
# ax2.set_xlim(xlims[0],xlims[1])

# c3 = plot_cfill(data.xf, data.zf, np.mean(data.thpert[0,:,:,:],1), 'thpert', ax3, datalims=pthlims)
# #c3 = plot_cfill(data.xf, data.zf, data.thpert[0,:,iy,:], 'thpert', ax3, datalims=pthlims)
# #ax3.quiver(data.xh[::qix], data.z[::qiz], data.uinterp[0,::qiz,iy,::qix], data.winterp[0,::qiz,iy,::qix], scale=250, width=0.004, pivot='middle')
# ax3.set_title("\u03B8'", fontsize=14)
# #c3 = plot_cfill(data.xf, data.zf, data.prspert[0,:,iy,:]/100, 'prspert', ax3, datalims=pplims)
# #ax3.set_title(f"p' (hPa)", fontsize=14)
# # c3 = plot_cfill(data.xf, data.zf, np.mean(data.uinterp[0,:,210:260,:],1), 'u', ax3, datalims=ulims)
# # ax3.set_title("$u$", fontsize=14)
# ax3.set_xlabel('x distance (km)')
# ax3.set_ylabel('Height (km)')
# ax3.set_ylim(zlims[0],zlims[1])
# ax3.set_xlim(xlims[0],xlims[1])

# c4 = plot_cfill(data.xf, data.zf, np.mean(data.zvort[0,:,:,:],1), 'zvort', ax4, datalims=svlims)
# #c4 = plot_cfill(data.xf, data.zf, data.zvort[0,:,iy,:], 'zvort', ax4, datalims=svlims)
# #ax4.quiver(data.xh[::qix], data.z[::qiz], data.uinterp[0,::qiz,iy,::qix], data.winterp[0,::qiz,iy,::qix], scale=250, width=0.004, pivot='middle')
# ax4.set_xlabel('x distance (km)')
# ax4.set_title("\u03B6", fontsize=14)
# ax4.set_ylim(zlims[0],zlims[1])
# ax4.set_xlim(xlims[0],xlims[1])

# plt.suptitle(f"$t={data.time[0]/60:.0f}$ min", fontsize=14)
# plt.show()
# if figsave:
#     plt.savefig(fp+f"dataV_{data.time[0]:.0f}s.png")





#%%

iz = np.where(data.z >= 0)[0][0]
iy = np.where(data.yh >= -52)[0][0]
figsave = True

fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(12,8))
xl = [-60, 0] #[-70, -10]
yl = [-60, 0] #[-160, -100]

c1 = plot_cfill(data.xf, data.yf, data.dbz[0,iz,:,:], 'dbz', ax1, datalims=dbzlims)
ax1.plot([xl[0],xl[1]], [data.yh[iy], data.yh[iy]], '--k', linewidth=1.5)
# ax1.plot([xl[0],xl[1]], [data.yh[iy2], data.yh[iy2]], '--k', linewidth=1.5)
# ax1.plot([xl[0],xl[1]], [data.yh[iy3], data.yh[iy3]], '--k', linewidth=1.5)
ax1.set_ylabel('y distance (km)')
ax1.set_title(f"$Z_H$ at {1000*data.z[iz]:.0f} m", fontsize=14)
ax1.set_xlim(xl[0], xl[1])
ax1.set_ylim(yl[0], yl[1])

c2 = plot_cfill(data.xf, data.yf, data.winterp[0,iz,:,:], 'w', ax2, datalims=wlims)
ax2.plot([xl[0],xl[1]], [data.yh[iy], data.yh[iy]], '--k', linewidth=1.5)
# ax2.plot([xl[0],xl[1]], [data.yh[iy2], data.yh[iy2]], '--k', linewidth=1.5)
# ax2.plot([xl[0],xl[1]], [data.yh[iy3], data.yh[iy3]], '--k', linewidth=1.5)
ax2.set_title(f"$w$ at {1000*data.z[iz]:.0f} m", fontsize=14)
ax2.set_xlim(xl[0], xl[1])
ax2.set_ylim(yl[0], yl[1])

c3 = plot_cfill(data.xf, data.yf, data.thpert[0,iz,:,:], 'thpert', ax3, datalims=pthlims)
#c3 = plot_cfill(data.xf, data.yf, data.prspert[0,iz,:,:]/100, 'prspert', ax3, datalims=pplims)
ax3.plot([xl[0],xl[1]], [data.yh[iy], data.yh[iy]], '--k', linewidth=1.5)
# ax3.plot([xl[0],xl[1]], [data.yh[iy2], data.yh[iy2]], '--k', linewidth=1.5)
# ax3.plot([xl[0],xl[1]], [data.yh[iy3], data.yh[iy3]], '--k', linewidth=1.5)
ax3.set_xlabel('x distance (km)')
ax3.set_ylabel('y distance (km)')
ax3.set_title(f"\u03B8' at {1000*data.z[iz]:.0f} m", fontsize=14)
#ax3.set_title(f"p' at {1000*data.z[0]:.0f} m", fontsize=14)
ax3.set_xlim(xl[0],xl[1])
ax3.set_ylim(yl[0],yl[1])

c4 = plot_cfill(data.xf, data.yf, data.zvort[0,iz,:,:], 'zvort', ax4, datalims=svlims)
ax4.plot([xl[0],xl[1]], [data.yh[iy], data.yh[iy]], '--k', linewidth=1.5)
# ax4.plot([xl[0],xl[1]], [data.yh[iy2], data.yh[iy2]], '--k', linewidth=1.5)
# ax4.plot([xl[0],xl[1]], [data.yh[iy3], data.yh[iy3]], '--k', linewidth=1.5)
ax4.set_xlabel('x distance (km)')
ax4.set_title(f"\u03B6 at {1000*data.z[iz]:.0f} m", fontsize=14)
ax4.set_xlim(xl[0],xl[1])
ax4.set_ylim(yl[0],yl[1])

plt.suptitle(f"$t={data.time[0]/60:.0f}$ min", fontsize=14)
plt.show()
if figsave:
    plt.savefig(fp+f"dataH_zoom_{data.time[0]/60:.0f}min_{round(1000*data.z[iz],0):.0f}m.png")


# iy = iy2
fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(12,8))

c1 = plot_cfill(data.xf, data.zf, data.dbz[0,:,iy,:], 'dbz', ax1, datalims=dbzlims)
ax1.set_ylabel('Height (km)')
ax1.set_title('$Z_H$', fontsize=14)
ax1.set_xlim(xl[0],xl[1])
ax1.set_ylim(0,10)

c2 = plot_cfill(data.xf, data.zf, data.winterp[0,:,iy,:], 'w', ax2, datalims=wlims)
ax2.set_title("$w$", fontsize=14)
ax2.set_xlim(xl[0],xl[1])
ax2.set_ylim(0,10)

c3 = plot_cfill(data.xf, data.zf, data.thpert[0,:,iy,:], 'thpert', ax3, datalims=pthlims)
#c3 = plot_cfill(data.xf, data.zf, data.prspert[0,:,iy,:]/100, 'prspert', ax3, datalims=pplims)
ax3.set_xlabel('x distance (km)')
ax3.set_ylabel('Height (km)')
ax3.set_title("\u03B8'", fontsize=14)
#ax3.set_title("p'", fontsize=14)
ax3.set_xlim(xl[0],xl[1])
ax3.set_ylim(0,10)

c4 = plot_cfill(data.xf, data.zf, data.zvort[0,:,iy,:], 'zvort', ax4, datalims=svlims)
ax4.set_xlabel('x distance (km)')
ax4.set_title("\u03B6", fontsize=14)
ax4.set_xlim(xl[0],xl[1])
ax4.set_ylim(0,10)

plt.suptitle(f"$t={data.time[0]/60:.0f}$ min", fontsize=14)
plt.show()
if figsave:
    plt.savefig(fp+f"dataV_zoom_{data.time[0]/60:.0f}min.png")



#%%


fig = plt.figure()
plt.plot(stats.time, stats.vortsfc, '-k', lw=1)
plt.xlabel('Model time (s)',size=14)
plt.ylabel('Max sfc vertical vorticity',size=14)
plt.grid(b=True)
plt.show()
if figsave:
    plt.savefig(fp+f"vortsfc.png")

fig = plt.figure()
plt.plot(stats.time, stats.wmax500, '-k', lw=1)
plt.xlabel('Model time (s)',size=14)
plt.ylabel('Max low-level w',size=14)
plt.grid(b=True)
plt.show()
if figsave:
    plt.savefig(fp+f"wmax500.png")



#%%

tapfrq = 300.
timax = 28800.
dtime = 150.*60.

nf = timax/tapfrq + 1

t = np.array([ tapfrq*i for i in np.arange(nf) ])

fnum = np.nonzero(t==dtime)[0][0] + 1
print(f"cm1out_{fnum:06d}")


#%%

fp = 'D:/cm1out/merger/supercell2-500m/'
#stats = read_cm1out(fp+'cm1out_stats.nc')

ff = glob(fp+'cm1out_0000*.nc')
files = ff[0:len(ff):2]


dsvars = ['time','zf','zvort','winterp']
data = [ structtype() for n in range(len(files)) ]
for n in range(len(files)):
    print(f"Reading file {files[n]} ...")
    data[n] = read_cm1out(files[n], dsvars)
print("Finished reading files!")

#%%

time = np.array([data[i].time[0] for i in range(len(data))])
zf = np.array(data[0].zf)
zvort_max = np.array([np.max(data[i].zvort[0,:,:,:], axis=(1,2)) for i in range(len(data))])
w_max = np.array([np.max(data[i].winterp[0,:,:,:], axis=(1,2)) for i in range(len(data))])

#%%
fig, (ax1,ax2) = plt.subplots(2,1,figsize=(12,8))

c1 = plot_cfill(time/60, zf, w_max[1::,:].transpose(), 'w', ax1, cmap='Reds')
ax1.set_ylabel('Height (km)')
ax1.set_title('Maximum $w$', fontsize=14)

c2 = plot_cfill(time/60, zf, zvort_max[1::,:].transpose(), 'zvort', ax2, cmap='Reds')
ax2.set_xlabel('Time (min)')
ax2.set_ylabel('Height (km)')
ax2.set_title("Maximum \u03B6", fontsize=14)


#plt.suptitle(f"$t={data.time[0]:.0f}$ s\n$h_{{CP}}=4$ km, ${chr(952)}'_{{CP}}=-8$ K", fontsize=14)
plt.show()
#if figsave:
    #plt.savefig(fp+f"dataV_zoom_{data.time[0]:.0f}s_hook.png")



