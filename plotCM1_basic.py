# -*- coding: utf-8 -*-
"""
Morgan Schneider

Test code for plotting CM1 output.
"""

####################
### Load modules ###
####################

from CM1utils import *
from matplotlib.animation import FuncAnimation


#%% Read data

figsave = True
plot4p = True

#fp = '/Users/morgan.schneider/Documents/cm1/merger/suptest/'
fp = 'D:/cm1out/merger/merger2-500m/'
stats = read_cm1out(fp+'cm1out_stats.nc')


ff = glob(fp+'cm1out_0000*.nc')
files = ff[0:len(ff):3]

#ic = read_cm1out(files[0], ['th0','prs0','qv0','u0','v0'])

dsvars = ['time','xh','yh','z','xf','yf','zf','thpert','dbz','zvort','winterp']
data = [ structtype() for n in range(len(files)) ]
for n in range(len(files)):
    print(f"Reading file {files[n]} ...")
    data[n] = read_cm1out(files[n], dsvars)
print("Finished reading files!")


#%% Make plots

#plt.close('all')


fig1 = plt.figure()
plt.plot(stats.time, stats.cflmax, '-k', lw=1)
plt.xlabel('Model time (s)',size=16)
plt.ylabel('Max CFL',size=16)
plt.title('Model stability',size=18)
plt.grid(b=True)
plt.show()

if figsave:
    plt.savefig(fp+'cfl.png')



#%% More plots

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

lev = 3.
iz = np.where(data[0].z >= lev)[0][0]
iy = np.where(data[0].yh >= 27.)[0][0]
zlims = [0, 12]

#%%

if False:
    ii = 12
    
    fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(12,8))
    
    c1 = plot_cfill(data[ii].xf, data[ii].yf, data[ii].dbz[0,0,:,:], 'dbz', ax1, datalims=dbzlims)
    ax1.set_ylabel('y distance (km)')
    ax1.set_title(f"$Z_H$ at {1000*data[ii].z[0]:.0f} m", fontsize=14)
    rect = patches.Rectangle((80,-20), 40, 40, linewidth=2, edgecolor='k', facecolor='none')
    ax1.add_patch(rect)
    
    c2 = plot_cfill(data[ii].xf, data[ii].yf, np.mean(data[ii].winterp[0,:,:,:],0), 'w', ax2, datalims=wlims)
    ax2.set_title(f"Column-mean $w$", fontsize=14)
    rect = patches.Rectangle((80,-20), 40, 40, linewidth=2, edgecolor='k', facecolor='none')
    ax2.add_patch(rect)
    
    c3 = plot_cfill(data[ii].xf, data[ii].yf, data[ii].thpert[0,0,:,:], 'thpert', ax3, datalims=pthlims)
    ax3.set_xlabel('x distance (km)')
    ax3.set_ylabel('y distance (km)')
    ax3.set_title(f"\u03B8' at {1000*data[ii].z[0]:.0f} m", fontsize=14)
    rect = patches.Rectangle((80,-20), 40, 40, linewidth=2, edgecolor='k', facecolor='none')
    ax3.add_patch(rect)
    
    c4 = plot_cfill(data[ii].xf, data[ii].yf, data[ii].zvort[0,iz,:,:], 'zvort', ax4, datalims=svlims)
    ax4.set_xlabel('x distance (km)')
    ax4.set_title(f"\u03B6 at {1000*data[ii].z[0]:.0f} m", fontsize=14)
    rect = patches.Rectangle((80,-20), 40, 40, linewidth=2, edgecolor='k', facecolor='none')
    ax4.add_patch(rect)
    
    plt.suptitle(f"$t={data[ii].time[0]:.0f}$ s\n$h_{{CP}}=4$ km, ${chr(952)}'_{{CP}}=-8$ K", fontsize=14)
    plt.show()
    #plt.savefig(fp+f"dataH_{data[ii].time[0]:.0f}s.png")
    
    
    
    fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(12,8))
    
    #c1 = plot_cfill(data[ii].xf, data[ii].zf, np.mean(data[ii].dbz[0,:,:,:],1), 'dbz', ax1, datalims=dbzlims)
    c1 = plot_cfill(data[ii].xf, data[ii].zf, data[ii].dbz[0,:,iy,:], 'dbz', ax1, datalims=dbzlims)
    ax1.set_ylabel('Height (km)')
    ax1.set_title('Meridional mean $Z_H$', fontsize=14)
    ax1.set_ylim(zlims[0],zlims[1])
    
    #c2 = plot_cfill(data[ii].xf, data[ii].zf, np.mean(data[ii].winterp[0,:,:,:],1), 'w', ax2, datalims=wlims)
    c2 = plot_cfill(data[ii].xf, data[ii].zf, data[ii].winterp[0,:,iy,:], 'w', ax2, datalims=wlims)
    ax2.set_title("Meridional mean $w$", fontsize=14)
    ax2.set_ylim(zlims[0],zlims[1])
    
    #c3 = plot_cfill(data[ii].xf, data[ii].zf, np.mean(data[ii].thpert[0,:,:,:],1), 'pth', ax3, datalims=pthlims)
    c3 = plot_cfill(data[ii].xf, data[ii].zf, data[ii].thpert[0,:,iy,:], 'thpert', ax3, datalims=pthlims)
    ax3.set_xlabel('x distance (km)')
    ax3.set_ylabel('Height (km)')
    ax3.set_title("Meridional mean \u03B8'", fontsize=14)
    ax3.set_ylim(zlims[0],zlims[1])
    
    #c4 = plot_cfill(data[ii].xf, data[ii].zf, np.mean(data[ii].zvort[0,:,:,:],1), 'zvort', ax4, datalims=svlims)
    c4 = plot_cfill(data[ii].xf, data[ii].zf, data[ii].zvort[0,:,iy,:], 'zvort', ax4, datalims=svlims)
    ax4.set_xlabel('x distance (km)')
    ax4.set_title("Meridional mean \u03B6", fontsize=14)
    ax4.set_ylim(zlims[0],zlims[1])
    
    plt.suptitle(f"$t={data[ii].time[0]:.0f}$ s\n$h_{{CP}}=4$ km, ${chr(952)}'_{{CP}}=-8$ K", fontsize=14)
    plt.show()
    #plt.savefig(fp+f"dataV_{data[ii].time[0]:.0f}s.png")
    
    
    
    fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(12,8))
    
    c1 = plot_cfill(data[ii].xf, data[ii].yf, data[ii].dbz[0,0,:,:], 'dbz', ax1, datalims=dbzlims)
    ax1.plot([95, 115], [data[ii].yh[iy], data[ii].yh[iy]], 'k', linewidth=2)
    ax1.set_ylabel('y distance (km)')
    ax1.set_title(f"$Z_H$ at {1000*data[ii].z[0]:.0f} m", fontsize=14)
    ax1.set_xlim(90,120)
    ax1.set_ylim(-20,20)
    
    c2 = plot_cfill(data[ii].xf, data[ii].yf, np.mean(data[ii].winterp[0,:,:,:],0), 'w', ax2, datalims=wlims)
    ax2.plot([95, 115], [data[ii].yh[iy], data[ii].yh[iy]], 'k', linewidth=2)
    ax2.set_title(f"Column-mean $w$", fontsize=14)
    ax2.set_xlim(90,120)
    ax2.set_ylim(-20,20)
    
    c3 = plot_cfill(data[ii].xf, data[ii].yf, data[ii].thpert[0,0,:,:], 'thpert', ax3, datalims=pthlims)
    ax3.plot([95, 115], [data[ii].yh[iy], data[ii].yh[iy]], 'k', linewidth=2)
    ax3.set_xlabel('x distance (km)')
    ax3.set_ylabel('y distance (km)')
    ax3.set_title(f"\u03B8' at {1000*data[ii].z[0]:.0f} m", fontsize=14)
    ax3.set_xlim(90,120)
    ax3.set_ylim(-20,20)
    
    c4 = plot_cfill(data[ii].xf, data[ii].yf, data[ii].zvort[0,0,:,:], 'zvort', ax4, datalims=svlims)
    ax4.plot([95, 115], [data[ii].yh[iy], data[ii].yh[iy]], 'k', linewidth=2)
    ax4.set_xlabel('x distance (km)')
    ax4.set_title(f"\u03B6 at {1000*data[ii].z[0]:.0f} m", fontsize=14)
    ax4.set_xlim(90,120)
    ax4.set_ylim(-20,20)
    
    plt.suptitle(f"$t={data[ii].time[0]:.0f}$ s\n$h_{{CP}}=4$ km, ${chr(952)}'_{{CP}}=-8$ K", fontsize=14)
    plt.show()
    
    
    
    fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(12,8))
    
    c1 = plot_cfill(data[ii].xf, data[ii].zf, data[ii].dbz[0,:,iy,:], 'dbz', ax1, datalims=dbzlims)
    ax1.set_ylabel('Height (km)')
    ax1.set_title('$Z_H$', fontsize=14)
    ax1.set_xlim(90,120)
    ax1.set_ylim(0,8)
    
    c2 = plot_cfill(data[ii].xf, data[ii].zf, data[ii].winterp[0,:,iy,:], 'w', ax2, datalims=wlims)
    ax2.set_title("$w$", fontsize=14)
    ax2.set_xlim(90,120)
    ax2.set_ylim(0,8)
    
    c3 = plot_cfill(data[ii].xf, data[ii].zf, data[ii].thpert[0,:,iy,:], 'thpert', ax3, datalims=pthlims)
    ax3.set_xlabel('x distance (km)')
    ax3.set_ylabel('Height (km)')
    ax3.set_title("\u03B8'", fontsize=14)
    ax3.set_xlim(90,120)
    ax3.set_ylim(0,8)
    
    c4 = plot_cfill(data[ii].xf, data[ii].zf, data[ii].zvort[0,:,iy,:], 'zvort', ax4, datalims=svlims)
    ax4.set_xlabel('x distance (km)')
    ax4.set_title("\u03B6", fontsize=14)
    ax4.set_xlim(90,120)
    ax4.set_ylim(0,8)
    
    plt.suptitle(f"$t={data[ii].time[0]:.0f}$ s\n$h_{{CP}}=4$ km, ${chr(952)}'_{{CP}}=-8$ K", fontsize=14)
    plt.show()
    
    
    
    
    
#%%
    
if False:
    th0 = ic.th0[0,:,0,0]
    th = data[0].th[0,:,0,0]
    prs0 = ic.prs0[0,:,0,0]
    qv0 = ic.qv0[0,:,0,0];
    u0 = ic.u0[0,:,0,0]
    v0 = ic.v0[0,:,0,0]
    
    T = th0*(prs0/100000.)**0.286
    e = (qv0*prs0/100)/(0.622+qv0)
    Td=243.5/((17.67/(np.log(e/6.112)))-1)
    T_parcel = mc.parcel_profile((prs0/100.)*units.hPa, T[0]*units.K, (Td[0]+273.15)*units.K)
    
    fig = plt.figure(figsize=(8,8))
    
    skew = SkewT(fig=fig)
    skew.plot(prs0/100., (T-273.15), '-r', linewidth=2)
    skew.plot(prs0/100., Td, '-g', linewidth=2)
    skew.plot(prs0/100., np.array(T_parcel.magnitude[:])-273.15, '-k', linewidth=2)
    skew.plot_dry_adiabats()
    skew.plot_moist_adiabats()
    skew.plot_mixing_lines()
    skew.ax.set_ylim(1000, 100)
    skew.ax.set_xlim(-40, 30)
    plt.title('Initial condition profile (domain-relative)')
    ax_hod = inset_axes(skew.ax, '42%', '42%', loc=1)
    H = Hodograph(ax_hod, component_range=40.)
    H.add_grid(increment=20)
    H.plot(u0, v0, color='k', linewidth=1.5)
    
    plt.show()
    plt.savefig(fp+'IC_profile.png', dpi=300)




#%% oops! all plots

def animateH(i):
    c1.set_array(data[i].dbz[0,0,:,:])
    c2.set_array(data[i].winterp[0,iz,:,:])
    c3.set_array(data[i].thpert[0,0,:,:])
    c4.set_array(data[i].zvort[0,iz,:,:])
    plt.suptitle(f"$t={data[i].time[0]/60:.0f}$ min", fontsize=14)

if plot4p:
    fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(12,8))
    
    c1 = plot_cfill(data[0].xf, data[0].yf, data[0].dbz[0,0,:,:], 'dbz', ax1, datalims=dbzlims)
    ax1.set_ylabel('y distance (km)')
    ax1.set_title(f"$Z_H$ at {1000*data[0].z[0]:.0f} m AGL", fontsize=12)
    
    #c2 = plot_cfill(data[0].xf, data[0].yf, np.mean(data[0].winterp[0,:,:,:],0), 'w', ax2, datalims=wlims)
    #ax2.set_title("Column-average $w$", fontsize=12)
    c2 = plot_cfill(data[0].xf, data[0].yf, data[0].winterp[0,iz,:,:], 'w', ax2, datalims=wlims)
    ax2.set_title(f"$w$ at {1000*data[0].z[iz]:.0f} m AGL", fontsize=12)
    
    c3 = plot_cfill(data[0].xf, data[0].yf, data[0].thpert[0,0,:,:], 'thpert', ax3, datalims=pthlims)
    ax3.set_xlabel('x distance (km)')
    ax3.set_ylabel('y distance (km)')
    ax3.set_title(f"\u03B8' at {1000*data[0].z[0]:.0f} m AGL", fontsize=12)
    
    c4 = plot_cfill(data[0].xf, data[0].yf, data[0].zvort[0,iz,:,:], 'zvort', ax4, datalims=svlims)
    ax4.set_xlabel('x distance (km)')
    ax4.set_title(f"\u03B6 at {1000*data[0].z[iz]:.0f} m AGL", fontsize=12)
    
    anim = FuncAnimation(fig, animateH, frames=len(data), interval=100, repeat=False, blit=False)
    if figsave:
        anim.save(fp+'animH_4p.gif')
    plt.show()



#%%

# figure out iy indices to slice through the supercell
def animateV(i):
    c1.set_array(np.mean(data[i].dbz[0,:,:,:],1))
    c2.set_array(np.mean(data[i].winterp[0,:,:,:],1))
    c3.set_array(np.mean(data[i].thpert[0,:,:,:],1))
    c4.set_array(np.mean(data[i].zvort[0,:,:,:],1))
    plt.suptitle(f"$t={data[i].time[0]/60:.0f}$ min", fontsize=14)

if plot4p:
    fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(12,8))
    
    c1 = plot_cfill(data[0].xf, data[0].zf, np.mean(data[0].dbz[0,:,:,:],1), 'dbz', ax1, datalims=dbzlims)
    ax1.set_ylabel('Height (km)')
    ax1.set_title('Meridional mean $Z_H$', fontsize=12)
    # ax1.set_xlim(xlims[0],xlims[1])
    ax1.set_ylim(zlims[0],zlims[1])
    
    c2 = plot_cfill(data[0].xf, data[0].zf, np.mean(data[0].winterp[0,:,:,:],1), 'w', ax2, datalims=wlims)
    ax2.set_title("Meridional mean $w$", fontsize=12)
    # ax2.set_xlim(xlims[0],xlims[1])
    ax2.set_ylim(zlims[0],zlims[1])
    
    c3 = plot_cfill(data[0].xf, data[0].zf, np.mean(data[0].thpert[0,:,:,:],1), 'thpert', ax3, datalims=pthlims)
    ax3.set_xlabel('x distance (km)')
    ax3.set_ylabel('Height (km)')
    ax3.set_title("Meridional mean \u03B8'", fontsize=12)
    # ax3.set_xlim(xlims[0],xlims[1])
    ax3.set_ylim(zlims[0],zlims[1])
    
    c4 = plot_cfill(data[0].xf, data[0].zf, np.mean(data[0].zvort[0,:,:,:],1), 'zvort', ax4, datalims=svlims)
    ax4.set_xlabel('x distance (km)')
    ax4.set_title("Meridional mean \u03B6", fontsize=12)
    # ax4.set_xlim(xlims[0],xlims[1])
    ax4.set_ylim(zlims[0],zlims[1])
    
    anim = FuncAnimation(fig, animateV, frames=len(data), interval=1000, repeat=False, blit=False)
    if figsave:
        anim.save(fp+'animV_4p_0.gif')
    plt.show()



def animateVort(i):
    c1.set_array(data[i].xvort[0,0,:,:])
    c2.set_array(data[i].yvort[0,0,:,:])
    c3.set_array(data[i].zvort[0,0,:,:])
    # c4.set_array(np.mean(data[i].xvort[0,:,:,:],1))
    # c5.set_array(np.mean(data[i].yvort[0,:,:,:],1))
    # c6.set_array(np.mean(data[i].zvort[0,:,:,:],1))
    plt.suptitle(f"t = {data[i].time[0]/60:.0f} min", fontsize=14)
    
if False:
    fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(16,4))
    
    c1 = plot_cfill(data[0].xf, data[0].yf, data[0].xvort[0,0,:,:], 'hvort', ax1, datalims=svlims)
    ax1.set_xlabel('x distance (km)')
    ax1.set_ylabel('y distance (km)')
    ax1.set_title(f"x vorticity at {1000*data[0].z[0]:.0f} m", fontsize=12)
    
    c2 = plot_cfill(data[0].xf, data[0].yf, data[0].yvort[0,0,:,:], 'hvort', ax2, datalims=svlims)
    ax2.set_xlabel('x distance (km)')
    #ax2.set_ylabel('y distance (km)')
    ax2.set_title(f"y vorticity at {1000*data[0].z[0]:.0f} m", fontsize=12)
    
    c3 = plot_cfill(data[0].xf, data[0].yf, data[0].zvort[0,0,:,:], 'zvort', ax3, datalims=svlims)
    ax3.set_xlabel('x distance (km)')
    #ax3.set_ylabel('y distance (km)')
    ax3.set_title(f"z vorticity at {1000*data[0].z[0]:.0f} m", fontsize=12)
    
    # c4 = plot_cfill(data[0].xf, data[0].zf, np.mean(data[0].xvort[0,:,:,:],1), 'vort', ax4, datalims=svlims)
    # ax4.set_xlabel('x distance (km)')
    # ax4.set_ylabel('height (km)')
    # ax4.set_title(f"meridional mean x vorticity", fontsize=12)
    
    # c5 = plot_cfill(data[0].xf, data[0].zf, np.mean(data[0].yvort[0,:,:,:],1), 'vort', ax5, datalims=svlims)
    # ax5.set_xlabel('x distance (km)')
    # #ax5.set_ylabel('height (km)')
    # ax5.set_title(f"meridional mean y vorticity", fontsize=12)
    
    # c6 = plot_cfill(data[0].xf, data[0].zf, np.mean(data[0].zvort[0,:,:,:],1), 'vort', ax6, datalims=svlims)
    # ax6.set_xlabel('x distance (km)')
    # #ax6.set_ylabel('height (km)')
    # ax6.set_title(f"meridional mean vertical vorticity", fontsize=12)
    
    anim = FuncAnimation(fig, animateVort, frames=len(data), interval=1000, repeat=False, blit=False)
    if figsave:
        anim.save(fp+'vort_anim.gif')
    plt.show()


#%%

plevs = [1, 2, 4, 6, 8, 10]
nlevs = [-10, -8, -6, -4, -2, -0.5]


if True:
    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(10,4))
    
    a1 = plot_cfill(data[0].xf, data[0].zf, np.mean(data[0].thpert[0,:,:,:],1), 'thpert', ax1, datalims=pthlims)
    b1 = ax1.contour(data[0].xh, data[0].z, np.mean(data[0].winterp[0,:,:,:],1), levels=plevs, colors='k', linewidths=1)
    c1 = ax1.contour(data[0].xh, data[0].z, np.mean(data[0].winterp[0,:,:,:],1), levels=nlevs, colors='k', linewidths=1, linestyles='dashed')
    ax1.set_title("Meridional average \u03B8'")
    ax1.set_xlabel('x distance (km)'); ax1.set_ylabel('Height (km)')
    ax1.set_xlim(-10,150)
    ax1.set_ylim(zlims[0],zlims[1])
    
    a2 = plot_cfill(data[0].xf, data[0].zf, np.mean(data[0].dbz[0,:,:,:],1), 'dbz', ax2, datalims=dbzlims)
    b2 = ax2.contour(data[0].xh, data[0].z, np.mean(data[0].winterp[0,:,:,:],1), levels=plevs, colors='k', linewidths=1)
    c2 = ax2.contour(data[0].xh, data[0].z, np.mean(data[0].winterp[0,:,:,:],1), levels=nlevs, colors='w', linewidths=1, linestyles='dashed')
    ax2.set_title("Meridional average $Z_H$ and w")
    ax2.set_xlabel('x distance (km)')
    ax2.set_xlim(-10,150)
    ax2.set_ylim(zlims[0],zlims[1])
    
    
    def anim2p(i):
        global a1,a2
        ax1.clear()
        ax2.clear()
        a1.colorbar.remove()
        a2.colorbar.remove()

        a1 = plot_cfill(data[0].xf, data[0].zf, np.mean(data[i].thpert[0,:,:,:],1), 'thpert', ax1, datalims=pthlims)
        b1 = ax1.contour(data[0].xh, data[0].z, np.mean(data[i].winterp[0,:,:,:],1), levels=plevs, colors='k', linewidths=1)
        c1 = ax1.contour(data[0].xh, data[0].z, np.mean(data[i].winterp[0,:,:,:],1), levels=nlevs, colors='k', linewidths=1, linestyles='dashed')
        ax1.set_title("Meridional average \u03B8'")
        ax1.set_xlabel('x distance (km)'); ax1.set_ylabel('Height (km)')
        ax1.set_xlim(-10,150)
        ax1.set_ylim(zlims[0],zlims[1])

        a2 = plot_cfill(data[0].xf, data[0].zf, np.mean(data[i].dbz[0,:,:,:],1), 'dbz', ax2, datalims=dbzlims)
        b2 = ax2.contour(data[0].xh, data[0].z, np.mean(data[i].winterp[0,:,:,:],1), levels=plevs, colors='k', linewidths=1)
        c2 = ax2.contour(data[0].xh, data[0].z, np.mean(data[i].winterp[0,:,:,:],1), levels=nlevs, colors='w', linewidths=1, linestyles='dashed')
        ax2.set_title("Meridional average $Z_H$ and w")
        ax2.set_xlabel('x distance (km)')
        ax2.set_xlim(-10,150)
        ax2.set_ylim(zlims[0],zlims[1])
        
        plt.suptitle(f"t = {data[i].time[0]/60:.0f} min", fontsize=14)
        
    
    anim = FuncAnimation(fig, anim2p, frames=len(data), interval=1000, repeat=False, blit=False)
    if figsave:
        anim.save(fp+'anim2p_theta+dbz.gif')
    plt.show()




#%% EXPLANATION OF VARIABLES:

'''
############################
###   cm1out_0000**.nc   ###
############################

### NOTE: horizontal wind speeds are grid-relative
### LML = lowest model level

# time:         model time in s
# ztop:         height of domain top in km
# umove:        zonal domain translation speed in m/s
# vmove:        meridional domain translation speed in m/s
# xh (nx,):     scalar x grid point locations in km (x grid centers)
# yh (ny,):     scalar y grid point locations in km (x grid centers)
# z (nz,):      scalar z grid point locations in km (x grid centers)
# xf (nx+1,):   staggered u grid point locations in km (x grid corners)
# yf (ny+1,):   staggered v grid point locations in km (y grid corners)
# zf (nz+1,):   staggered w grid point locations in km (z grid corners)

# u0 (time, zh, yh, xf):        base state u in m/s
# v0 (time, zh, yf, xh):        base state v in m/s
# w0 (time, zf, yh, xh):        base state w in m/s
# prs0 (time, zh, yh, xh):      base state pressure in Pa
# th0 (time, zh, yh, xh):       base state potential temperature in K
# pi0 (time, zh, yh, xh):       base state nondimensional pressure
# qv0 (time, zh, yh, zh):       base state water vapor mixing ratio in kg/kg

# u (time, zh, yh, xf):         u velocity in m/s
# v (time, zh, yf, xh):         v velocity in m/s
# w (time, zf, yh, xh):         w velocity in m/s
# uinterp (time, zh, yh, xh):   u interpolated to scalar grid
# vinterp (time, zh, yh, xh):   v interpolated to scalar grid
# winterp (time, zh, yh, xh):   w interpolated to scalar grid

# dbz (time, zh, yh, xh):       reflectivity in dBZ
# cref (time, yh, xh):          composite reflectivity (max in column) in dBZ
# th (time, zh, yh, xh):        potential temperature in K
# thpert(time, zh, yh, xh):     potential temperature perturbation in K
# prs (time, zh, yh, xh):       pressure in Pa
# prspert (time, zh, yh, xh):   pressure perturbation in Pa
# pi (time, zh, yh, xh):        nondimensional pressure (Exner function)
# pipert (time, zh, yh, xh):    nondimensional pressure perturbation
# rho (time, zh, yh, xh):       dry air density in kg/m^3
# rhopert (time, zh, yh, xh):   dry air density perturbation in kg/m^3
# xvort (time, zh, yh, xh):     horizontal x vorticity in 1/s
# yvort (time, zh, yh, xh):     horizontal y vorticity in 1/s
# zvort (time, zh, yh, xh):     vertical vorticity in 1/s

# qv (time, zh, yh, xh):        water vapor mixing ratio in kg/kg
# qvpert (time, zh, yh, xh):    water vapor mixing ratio perturbation in kg/kg
# qc (time, zh, yh, xh):        cloud water mixing ratio in kg/kg
# qr (time, zh, yh, xh):        rain water mixing ratio in kg/kg
# qi (time, zh, yh, xh):        cloud ice mixing ratio in kg/kg
# qs (time, zh, yh, xh):        snow mixing ratio in kg/kg
# qg (time, zh, yh, xh):        graupel mixing ratio in kg/kg
# nci (time, zh, yh, xh):       cloud ice number concentration in #/m^3
# ncs (time, zh, yh, xh):       snow number concentration in #/m^3
# ncr (time, zh, yh, xh):       rain drop number concentration in #/m^3
# ncg (time, zh, yh, xh):       graupel number concentration in #/m^3

# tke (time, zf, yh, xh):       subgrid TKE in m2/s2
# kmh (time, zf, yh, xh):       horizontal eddy viscosity for momentum in m2/s 
# kmv (time, zf, yh, xh):       vertical eddy viscosity for momentum in m2/s
# khh (time, zf, yh, xh):       horizontal eddy diffusivity for scalars in m2/s
# khv (time, zf, yh, xh):       vertical eddy diffusivity for scalars in m2/s

# rain (time, yh, xh):          accumulated sfc rainfall in cm
# prate (time, yh, xh):         sfc precip rate in kg/m2/s
# sws (time, yh, xh):           max LML horizontal wind speed in m/s
# svs (time, yh, xh):           max LML vertical vorticity in 1/s
# sps (time, yh, xh):           min LML pressure in Pa
# srs (time, yh, xh):           max LML rain water mixing ratio in kg/kg
# sgs (time, yh, xh):           max LML graupel mixing ratio in kg/kg
# sus (time, yh, xh):           max w at 5 km AGL in m/s
# shs (time, yh, xh):           max integrated updraft helicity in m2/s2
# cpc (time, yh, xh):           cold pool intensity C in m/s
# cph (time, yh, xh):           cold pool depth H in m AGL

# rain2 (time, yh, xh):         translated rain (imove = 1)
# prate2 (time, yh, xh):        translated prate
# sws2 (time, yh, xh):          translated sws
# svs2 (time, yh, xh):          translated svs
# sps2 (time, yh, xh):          translated sps
# srs2 (time, yh, xh):          translated srs
# sgs2 (time, yh, xh):          translated sgs
# sus2 (time, yh, xh):          translated sus
# shs2 (time, yh, xh):          translated shs







###########################
###   cm1out_stats.nc   ###
###########################


### NOTE: all variables have dimension (time,) except xh/yh/zh which are scalar
### NOTE: horizontal wind speeds are grid-relative
### LML = lowest model level

# time:     time in s? not sure how this is different from mtime. maybe for restarts?
# mtime:    model time in s
# xh:       zonal location in degrees E
# yh:       meridional location in degrees N
# z:        height in m

# umax/umin:            max/min u in m/s
# vmax/vmin:            max/min v in m/s
# wmax*/wmin*:          max/min w in m/s      (also wmax/wmin 500/1000/2500/5000/10k)
# zwmax/zwmin:          level of max/min w in m AGL
# sumax/sumin:          max/min LML u in m/s
# svmax/svmin:          max/min LML v in m/s

# ppimax/ppimin:        max/min nondimensional p'
# ppmax/ppmin:          max/min p' in Pa
# thpmax/thpmin:        max/min theta' in K
# sthpmax/sthpmin:      max/min LML theta' in K
# maxq*/minq*:          max/min mixing ratio in kg/kg     (qv,qc,qr,qi,qs,qg)
# maxnc*/minnc*:        max/min number concentration in #/kg    (nci,ncs,ncr,ncg)
# pratemax/pratemin:    max/min sfc precip rate in kg/m2/s
# tkemax/tkemin:        max/min subgrid TKE in m2/s2
# km*max/km*min:        max/min eddy viscosity for momentum in m2/s (kmh,kmv)
# kh*max/kh*min:        max/min eddy diffusivity for scalars in m2/s (khh,khv)

# divmax/divmin:        max/min 3D divergence in 1/s
# rhmax/rhmin:          max/min RH w.r.t. liquid
# rhimax/rhimin:        max/min RH w.r.t. ice
# themax/themin:        max/min theta-e below 10 km in K
# sthemax/sthemin:      max/min LML theta-e in K
# sprsmax/sprsmin:      max/min LML pressure in Pa
# psfcmax/psfcmin:      max/min sfc pressure in Pa
# wspmax/wspmin:        max/min horiz wind speed in m/s
# zwspmax/zwspmin:      level of max/min horiz wind speed in m AGL
# swspmax/swspmin:      max/min LML horiz wind speed in m/s

# vortsfc:      max LML vertical vorticity in 1/s
# vort*km:      max vertical vorticity in 1/s   (1km,2km,3km,4km,5km)
# qctop:        max cloud top height in m
# qcbot:        min cloud base height in m

# cflmax:       max Courant number
    # The Courant number is a dimensionless value representing the time a
    # particle stays in one cell of the mesh. It must be below 1 and should
    # ideally be below 0.7. If the Courant number exceeds 1, the time step is
    # too large to see the particle in one cell, it “skips” the cell.

'''



