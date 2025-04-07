# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 12:34:41 2023
@author: morgan.schneider

Plotting RaXPol data
"""

####################
### Load modules ###
####################

from RaxpolUtils import *
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#%% IOP2 single time with mobile mesonet

figsave = False

fp = '/Users/morgan.schneider/Documents/perils2023/iop2/'
filetime = '081928'
fn = glob(fp+f"raxpol/CFradial_cal/*_{filetime}_*.nc")[0]

rax = read_raxpol(fn)
raxpol = pyart.io.read(fn)

P1 = read_MM_dat(fp+'mesonet/Probe_1_IOP2_QC_all.dat')
P2 = read_MM_dat(fp+'mesonet/Probe_2_IOP2_QC_all.dat')
xx1,yy1 = latlon2xy(P1['lat'], P1['lon'], raxpol.latitude['data'][0], raxpol.longitude['data'][0])
xx2,yy2 = latlon2xy(P2['lat'], P2['lon'], raxpol.latitude['data'][0], raxpol.longitude['data'][0])
P1.update({'xx':xx1, 'yy':yy1})
P2.update({'xx':xx2, 'yy':yy2})

elev = np.mean(raxpol.elevation['data']).round(1)

i1 = np.where(P1['time'] == int(filetime))[0][0]
i2 = np.where(P2['time'] == int(filetime))[0][0]

T_lims = [18,21] # Temperature
Td_lims = [17,19] # Dewpoint
datalims = T_lims


c = pyart.graph.RadarDisplay(raxpol)

fig = plt.figure(figsize=(9,6))
ax = fig.add_subplot(111)
# c.plot('DBZ', 0, vmin=0, vmax=70, cmap='pyart_NWSRef')
c.plot('VEL', 0, vmin=-30, vmax=30, cmap='pyart_Carbone42')
# c.plot_range_rings([2.5, 5, 7.5, 10])
ax.set_xlim([-4,3])
ax.set_ylim([-1,6])
plt.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')
# b1 = ax.barbs(P1['xx'][i1], P1['yy'][i1], P1['uwind'][i1], P1['vwind'][i1], barbcolor='k', length=9)
# b2 = ax.barbs(P2['xx'][i2], P2['yy'][i2], P2['uwind'][i2], P2['vwind'][i2], barbcolor='k', length=9)
s1 = ax.scatter(P1['xx'][i1], P1['yy'][i1], s=150, c=P1['Utube'][i1], cmap=cmaps['temp']['cm'],
                vmin=datalims[0], vmax=datalims[1], marker='s', edgecolors='k')
s2 = ax.scatter(P2['xx'][i2], P2['yy'][i2], s=175, c=P2['Utube'][i2], cmap=cmaps['temp']['cm'],
                vmin=datalims[0], vmax=datalims[1], marker='^', edgecolors='k')
ax.scatter(0, 0, s=30, c='k')
plt.colorbar(s1,label='Sfc temperature (C)')
plt.legend([s1,s2], ['Probe 1', 'Probe 2'], loc='upper right')
plt.show()

if figsave:
    plt.savefig(fp+f"figs/dbz+T_{filetime}.png", dpi=400)




#%% IOP5 single time with skyler

# fp = 'Documents/perils2023/iop5/raxpol/CFradial/'
# filetime = '160606'
# fn = glob(fp+f"*_{filetime}_*.nc")[0]

# 163100 - 1065
# 163500 - 1140

fp = '/Users/morgan.schneider/Documents/perils2023/iop5/raxpol/CFradial/'
files = sorted(glob(fp+'*.nc'))
# fn = files[1140]
fn = fp+'cfrad.20230405_160705_RaXPol_v104_s4.nc'


raxpol = pyart.io.read(fn)

sky_lat = 35.4110
sky_lon = -90.9549
sky_xx,sky_yy = latlon2xy(sky_lat, sky_lon, raxpol.latitude['data'][0], raxpol.longitude['data'][0])

c = pyart.graph.RadarDisplay(raxpol)

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111)
c.plot('DBZ', 0, vmin=0, vmax=70, cmap='pyart_NWSRef')
c.plot_range_rings([10,20,30])
ax.set_xlim([-30,30])
ax.set_ylim([-30,30])
ax.scatter(0, 0, s=50, c='k')
ax.scatter(sky_xx, sky_yy, s=50, c='k')
plt.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')
plt.text(sky_xx-1, sky_yy+0.4, 'Skyler', fontsize=10, fontweight='bold')
plt.show()


fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111)
c.plot('VEL', 0, vmin=-30, vmax=30, cmap='pyart_Carbone42')
c.plot_range_rings([10,20,30])
ax.set_xlim([-30,30])
ax.set_ylim([-30,30])
ax.scatter(0, 0, s=50, c='k')
ax.scatter(sky_xx, sky_yy, s=50, c='k')
plt.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')
plt.text(sky_xx-1, sky_yy+0.4, 'Skyler', fontsize=10, fontweight='bold')
plt.show()






#%% IOP2 load data into volume dicts

# from RaXPolUtils import *

# Pre-transect:       File 381 (index 380) 0.9 deg, 081758 UTC
# Transect start:     File 411 (index 410) 0.9 deg, 081858 UTC
# Outflow at P1:      File 462 (index 461) 2.6 deg, 082100 UTC
# Reflectivity at P1: File 495 (index 494) 0.9 deg, 082228 UTC
# Transect end:       File 600 (index 599) 0.9 deg, 082558 UTC
# Post-transect:      File (index ) 0.9 deg, 082658 UTC

vol_nums = [78, 81, 84, 87, 90, 93, 96, 102, 105, 108, 111, 112, 115, 118, 121, 124, 127, 130, 133, 136, 139]
# v78 320-331 # 081558-081620
# v81 335-346 # 081628-081650
# v84 350-361 # 081658-081720
# v87 365-376 # 081728-181750
# v90 380-391 # 081758-081820
# v93 395-406 # 081828-181850
# v96 410-421 # 081858-181920 *Transect start - vol[6]
# v102 430-441 # 081958-082020
# v105 445-456 # 082028-082050
# v108 460-471 # 082058-082120 *Outflow at P1 - vol[9]
# v112 480-490 # 082158-082220
# v115 494-505 # 082228-082250 *Reflectivity at P1 - vol[11]
# v118 509-520 # 082258-082320
# v121 524-535 # 082328-082350
# v124 539-550 # 082358-082420
# v127 554-565 # 082428-082450
# v130 569-580 # 082458-082520
# v133 584-595 # 082528-082550
# v136 599-610 # 082558-082620 *Transect end - vol[18]
# v139 614-625 # 082628-082650

t_start = '081758'
transect_start = '081858'
outflow_at_MM = '082100'
precip_at_MM = '082228'
transect_end = '082558'
t_end = '082658'

fp = '/Users/morgan.schneider/Documents/perils2023/iop2/'
files = sorted(glob(fp+'raxpol/CFradial_cal/*.nc'))

vol = [dict() for i in range(len(vol_nums))]

# # for vol 112
# els_actual = [0.9, 2.4, 3.8, 5.4, 6.9,  8.4,  9.8,  11.4, 12.8, 14.4, 15.8, 17.3, 18.8]
# els_bad    = [6.4, 6.9, 8.0, 9.5, 11.2, 12.4, 13.7, 16.0,  6.1, 15.1,       17.3, 18.8]

for vn in range(len(vol_nums)):
    print(vn)
    inds = [i for i,s in enumerate(files) if f"v{vol_nums[vn]}" in s]
    # if vol_nums[vn] == 111:
    #     inds = [i for i,s in enumerate(files) if f"v{vol_nums[vn]}" in s or f"v{vol_nums[vn]+1}" in s]
    a = slice(inds[1], inds[-1]+2)
    # if vol_nums[vn] == 112:
    #     a = slice(inds[0], inds[-1]+2)
    f = files[a]
    # if vol_nums[vn] == 111:
    #     f = np.append(np.append(f[0:4], f[6:12]), f[14:])
    # if vol_nums[vn] == 112:
    #     f = np.append(np.append(f[0], f[2:6]), f[8:])
    
    dbz_tmp = np.zeros(shape=(len(f),360,1246,), dtype=float)
    vel_tmp = np.zeros(shape=(len(f),360,1246,), dtype=float)
    sw_tmp = np.zeros(shape=(len(f),360,1246,), dtype=float)
    zdr_tmp = np.zeros(shape=(len(f),360,1246,), dtype=float)
    rhohv_tmp = np.zeros(shape=(len(f),360,1246,), dtype=float)
    xx_tmp = np.zeros(shape=(len(f),360,1246,), dtype=float)
    yy_tmp = np.zeros(shape=(len(f),360,1246,), dtype=float)
    zz_tmp = np.zeros(shape=(len(f),360,1246,), dtype=float)
    azvort_tmp = np.zeros(shape=(len(f),360,1246,), dtype=float) # vertical pseudovorticity
    elvort_tmp = np.zeros(shape=(len(f),360,1246,), dtype=float) # horizontal pseudovorticity
    div_tmp = np.zeros(shape=(len(f),360,1246,), dtype=float) # horizontal pseudodivergence
    az_tmp = np.zeros(shape=(len(f),360,), dtype=float)
    el_tmp = np.zeros(shape=(len(f),), dtype=float)
    time_tmp = []
    fname_tmp = []
    
    for ii in range(len(f)):
        fn = f[ii]
        d = read_raxpol(fn)
        
        r = d['r']
        el_tmp[ii] = d['elev']
        time_tmp.append(fn[85:91])
        fname_tmp.append(fn)
        
        for ix in np.linspace(0,359,360):
            if ix in d['az'].round(0):
                ind = np.where(d['az'].round(0) == ix)[0][0]
                
                az_tmp[ii,int(ix)] = d['az'][ind]
                dbz_tmp[ii,int(ix),:] = d['dbz'][ind,:]
                vel_tmp[ii,int(ix),:] = d['vel'][ind,:]
                sw_tmp[ii,int(ix),:] = d['sw'][ind,:]
                zdr_tmp[ii,int(ix),:] = d['zdr'][ind,:]
                rhohv_tmp[ii,int(ix),:] = d['rhohv'][ind,:]
                xx_tmp[ii,int(ix),:] = d['xx'][ind,:]
                yy_tmp[ii,int(ix),:] = d['yy'][ind,:]
                zz_tmp[ii,int(ix),:] = d['zz'][ind,:]
                
                if ix == 0:
                    v1 = d['vel'][-1,:]
                    v2 = d['vel'][ind+1,:]
                elif ix == 359:
                    v1 = d['vel'][ind-1,:]
                    v2 = d['vel'][0,:]
                else:
                    v1 = d['vel'][ind-1,:]
                    v2 = d['vel'][ind+1,:]
                azvort_tmp[ii,int(ix),:] = 1/(r*1000) * (v2-v1)/2
            else:
                az_tmp[ii,int(ix)] = ix
                xx_tmp[ii,int(ix),:] = r * np.sin(ix*np.pi/180) * np.cos(d['elev']*np.pi/180)
                yy_tmp[ii,int(ix),:] = r * np.cos(ix*np.pi/180) * np.cos(d['elev']*np.pi/180)
                zz_tmp[ii,int(ix),:] = r * np.sin(d['elev']*np.pi/180)
    
    elvort_tmp = 1/(r*1000) * np.gradient(vel_tmp, el_tmp, axis=0)
    div_tmp = np.gradient(vel_tmp, r*1000, axis=2)
    
    vol[vn].update({'dbz':dbz_tmp, 'vel':vel_tmp, 'sw':sw_tmp, 'zdr':zdr_tmp, 'rhohv':rhohv_tmp,
                    'xx':xx_tmp, 'yy':yy_tmp, 'zz':zz_tmp, 'az':az_tmp, 'elev':el_tmp, 'r':r,
                    'zvort':azvort_tmp, 'hvort':elvort_tmp, 'div':div_tmp,
                    'scan_time':time_tmp, 'vol_num':vol_nums[vn], 'filename':fname_tmp, 'va':d['va'],
                    'lat':d['lat'], 'lon':d['lon']})
    
    
    

r = d['r']
rax_lat = d['lat']
rax_lon = d['lon']
va = d['va']

P1 = read_MM_dat(fp+'mesonet/Probe_1_IOP2_QC_all.dat')
P2 = read_MM_dat(fp+'mesonet/Probe_2_IOP2_QC_all.dat')
xx1,yy1 = latlon2xy(P1['lat'], P1['lon'], rax_lat, rax_lon)
xx2,yy2 = latlon2xy(P2['lat'], P2['lon'], rax_lat, rax_lon)
P1.update({'xx':xx1, 'yy':yy1})
P2.update({'xx':xx2, 'yy':yy2})



dbfile = open('/Users/morgan.schneider/Documents/perils2023/iop2/circuit_data.pkl', 'wb')
data = {'Rax':vol, 'P1':P1, 'P2':P2}
pickle.dump(data, dbfile)
dbfile.close()


#%% IOP2 plot reconstructed RHIs and azimuth-height cross sections

# Scatter plot of rotor(s) location on top of PPIs for each volume?
# Azimuth-height pcolor plot of velocity through MV/rotor
# Do I need to advection correct for the reconstructed RHIs?? Check storm motion
# If couplet translation is due to advection, storm motion is about 15 m/s to the ~NE

ip = '/Users/morgan.schneider/Documents/perils2023/iop2/figs/'

# saved raxpol and probe data from circuit
if 'vol' not in locals():
    dbfile = open('/Users/morgan.schneider/Documents/perils2023/iop2/circuit_data.pkl', 'rb')
    dat = pickle.load(dbfile)
    vol = dat['Rax']
    # P1 = dat['P1']
    # P2 = dat['P2']
    dbfile.close()


dbfile = open('/Users/morgan.schneider/Documents/perils2023/iop2/raxpol_vortex_locs.pkl', 'rb')
locs = pickle.load(dbfile)
dbfile.close()

# i1 = np.where(P1['time'] == int(vol[vi]['scan_time'][eli]))[0][0]
# i2 = np.where(P2['time'] == int(vol[vi]['scan_time'][eli]))[0][0]

T_lims = [18,21] # Temperature
Td_lims = [17,19] # Dewpoint
datalims = T_lims
vort_lim = 0.003
div_lim = 0.1



vi = 7
eli = 1
filetime = vol[vi]['scan_time'][eli]
vortex_num = 1
# vortex 1: vi = 7-9
# vortex 2: vi = 8-13
# vortex 3: vi = 8-17
# vortex 4: vi = 9-12
# rotor: vi = 8

az_rot = locs[filetime][f"vortex{vortex_num}"]['az']
r_rot = locs[filetime][f"vortex{vortex_num}"]['r']
z_rot = locs[filetime][f"vortex{vortex_num}"]['z']

# azimuth = round(np.mean(az_rot))
azimuth = az_rot[round(len(az_rot)/2)]
# if az_rot[-1] < az_rot[0]:
#     az_tmp = [az_rot[i]+360 if az_rot[i] < 180 else az_rot[i] for i in range(len(az_rot))]

# if azimuth not in az_rot:
#     print(f"Invalid azimuth, must be between {az_rot[0]}-{az_rot[-1]} degrees")
#     azimuth = az_rot[0]

azi = np.where(vol[ii]['az'][eli,:].round(0) == azimuth)[0][0]
rr = (vol[vi]['xx'][:,azi,:]**2 + vol[vi]['yy'][:,azi,:]**2)**0.5


x_rot = r_rot * np.sin(az_rot*np.pi/180)
y_rot = r_rot * np.cos(az_rot*np.pi/180)
irot = np.where(np.isclose(az_rot, azimuth))[0][0]


if vi <= 12:
    rlim = 6
    zlim = 2
    xl = [-6,0]
    yl = [-3,3]
else:
    rlim = 8
    zlim = 2.5
    xl = [-3,3]
    yl = [0,6]




# dbz/vel PPIs, dbz/vel RHIs, dbz/vel cross sections, vort PPIs, vort RHIs, vort cross sections
plot_flag = [1,1,1,1,1,1]

figsave = False



if plot_flag[0]:
    # xl = [-7, 3] # was [-rlim/2, 0]. [-6,0]/[-3,0] until vi 12 | [0,3] for vi 17
    # yl = [-3, 7] # was [0, rlim/2]. [-3,3]/[0,3] until vi 12 | [4,7] for vi 17
    
    fig,(ax1,ax2) = plt.subplots(1,2,figsize=(10,4), sharex=True, sharey=True, subplot_kw=dict(box_aspect=1), layout='constrained')
    
    plot_cfill(vol[vi]['xx'][eli,:,:], vol[vi]['yy'][eli,:,:], vol[vi]['dbz'][eli,:,:], 'dbz', ax1, datalims=[0,70], xlims=xl, ylims=yl)
    ax1.scatter(x_rot, y_rot, s=30, c='k', marker='.')
    ax1.set_title(f"{filetime} UTC {vol[vi]['elev'][eli].round(1)}\N{DEGREE SIGN} reflectivity", fontsize=14)
    ax1.set_xlabel('E-W distance from radar (km)', fontsize=12)
    ax1.set_ylabel('N-S distance from radar (km)', fontsize=12)
    # s1 = ax1.scatter(P1['xx'][i1], P1['yy'][i1], s=50, c=P1['Utube'][i1], cmap=cmaps['temp']['cm'],
    #                 vmin=datalims[0], vmax=datalims[1], marker='s', edgecolors='k')
    # s2 = ax1.scatter(P2['xx'][i2], P2['yy'][i2], s=50, c=P2['Utube'][i2], cmap=cmaps['temp']['cm'],
    #                 vmin=datalims[0], vmax=datalims[1], marker='^', edgecolors='k')
    # ax1.scatter(0, 0, s=50, c='k')
    # ax1.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')
    # plt.colorbar(s1,label='Sfc temperature (C)')
    # plt.legend([s1,s2], ['Probe 1', 'Probe 2'], loc='upper right')
    ax1.plot(vol[vi]['xx'][eli,azi,:], vol[vi]['yy'][eli,azi,:], '--k', linewidth=1.25)
    
    plot_cfill(vol[vi]['xx'][eli,:,:], vol[vi]['yy'][eli,:,:], vol[vi]['vel'][eli,:,:], 'vel', ax2, datalims=[-va,va], xlims=xl, ylims=yl)
    ax2.scatter(x_rot, y_rot, s=30, c='k', marker='.')
    ax2.set_title(f"{filetime} UTC {vol[vi]['elev'][eli].round(1)}\N{DEGREE SIGN} radial velocity", fontsize=14)
    ax2.set_xlabel('E-W distance from radar (km)', fontsize=12)
    ax2.set_ylabel('N-S distance from radar (km)', fontsize=12)
    # s1 = ax2.scatter(P1['xx'][i1], P1['yy'][i1], s=50, c=P1['Utube'][i1], cmap=cmaps['temp']['cm'],
    #                 vmin=datalims[0], vmax=datalims[1], marker='s', edgecolors='k')
    # s2 = ax2.scatter(P2['xx'][i2], P2['yy'][i2], s=50, c=P2['Utube'][i2], cmap=cmaps['temp']['cm'],
    #                 vmin=datalims[0], vmax=datalims[1], marker='^', edgecolors='k')
    # ax2.scatter(0, 0, s=50, c='k')
    # ax2.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')
    # plt.colorbar(s1,label='Sfc temperature (C)')
    # plt.legend([s1,s2], ['Probe 1', 'Probe 2'], loc='upper right')
    ax2.plot(vol[vi]['xx'][eli,azi,:], vol[vi]['yy'][eli,azi,:], '--k', linewidth=1.25)
    
    # plt.suptitle(f"{filetime} UTC, azimuth = {azimuth}\N{DEGREE SIGN}", fontsize=14)
    # plt.suptitle(f"{filetime} UTC", fontsize=14)
    if figsave:
        plt.savefig(ip+f"vol{vi}_{filetime}_PPI_vortex{vortex_num}.png", dpi=300)
    


if plot_flag[1]:
    fig,((ax1),(ax2)) = plt.subplots(2,1,figsize=(8,6), sharex=True, sharey=True, layout='constrained')
    
    plot_cfill(rr, vol[vi]['zz'][:,azi,:], vol[vi]['dbz'][:,azi,:], 'dbz', ax1, datalims=[0,70], xlims=[0,rlim], ylims=[0,zlim])
    ax1.scatter(r_rot[irot], z_rot[irot], s=10, c='k', marker='x')
    ax1.set_ylabel('Height ARL (km)', fontsize=14)
    ax1.set_title(f"{filetime} UTC (Azimuth = {azimuth}\N{DEGREE SIGN})\n Reflectivity", fontsize=14, fontweight='bold')
    ax1.invert_xaxis()
    ax1_ppi = inset_axes(ax1, '20%', '52%', loc=1)
    c = plot_cfill(vol[vi]['xx'][eli,:,:], vol[vi]['yy'][eli,:,:], vol[vi]['dbz'][eli,:,:], 'dbz', ax1_ppi, datalims=[0,70], xlims=xl, ylims=yl, cbar=False)
    ax1_ppi.plot(vol[vi]['xx'][eli,azi,:], vol[vi]['yy'][eli,azi,:], '--k', linewidth=1.25)
    ax1_ppi.set_xticks([])
    ax1_ppi.set_yticks([])
    
    
    plot_cfill(rr, vol[vi]['zz'][:,azi,:], vol[vi]['vel'][:,azi,:], 'vel', ax2, datalims=[-va,va], xlims=[0,rlim], ylims=[0,zlim])
    ax2.scatter(r_rot[irot], z_rot[irot], s=10, c='k', marker='x')
    ax2.set_xlabel('Range from radar (km)', fontsize=14)
    ax2.set_ylabel('Height ARL (km)', fontsize=14)
    ax2.set_title(f"Radial velocity", fontsize=14, fontweight='bold')
    ax2.invert_xaxis()
    ax2_ppi = inset_axes(ax2, '20%', '52%', loc=1)
    c = plot_cfill(vol[vi]['xx'][eli,:,:], vol[vi]['yy'][eli,:,:], vol[vi]['vel'][eli,:,:], 'vel', ax2_ppi, datalims=[-va,va], xlims=xl, ylims=yl, cbar=False)
    ax2_ppi.plot(vol[vi]['xx'][eli,azi,:], vol[vi]['yy'][eli,azi,:], '--k', linewidth=1.25)
    ax2_ppi.set_xticks([])
    ax2_ppi.set_yticks([])
    
    # plt.suptitle(f"{filetime} UTC, azimuth = {azimuth}\N{DEGREE SIGN}")
    if figsave:
        plt.savefig(ip+f"vol{vi}_{filetime}_az{azimuth}_RHI_vortex{vortex_num}.png", dpi=300, bbox_inches='tight')

#% rest of plots

if plot_flag[2]:
    ir = np.where(np.isclose(np.mean(r_rot), rr[eli,:], rtol=0.01))[0][0]
    ir1 = np.where(np.isclose(r_rot[0], rr[eli,:], rtol=0.01))[0][0]
    ir2 = np.where(np.isclose(r_rot[-1], rr[eli,:], rtol=0.01))[0][0]
    ia1 = np.where(np.isclose(az_rot[0], vol[vi]['az'][eli,:], atol=0.1))[0][0]
    ia2 = np.where(np.isclose(az_rot[-1], vol[vi]['az'][eli,:], atol=0.1))[0][0]
    
    if ia2 > ia1:
        raz = wrf.xy(vol[vi]['dbz'], start_point=(ir,ia1), end_point=(ir,ia2))
        dbz_cs = wrf.interp2dxy(vol[vi]['dbz'], raz)
        dbz_cross = dbz_cs.data
        
        raz = wrf.xy(vol[vi]['vel'], start_point=(ir,ia1), end_point=(ir,ia2))
        vel_cs = wrf.interp2dxy(vol[vi]['vel'], raz)
        vel_cross = vel_cs.data
        
        raz = wrf.xy(vol[vi]['zdr'], start_point=(ir,ia1), end_point=(ir,ia2))
        zdr_cs = wrf.interp2dxy(vol[vi]['zdr'], raz)
        zdr_cross = zdr_cs.data
        
        zz_rot = np.mean(vol[vi]['zz'][:,ia1,slice(ir1,ir2+1)], axis=1)
        az_rot2 = az_rot
        if az_rot.shape[0] < dbz_cross.shape[1]:
            az_rot2 = np.linspace(az_rot[0]-0.5, az_rot[-1]+0.5, len(az_rot)+1)
    
    else:
        raz1 = wrf.xy(vol[vi]['dbz'], start_point=(ir,ia1), end_point=(ir,359))
        dbz_cs1 = wrf.interp2dxy(vol[vi]['dbz'], raz1)
        raz2 = wrf.xy(vol[vi]['dbz'], start_point=(ir,0), end_point=(ir,ia2))
        dbz_cs2 = wrf.interp2dxy(vol[vi]['dbz'], raz2)
        dbz_cross = np.append(dbz_cs1.data, dbz_cs2.data, axis=1)
        
        raz1 = wrf.xy(vol[vi]['vel'], start_point=(ir,ia1), end_point=(ir,359))
        vel_cs1 = wrf.interp2dxy(vol[vi]['vel'], raz1)
        raz2 = wrf.xy(vol[vi]['vel'], start_point=(ir,0), end_point=(ir,ia2))
        vel_cs2 = wrf.interp2dxy(vol[vi]['vel'], raz2)
        vel_cross = np.append(vel_cs1.data, vel_cs2.data, axis=1)
        
        raz1 = wrf.xy(vol[vi]['zdr'], start_point=(ir,ia1), end_point=(ir,359))
        zdr_cs1 = wrf.interp2dxy(vol[vi]['zdr'], raz1)
        raz2 = wrf.xy(vol[vi]['zdr'], start_point=(ir,0), end_point=(ir,ia2))
        zdr_cs2 = wrf.interp2dxy(vol[vi]['zdr'], raz2)
        zdr_cross = np.append(zdr_cs1.data, zdr_cs2.data, axis=1)
        
        zz_rot = np.mean(vol[vi]['zz'][:,ia1,slice(ir1,ir2+1)], axis=1)
        az_rot2 = np.array([az_rot[i]+360 if az_rot[i]<az_rot[0] else az_rot[i] for i in range(len(az_rot))])
        if az_rot2.shape[0] < dbz_cross.shape[1]:
            az_rot2 = np.linspace(az_rot2[0]-0.5, az_rot2[-1]+0.5, len(az_rot2)+1)
    
    
    zl = [0, np.max(zz_rot)]
    # if vi == 7:
    #     zl = [0, 0.7]
    # elif vi == 8:
    #     zl = [0, 0.75]
    # elif vi == 9:
    #     zl = [0, 0.6]
    # elif vi == 10:
    #     zl = [0, 0.6]
    # elif vi == 12:
    #     zl = [0, 0.6]
    # elif vi == 13:
    #     zl = [0, 1.25]
    # elif vi == 14:
    #     zl = [0, 1.4]
    # elif vi == 15:
    #     zl = [0, 1.5]
    # elif vi == 16:
    #     zl = [0, 1.6]
    # elif vi == 17:
    #     zl = [0, 2]
    
    
    
    fig,(ax1,ax2) = plt.subplots(2,1, figsize=(9,7), sharex=True, layout='constrained')
    
    plot_cfill(az_rot2, zz_rot, dbz_cross, 'dbz', ax1, datalims=[0,70], xlims=[az_rot2[0],az_rot2[-1]], ylims=zl)
    ax1.set_title(f"{filetime} UTC Reflectivity cross-section", fontsize=14, fontweight='bold')
    ax1.set_ylabel('Height ARL (km)', fontsize=12)
    ax1_ppi = inset_axes(ax1, '21%', '52%', loc=2) # 21, 52 OR 23, 56
    c = plot_cfill(vol[vi]['xx'][eli,:,:], vol[vi]['yy'][eli,:,:], vol[vi]['dbz'][eli,:,:], 'dbz', ax1_ppi, datalims=[0,70], xlims=xl, ylims=yl, cbar=False)
    ax1_ppi.scatter(x_rot, y_rot, s=4, c='k', marker='.')
    # ax1_ppi.plot([vol[vi]['xx'][eli,ia1,ir1], vol[vi]['xx'][eli,ia2,ir2]], [vol[vi]['yy'][eli,ia1,ir1], vol[vi]['yy'][eli,ia2,ir2]], '-k', linewidth=1.25)
    ax1_ppi.plot([0,x_rot[0]], [0,y_rot[0]], '--k', linewidth=1)
    ax1_ppi.plot([0,x_rot[-1]], [0,y_rot[-1]], '--k', linewidth=1)
    ax1_ppi.set_xticks([])
    ax1_ppi.set_yticks([])
    
    plot_cfill(az_rot2, zz_rot, vel_cross, 'vel', ax2, datalims=[-va,va], xlims=[az_rot2[0],az_rot2[-1]], ylims=zl)
    ax2.set_title(f"{filetime} UTC Radial velocity cross-section", fontsize=14, fontweight='bold')
    ax2.set_xlabel('Azimuth (deg)', fontsize=12)
    ax2.set_ylabel('Height ARL (km)', fontsize=12)
    ax2.set_xticklabels([f"{az_rot[i]}" for i in range(len(az_rot))])
    ax2_ppi = inset_axes(ax2, '21%', '52%', loc=2)
    c = plot_cfill(vol[vi]['xx'][eli,:,:], vol[vi]['yy'][eli,:,:], vol[vi]['vel'][eli,:,:], 'vel', ax2_ppi, datalims=[-va,va], xlims=xl, ylims=yl, cbar=False)
    ax2_ppi.scatter(x_rot, y_rot, s=4, c='k', marker='.')
    ax2_ppi.plot([0,x_rot[0]], [0,y_rot[0]], '--k', linewidth=1)
    ax2_ppi.plot([0,x_rot[-1]], [0,y_rot[-1]], '--k', linewidth=1)
    ax2_ppi.set_xticks([])
    ax2_ppi.set_yticks([])
    
    if figsave:
        plt.savefig(ip+f"vol{vi}_{filetime}_cross-section_vortex{vortex_num}.png", dpi=300, bbox_inches='tight')
    




# Pseudovorticity plots
if plot_flag[3]:
    # Column max PPIs
    fig,(ax1,ax2) = plt.subplots(1,2,figsize=(10,4), sharex=True, sharey=True, subplot_kw=dict(box_aspect=1), layout='constrained')
    
    plot_cfill(vol[vi]['xx'][eli,:,:], vol[vi]['yy'][eli,:,:], np.max(vol[vi]['zvort'][2:,:,:],axis=0), 'vort', ax1, datalims=[0,vort_lim], xlims=xl, ylims=yl)
    ax1.scatter(x_rot, y_rot, s=20, facecolor='w', edgecolor='k', marker='o', linewidth=0.5)
    ax1.set_title(f"{filetime} UTC Max vertical pseudovorticity", fontsize=14)
    ax1.set_xlabel('E-W distance from radar (km)', fontsize=12)
    ax1.set_ylabel('N-S distance from radar (km)', fontsize=12)
    ax1.plot(vol[vi]['xx'][eli,azi,:], vol[vi]['yy'][eli,azi,:], '--k', linewidth=1.25)
    
    plot_cfill(vol[vi]['xx'][eli,:,:], vol[vi]['yy'][eli,:,:], np.max(np.abs(vol[vi]['hvort'][2:,:,:]),axis=0), 'vort', ax2, datalims=[0,vort_lim], xlims=xl, ylims=yl)
    ax2.scatter(x_rot, y_rot, s=20, facecolor='w', edgecolor='k', marker='o', linewidth=0.5)
    ax2.set_title(f"{filetime} UTC Max horizontal pseudovorticity", fontsize=14)
    ax2.set_xlabel('E-W distance from radar (km)', fontsize=12)
    ax2.set_ylabel('N-S distance from radar (km)', fontsize=12)
    ax2.plot(vol[vi]['xx'][eli,azi,:], vol[vi]['yy'][eli,azi,:], '--k', linewidth=1.25)
    
    # plot_cfill(vol[vi]['xx'][eli,:,:], vol[vi]['yy'][eli,:,:], np.min(vol[vi]['div'], axis=0), 'div', ax1, 
    #            datalims=[-0.16,0], xlims=xl, ylims=yl, cmap='pyart_ChaseSpectral_r')
    # ax1.scatter(x_rot, y_rot, s=30, c='k', marker='.')
    # ax1.set_title(f"{filetime} UTC {vol[vi]['elev'][eli].round(1)}\N{DEGREE SIGN} Pseudodivergence", fontsize=14)
    # ax1.set_xlabel('E-W distance from radar (km)', fontsize=12)
    # ax1.set_ylabel('N-S distance from radar (km)', fontsize=12)
    
    if figsave:
        plt.savefig(ip+f"vol{vi}_{filetime}_PPI_vort_vortex{vortex_num}.png", dpi=300)
    
    
    
if plot_flag[4]:
    # Reconstructed PPIs
    # Positive hvort is LEFT OF RADIAL, negative hvort is RIGHT OF RADIAL
    fig,((ax1),(ax2)) = plt.subplots(2,1,figsize=(8,6), sharex=True, sharey=True, layout='constrained')
    
    plot_cfill(rr, vol[vi]['zz'][:,azi,:], vol[vi]['zvort'][:,azi,:], 'vort', ax1, datalims=[-vort_lim,vort_lim], xlims=[0,rlim], ylims=[0,zlim])
    ax1.scatter(r_rot[irot], z_rot[irot], s=10, c='k', marker='x')
    ax1.set_ylabel('Height ARL (km)', fontsize=14)
    ax1.set_title(f"{filetime} UTC (Azimuth = {azimuth}\N{DEGREE SIGN})\n Vertical pseudovorticity", fontsize=14, fontweight='bold')
    ax1.invert_xaxis()
    ax1_ppi = inset_axes(ax1, '20%', '52%', loc=1)
    c = plot_cfill(vol[vi]['xx'][eli,:,:], vol[vi]['yy'][eli,:,:], np.max(vol[vi]['zvort'][2:,:,:],axis=0), 'vort', ax1_ppi, datalims=[0,vort_lim], xlims=xl, ylims=yl, cbar=False)
    ax1_ppi.plot(vol[vi]['xx'][eli,azi,:], vol[vi]['yy'][eli,azi,:], '--k', linewidth=1.25)
    ax1_ppi.set_xticks([])
    ax1_ppi.set_yticks([])
    
    plot_cfill(rr, vol[vi]['zz'][:,azi,:], vol[vi]['hvort'][:,azi,:], 'vort', ax2, datalims=[-vort_lim,vort_lim], xlims=[0,rlim], ylims=[0,zlim])
    ax2.scatter(r_rot[irot], z_rot[irot], s=10, c='k', marker='x')
    ax2.set_xlabel('Range from radar (km)', fontsize=14)
    ax2.set_ylabel('Height ARL (km)', fontsize=14)
    ax2.set_title(f"Horizontal pseudovorticity", fontsize=14, fontweight='bold')
    ax2.invert_xaxis()
    ax2_ppi = inset_axes(ax2, '20%', '52%', loc=1)
    c = plot_cfill(vol[vi]['xx'][eli,:,:], vol[vi]['yy'][eli,:,:], np.max(np.abs(vol[vi]['hvort'][2:,:,:]),axis=0), 'vort', ax2_ppi, datalims=[0,vort_lim], xlims=xl, ylims=yl, cbar=False)
    ax2_ppi.plot(vol[vi]['xx'][eli,azi,:], vol[vi]['yy'][eli,azi,:], '--k', linewidth=1.25)
    ax2_ppi.set_xticks([])
    ax2_ppi.set_yticks([])
    
    # plot_cfill(rr, vol[vi]['zz'][:,azi,:], vol[vi]['div'][:,azi,:], 'div', ax1, datalims=[-0.1,0.1], xlims=[0,rlim], ylims=[0,zlim])
    # # ax1.scatter(r_rot[irot], z_rot[irot], s=10, c='k', marker='x')
    # ax1.set_ylabel('Height ARL (km)', fontsize=14)
    # ax1.set_title(f"{filetime} UTC (Azimuth = {azimuth}\N{DEGREE SIGN})\n Pseudodivergence", fontsize=14, fontweight='bold')
    # ax1.invert_xaxis()
    # ax1_ppi = inset_axes(ax1, '20%', '52%', loc=1)
    # c = plot_cfill(vol[vi]['xx'][eli,:,:], vol[vi]['yy'][eli,:,:], vol[vi]['dbz'][eli,:,:], 'dbz', ax1_ppi, datalims=[0,70], xlims=xl, ylims=yl, cbar=False)
    # ax1_ppi.plot(vol[vi]['xx'][eli,azi,:], vol[vi]['yy'][eli,azi,:], '--k', linewidth=1.25)
    # ax1_ppi.set_xticks([])
    # ax1_ppi.set_yticks([])
    
    if figsave:
        plt.savefig(ip+f"vol{vi}_{filetime}_RHI_az{azimuth}_vort_vortex{vortex_num}.png", dpi=300, bbox_inches='tight')
    
    
if plot_flag[5]:
    # Along-radial cross sections
    ir = np.where(np.isclose(np.mean(r_rot), rr[eli,:], rtol=0.01))[0][0]
    ir1 = np.where(np.isclose(r_rot[0], rr[eli,:], rtol=0.01))[0][0]
    ir2 = np.where(np.isclose(r_rot[-1], rr[eli,:], rtol=0.01))[0][0]
    ia1 = np.where(np.isclose(az_rot[0], vol[vi]['az'][eli,:], atol=0.1))[0][0]
    ia2 = np.where(np.isclose(az_rot[-1], vol[vi]['az'][eli,:], atol=0.1))[0][0]
    
    if ia2 > ia1:
        raz = wrf.xy(vol[vi]['hvort'], start_point=(ir,ia1), end_point=(ir,ia2))
        hvort_cs = wrf.interp2dxy(vol[vi]['hvort'], raz)
        hvort_cross = hvort_cs.data
        
        raz = wrf.xy(vol[vi]['zvort'], start_point=(ir,ia1), end_point=(ir,ia2))
        zvort_cs = wrf.interp2dxy(vol[vi]['zvort'], raz)
        zvort_cross = zvort_cs.data
        
        raz = wrf.xy(vol[vi]['div'], start_point=(ir,ia1), end_point=(ir,ia2))
        div_cs = wrf.interp2dxy(vol[vi]['div'], raz)
        div_cross = div_cs.data
        
        zz_rot = np.mean(vol[vi]['zz'][:,ia1,slice(ir1,ir2+1)], axis=1)
        if az_rot.shape[0] < dbz_cross.shape[1]:
            az_rot2 = np.linspace(az_rot[0]-0.5, az_rot[-1]+0.5, len(az_rot)+1)
    else:
        raz1 = wrf.xy(vol[vi]['hvort'], start_point=(ir,ia1), end_point=(ir,359))
        hvort_cs1 = wrf.interp2dxy(vol[vi]['hvort'], raz1)
        raz2 = wrf.xy(vol[vi]['hvort'], start_point=(ir,0), end_point=(ir,ia2))
        hvort_cs2 = wrf.interp2dxy(vol[vi]['hvort'], raz2)
        hvort_cross = np.append(hvort_cs1.data, hvort_cs2.data, axis=1)
        
        raz1 = wrf.xy(vol[vi]['zvort'], start_point=(ir,ia1), end_point=(ir,359))
        zvort_cs1 = wrf.interp2dxy(vol[vi]['zvort'], raz1)
        raz2 = wrf.xy(vol[vi]['zvort'], start_point=(ir,0), end_point=(ir,ia2))
        zvort_cs2 = wrf.interp2dxy(vol[vi]['zvort'], raz2)
        zvort_cross = np.append(zvort_cs1.data, zvort_cs2.data, axis=1)
        
        raz1 = wrf.xy(vol[vi]['div'], start_point=(ir,ia1), end_point=(ir,359))
        div_cs1 = wrf.interp2dxy(vol[vi]['div'], raz1)
        raz2 = wrf.xy(vol[vi]['div'], start_point=(ir,0), end_point=(ir,ia2))
        div_cs2 = wrf.interp2dxy(vol[vi]['div'], raz2)
        div_cross = np.append(div_cs1.data, div_cs2.data, axis=1)
        
        zz_rot = np.mean(vol[vi]['zz'][:,ia1,slice(ir1,ir2+1)], axis=1)
        az_rot2 = np.array([az_rot[i]+360 if az_rot[i]<az_rot[0] else az_rot[i] for i in range(len(az_rot))])
        if az_rot2.shape[0] < dbz_cross.shape[1]:
            az_rot2 = np.linspace(az_rot2[0]-0.5, az_rot2[-1]+0.5, len(az_rot2)+1)
    
    zl = [0, np.max(zz_rot)]
    
    # if vi == 7:
    #     zl = [0, 0.7]
    # elif vi == 8:
    #     zl = [0, 0.75]
    # elif vi == 9:
    #     zl = [0, 0.6]
    # elif vi == 10:
    #     zl = [0, 0.6]
    # elif vi == 12:
    #     zl = [0, 0.6]
    # elif vi == 13:
    #     zl = [0, 1.25]
    # elif vi == 14:
    #     zl = [0, 1.4]
    # elif vi == 15:
    #     zl = [0, 1.5]
    # elif vi == 16:
    #     zl = [0, 1.6]
    # elif vi == 17:
    #     zl = [0, 2]
    
    
    
    fig,(ax1,ax2) = plt.subplots(2,1, figsize=(9,7), sharex=True, layout='constrained')
    
    plot_cfill(az_rot2, zz_rot, zvort_cross, 'vort', ax1, datalims=[-vort_lim,vort_lim], xlims=[az_rot2[0],az_rot2[-1]], ylims=zl)
    ax1.set_title(f"{filetime} UTC Vertical pseudovorticity cross-section", fontsize=14, fontweight='bold')
    ax1.set_ylabel('Height ARL (km)', fontsize=12)
    ax1_ppi = inset_axes(ax1, '21%', '52%', loc=2) # 21, 52 OR 23, 56
    c = plot_cfill(vol[vi]['xx'][eli,:,:], vol[vi]['yy'][eli,:,:], np.max(vol[vi]['zvort'][2:,:,:],axis=0), 'vort', ax1_ppi, datalims=[0,vort_lim], xlims=xl, ylims=yl, cbar=False)
    ax1_ppi.scatter(x_rot, y_rot, s=4, c='w', marker='.')
    # ax1_ppi.plot([vol[vi]['xx'][eli,ia1,ir1], vol[vi]['xx'][eli,ia2,ir2]], [vol[vi]['yy'][eli,ia1,ir1], vol[vi]['yy'][eli,ia2,ir2]], '-k', linewidth=1.25)
    ax1_ppi.plot([0,x_rot[0]], [0,y_rot[0]], '--w', linewidth=1)
    ax1_ppi.plot([0,x_rot[-1]], [0,y_rot[-1]], '--w', linewidth=1)
    ax1_ppi.set_xticks([])
    ax1_ppi.set_yticks([])
    
    plot_cfill(az_rot2, zz_rot, hvort_cross, 'vort', ax2, datalims=[-vort_lim,vort_lim], xlims=[az_rot2[0],az_rot2[-1]], ylims=zl)
    ax2.set_title(f"{filetime} UTC Horizontal pseudovorticity cross-section", fontsize=14, fontweight='bold')
    ax2.set_xlabel('Azimuth (deg)', fontsize=12)
    ax2.set_ylabel('Height ARL (km)', fontsize=12)
    ax2.set_xticklabels([f"{az_rot[i]}" for i in range(len(az_rot))])
    ax2_ppi = inset_axes(ax2, '21%', '52%', loc=2)
    c = plot_cfill(vol[vi]['xx'][eli,:,:], vol[vi]['yy'][eli,:,:], np.max(np.abs(vol[vi]['hvort'][2:,:,:]),axis=0), 'vort', ax2_ppi, datalims=[0,vort_lim], xlims=xl, ylims=yl, cbar=False)
    ax2_ppi.scatter(x_rot, y_rot, s=4, c='w', marker='.')
    ax2_ppi.plot([0,x_rot[0]], [0,y_rot[0]], '--w', linewidth=1)
    ax2_ppi.plot([0,x_rot[-1]], [0,y_rot[-1]], '--w', linewidth=1)
    ax2_ppi.set_xticks([])
    ax2_ppi.set_yticks([])
    
    # plot_cfill(az_rot, zz_rot, div_cross, 'div', ax1, datalims=[-0.1,0.1], xlims=[az_rot[0],az_rot[-1]], ylims=zl)
    # # plot_contourf(az_rot, zz_rot, div_cross, 'vort', ax1, levels=np.linspace(-0.16,0.16,33), datalims=[-0.16,0.16], xlims=[az_rot[0],az_rot[-1]], ylims=yl)
    # ax1.set_title(f"{filetime} UTC Pseudodivergence cross-section", fontsize=14, fontweight='bold')
    # ax1.set_ylabel('Height ARL (km)', fontsize=12)
    # ax1_ppi = inset_axes(ax1, '21%', '52%', loc=2) # 21, 52 OR 23, 56
    # c = plot_cfill(vol[vi]['xx'][eli,:,:], vol[vi]['yy'][eli,:,:], vol[vi]['dbz'][eli,:,:], 'dbz', ax1_ppi, datalims=[0,70], xlims=xl, ylims=yl, cbar=False)
    # ax1_ppi.scatter(x_rot, y_rot, s=4, c='k', marker='.')
    # # ax1_ppi.plot([vol[vi]['xx'][eli,ia1,ir1], vol[vi]['xx'][eli,ia2,ir2]], [vol[vi]['yy'][eli,ia1,ir1], vol[vi]['yy'][eli,ia2,ir2]], '-k', linewidth=1.25)
    # ax1_ppi.plot([0,x_rot[0]], [0,y_rot[0]], '--k', linewidth=1)
    # ax1_ppi.plot([0,x_rot[-1]], [0,y_rot[-1]], '--k', linewidth=1)
    # ax1_ppi.set_xticks([])
    # ax1_ppi.set_yticks([])
    
    if figsave:
        plt.savefig(ip+f"vol{vi}_{filetime}_rotor-cross-section_vort_vortex{vortex_num}.png", dpi=300, bbox_inches='tight')
    
    
#%% mobile mesonet time series

figsave = False

ip = '/Users/morgan.schneider/Documents/perils2023/iop2/figs/'

if True:
    import matplotlib.dates as mdates
    from datetime import datetime,timedelta
    
    i1_p1 = np.where(P1['time'] >= int(vol[0]['scan_time'][0]))[0][0]
    i2_p1 = np.where(P1['time'] >= int(vol[-1]['scan_time'][-1]))[0][1]
    i1_p2 = np.where(P2['time'] >= int(vol[0]['scan_time'][0]))[0][0]
    i2_p2 = np.where(P2['time'] >= int(vol[-1]['scan_time'][-1]))[0][1]
    i1 = slice(i1_p1,i2_p1)
    i2 = slice(i1_p2,i2_p2)
    
    
    P1_h = [float(str(P1['time'][i1][i])[0]) for i in range(len(P1['time'][i1]))]
    P1_m = [float(str(P1['time'][i1][i])[1:3]) for i in range(len(P1['time'][i1]))]
    P1_s = [float(str(P1['time'][i1][i])[3:]) for i in range(len(P1['time'][i1]))]
    P2_h = [float(str(P2['time'][i2][i])[0]) for i in range(len(P2['time'][i2]))]
    P2_m = [float(str(P2['time'][i2][i])[1:3]) for i in range(len(P2['time'][i2]))]
    P2_s = [float(str(P2['time'][i2][i])[3:]) for i in range(len(P2['time'][i2]))]
    
    dt = datetime(2023,3,3)
    
    P1_datetimes = [dt + timedelta(hours=P1_h[i], minutes=P1_m[i], seconds=P1_s[i]) for i in range(len(P1['time'][i1]))]
    P2_datetimes = [dt + timedelta(hours=P2_h[i], minutes=P2_m[i], seconds=P2_s[i]) for i in range(len(P2['time'][i2]))]
    
    P1_times = mdates.date2num(P1_datetimes)
    P2_times = mdates.date2num(P2_datetimes)
    
    fig,ax = plt.subplots(4,1,figsize=(8,9), sharex=True, layout='constrained')
    
    ax[0].plot_date(P1_times, P1['Utube'][i1], 'b', linewidth=1.5)
    ax[0].plot_date(P2_times, P2['Utube'][i2], 'r', linewidth=1.5)
    ax[0].set_ylabel('Temperature (C)', fontsize=14)
    ax[0].set_title('Mobile mesonet temperature', fontsize=16)
    ax[0].set_ylim([18.5,20.5])
    ax[0].legend(['Probe 1','Probe 2'], loc='upper right', fontsize=14)
    ax[0].tick_params(axis='y', which='major', labelsize=12)
    
    ax[1].plot_date(P1_times, P1['Dewpt'][i1], 'b', linewidth=1.5)
    ax[1].plot_date(P2_times, P2['Dewpt'][i2], 'r', linewidth=1.5)
    ax[1].set_ylabel('Dewpoint (C)', fontsize=14)
    ax[1].set_title('Mobile mesonet dewpoint', fontsize=16)
    ax[1].set_ylim([18,19])
    ax[1].tick_params(axis='y', which='major', labelsize=12)
    
    ax[2].plot_date(P1_times, P1['wspd_10s'][i1], 'b', linewidth=1.5)
    ax[2].plot_date(P2_times, P2['wspd_10s'][i2], 'r', linewidth=1.5)
    ax[2].set_ylabel('Wind speed (m/s)', fontsize=14)
    ax[2].set_title('Mobile mesonet 10-s wind speed', fontsize=16)
    ax[2].set_ylim([0,15])
    ax[2].tick_params(axis='y', which='major', labelsize=12)
    
    ax[3].plot_date(P1_times, P1['Pres_3s'][i1], 'b', linewidth=1.5)
    ax[3].plot_date(P2_times, P2['Pres_3s'][i2], 'r', linewidth=1.5)
    ax[3].set_xlabel('Time (UTC)', fontsize=14)
    ax[3].set_ylabel('Pressure (mb)', fontsize=14)
    ax[3].set_title('Mobile mesonet 3-s pressure', fontsize=16)
    ax[3].set_ylim([984,988])
    
    ax[3].set_xlim(datetime(2023,3,3,8,16,0), datetime(2023,3,3,8,27,0))
    ax[3].xaxis.set_major_locator(mdates.MinuteLocator())
    ax[3].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    ax[3].tick_params(which='major', labelsize=12)
    
    fig.autofmt_xdate()
    
    if figsave:
        plt.savefig(ip+f"MM_sfc_timeseries.png", dpi=300, bbox_inches='tight')
    
    
if True:
    fig,ax = plt.subplots(3,1,figsize=(8,8), sharex=True, layout='constrained')
    
    ax[0].plot_date(P1_times, P1['Theta'][i1], 'b', linewidth=1.5)
    ax[0].plot_date(P2_times, P2['Theta'][i2], 'r', linewidth=1.5)
    ax[0].set_ylabel('Theta (K)', fontsize=14)
    ax[0].set_title('Mobile mesonet Theta', fontsize=16)
    ax[0].set_ylim([293,295])
    ax[0].legend(['Probe 1','Probe 2'], loc='upper right', fontsize=14)
    ax[0].tick_params(axis='y', which='major', labelsize=12)
    
    ax[1].plot_date(P1_times, P1['ThetaV'][i1], 'b', linewidth=1.5)
    ax[1].plot_date(P2_times, P2['ThetaV'][i2], 'r', linewidth=1.5)
    ax[1].set_ylabel('ThetaV (K)', fontsize=14)
    ax[1].set_title('Mobile mesonet ThetaV', fontsize=16)
    ax[1].set_ylim([295.5,297.5])
    ax[1].tick_params(axis='y', which='major', labelsize=12)
    
    ax[2].plot_date(P1_times, P1['ThetaE'][i1], 'b', linewidth=1.5)
    ax[2].plot_date(P2_times, P2['ThetaE'][i2], 'r', linewidth=1.5)
    ax[2].set_xlabel('Time (UTC)', fontsize=14)
    ax[2].set_ylabel('ThetaE (K)', fontsize=14)
    ax[2].set_title('Mobile mesonet ThetaE', fontsize=16)
    ax[2].set_ylim([293,295])
    
    ax[2].set_xlim(datetime(2023,3,3,8,16,0), datetime(2023,3,3,8,27,0))
    ax[2].xaxis.set_major_locator(mdates.MinuteLocator())
    ax[2].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    ax[2].tick_params(which='major', labelsize=12)
    
    fig.autofmt_xdate()
    
    if figsave:
        plt.savefig(ip+f"MM_theta_timeseries.png", dpi=300, bbox_inches='tight')

if True:
    qix = 20
    fig,ax = plt.subplots(2,1, figsize=(8,6), sharex=True, layout='constrained')
    
    ax[0].plot_date(P1_times, P1['Utube'][i1], 'dodgerblue', linewidth=1)
    ax[0].plot_date(P2_times, P2['Utube'][i2], 'hotpink', linewidth=1)
    ax[0].barbs(P1_times[::qix], P1['Utube'][i1][::qix], P1['uwind'][i1][::qix], P1['vwind'][i1][::qix], barbcolor='b', length=6)
    ax[0].barbs(P2_times[::qix], P2['Utube'][i2][::qix], P2['uwind'][i2][::qix], P2['vwind'][i2][::qix], barbcolor='r', length=6)
    # ax[0].set_xlabel('Time (UTC)', fontsize=14)
    ax[0].set_ylabel('Temperature (C)', fontsize=14)
    ax[0].set_ylim([18.5,20.5])
    ax[0].legend(['Probe 1','Probe 2'], loc='upper right', fontsize=14)
    ax[0].set_title('Mobile mesonet temperature', fontsize=16)
    
    ax[1].plot_date(P1_times, P1['wspd_10s'][i1], 'dodgerblue', linewidth=1)
    ax[1].plot_date(P2_times, P2['wspd_10s'][i2], 'hotpink', linewidth=1)
    ax[1].barbs(P1_times[::qix], P1['wspd_10s'][i1][::qix], P1['uwind'][i1][::qix], P1['vwind'][i1][::qix], barbcolor='b', length=6)
    ax[1].barbs(P2_times[::qix], P2['wspd_10s'][i2][::qix], P2['uwind'][i2][::qix], P2['vwind'][i2][::qix], barbcolor='r', length=6)
    ax[1].set_xlabel('Time (UTC)', fontsize=14)
    ax[1].set_ylabel('Wind speed (m/s)', fontsize=14)
    ax[1].set_ylim([0,15])
    ax[1].set_title('Mobile mesonet 10-s wind speed', fontsize=16)
    
    ax[1].set_xlim(datetime(2023,3,3,8,16,0), datetime(2023,3,3,8,27,0))
    ax[1].xaxis.set_major_locator(mdates.MinuteLocator())
    ax[1].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    ax[1].tick_params(which='major', labelsize=12)
    
    fig.autofmt_xdate()
    
    if figsave:
        plt.savefig(ip+f"MM_timeseries_withbarbs.png", dpi=300, bbox_inches='tight')




#%% Big raxpol overview figure like in AMS Radar poster

vi = [6, 9, 12, 16] # 0819, 0821, 0823, 0825
tmin = 19
tmax = 20.4


fig,ax = plt.subplots(4, 4, figsize=(13.1,12), sharex=True, sharey=True, subplot_kw=dict(box_aspect=1), layout='constrained')

for i in range(len(vi)):
    if vi[i] == 16:
        cb_flag = True
    else:
        cb_flag = False
    
    ip1 = [np.where(P1['time'] >= int(vol[vi[i]]['scan_time'][1]))[0][0],
          np.where(P1['time'] >= int(vol[vi[i]]['scan_time'][6]))[0][0]]
    ip2 = [np.where(P2['time'] >= int(vol[vi[i]]['scan_time'][1]))[0][0],
          np.where(P2['time'] >= int(vol[vi[i]]['scan_time'][6]))[0][0]]
    
    
    plot_cfill(vol[vi[i]]['xx'][1,:,:], vol[vi[i]]['yy'][1,:,:], np.ma.masked_array(vol[vi[i]]['dbz'][1,:,:], vol[vi[i]]['dbz'][1,:,:]<1),
               'dbz', ax[0,i], datalims=[0,70], xlims=[-5,5], ylims=[-2,8], cbar=cb_flag)
    b1 = ax[0,i].barbs(P1['xx'][ip1[0]], P1['yy'][ip1[0]], P1['uwind'][ip1[0]], P1['vwind'][ip1[0]], barbcolor='k', length=9)
    b2 = ax[0,i].barbs(P2['xx'][ip2[0]], P2['yy'][ip2[0]], P2['uwind'][ip2[0]], P2['vwind'][ip2[0]], barbcolor='k', length=9)
    s1 = ax[0,i].scatter(P1['xx'][ip1[0]], P1['yy'][ip1[0]], s=150, c=P1['Utube'][ip1[0]], cmap=cmaps['temp']['cm'],
                    vmin=tmin, vmax=tmax, marker='s', edgecolors='k')
    s2 = ax[0,i].scatter(P2['xx'][ip2[0]], P2['yy'][ip2[0]], s=175, c=P2['Utube'][ip2[0]], cmap=cmaps['temp']['cm'],
                    vmin=tmin, vmax=tmax, marker='^', edgecolors='k')
    ax[0,i].scatter(0, 0, s=75, c='k')
    ax[0,i].text(0, 0.4, 'RaXPol', fontsize=16, fontweight='bold')
    ax[0,i].legend([s1,s2], ['Probe 1', 'Probe 2'], loc='lower right')
    ax[0,i].set_title(f"{vol[vi[i]]['scan_time'][1]} UTC - {vol[vi[i]]['elev'][1]:.1f}\N{DEGREE SIGN}", fontsize=16, fontweight='bold')
    if cb_flag:
        plt.colorbar(s1, ax=ax[0,i], label='Sfc temperature (C)', extend='both')#, ticks=np.linspace(tmin,tmax,7))
    
    plot_cfill(vol[vi[i]]['xx'][1,:,:], vol[vi[i]]['yy'][1,:,:], np.ma.masked_array(vol[vi[i]]['vel'][1,:,:], vol[vi[i]]['dbz'][1,:,:]<1),
               'vel', ax[1,i], datalims=[-va,va], xlims=[-5,5], ylims=[-2,8], cbar=cb_flag)
    b1 = ax[1,i].barbs(P1['xx'][ip1[0]], P1['yy'][ip1[0]], P1['uwind'][ip1[0]], P1['vwind'][ip1[0]], barbcolor='k', length=9)
    b2 = ax[1,i].barbs(P2['xx'][ip2[0]], P2['yy'][ip2[0]], P2['uwind'][ip2[0]], P2['vwind'][ip2[0]], barbcolor='k', length=9)
    s1 = ax[1,i].scatter(P1['xx'][ip1[0]], P1['yy'][ip1[0]], s=150, c=P1['Utube'][ip1[0]], cmap=cmaps['temp']['cm'],
                    vmin=tmin, vmax=tmax, marker='s', edgecolors='k')
    s2 = ax[1,i].scatter(P2['xx'][ip2[0]], P2['yy'][ip2[0]], s=175, c=P2['Utube'][ip2[0]], cmap=cmaps['temp']['cm'],
                    vmin=tmin, vmax=tmax, marker='^', edgecolors='k')
    ax[1,i].scatter(0, 0, s=75, c='k')
    ax[1,i].text(0, 0.4, 'RaXPol', fontsize=16, fontweight='bold')
    ax[1,i].legend([s1,s2], ['Probe 1', 'Probe 2'], loc='lower right')
    if cb_flag:
        plt.colorbar(s1, ax=ax[1,i], label='Sfc temperature (C)', extend='both')#, ticks=np.linspace(tmin,tmax,7))
    
    
    plot_cfill(vol[vi[i]]['xx'][6,:,:], vol[vi[i]]['yy'][6,:,:], np.ma.masked_array(vol[vi[i]]['dbz'][6,:,:], vol[vi[i]]['dbz'][6,:,:]<1),
               'dbz', ax[2,i], datalims=[0,70], xlims=[-5,5], ylims=[-2,8], cbar=cb_flag)
    b1 = ax[2,i].barbs(P1['xx'][ip1[1]], P1['yy'][ip1[1]], P1['uwind'][ip1[1]], P1['vwind'][ip1[1]], barbcolor='k', length=9)
    b2 = ax[2,i].barbs(P2['xx'][ip2[1]], P2['yy'][ip2[1]], P2['uwind'][ip2[1]], P2['vwind'][ip2[1]], barbcolor='k', length=9)
    s1 = ax[2,i].scatter(P1['xx'][ip1[1]], P1['yy'][ip1[1]], s=150, c=P1['Utube'][ip1[1]], cmap=cmaps['temp']['cm'],
                    vmin=tmin, vmax=tmax, marker='s', edgecolors='k')
    s2 = ax[2,i].scatter(P2['xx'][ip2[1]], P2['yy'][ip2[1]], s=175, c=P2['Utube'][ip2[1]], cmap=cmaps['temp']['cm'],
                    vmin=tmin, vmax=tmax, marker='^', edgecolors='k')
    ax[2,i].scatter(0, 0, s=75, c='k')
    ax[2,i].text(0, 0.4, 'RaXPol', fontsize=16, fontweight='bold')
    ax[2,i].legend([s1,s2], ['Probe 1', 'Probe 2'], loc='lower right')
    ax[2,i].set_title(f"{vol[vi[i]]['scan_time'][6]} UTC - {vol[vi[i]]['elev'][6]:.1f}\N{DEGREE SIGN}", fontsize=16, fontweight='bold')
    if cb_flag:
        plt.colorbar(s1, ax=ax[2,i], label='Sfc temperature (C)', extend='both')#, ticks=np.linspace(tmin,tmax,7))
    
    plot_cfill(vol[vi[i]]['xx'][6,:,:], vol[vi[i]]['yy'][6,:,:], np.ma.masked_array(vol[vi[i]]['vel'][6,:,:], vol[vi[i]]['dbz'][6,:,:]<1),
               'vel', ax[3,i], datalims=[-va,va], xlims=[-5,5], ylims=[-2,8], cbar=cb_flag)
    b1 = ax[3,i].barbs(P1['xx'][ip1[1]], P1['yy'][ip1[1]], P1['uwind'][ip1[1]], P1['vwind'][ip1[1]], barbcolor='k', length=9)
    b2 = ax[3,i].barbs(P2['xx'][ip2[1]], P2['yy'][ip2[1]], P2['uwind'][ip2[1]], P2['vwind'][ip2[1]], barbcolor='k', length=9)
    s1 = ax[3,i].scatter(P1['xx'][ip1[1]], P1['yy'][ip1[1]], s=150, c=P1['Utube'][ip1[1]], cmap=cmaps['temp']['cm'],
                    vmin=tmin, vmax=tmax, marker='s', edgecolors='k')
    s2 = ax[3,i].scatter(P2['xx'][ip2[1]], P2['yy'][ip2[1]], s=175, c=P2['Utube'][ip2[1]], cmap=cmaps['temp']['cm'],
                    vmin=tmin, vmax=tmax, marker='^', edgecolors='k')
    ax[3,i].scatter(0, 0, s=75, c='k')
    ax[3,i].text(0, 0.4, 'RaXPol', fontsize=16, fontweight='bold')
    ax[3,i].legend([s1,s2], ['Probe 1', 'Probe 2'], loc='lower right')
    if cb_flag:
        plt.colorbar(s1, ax=ax[3,i], label='Sfc temperature (C)', extend='both')#, ticks=np.linspace(tmin,tmax,7))

for i in range(len(vi)):
    ax[3,i].set_xlabel('E-W distance from radar (km)', fontsize=12)
    ax[i,0].set_ylabel('N-S distance from radar (km)', fontsize=12)


figsave = False
if figsave:
    plt.savefig(f"/Users/morgan.schneider/Documents/perils2023/iop2/figs/PPI_overview_with_mesonet_temp.png", dpi=300, bbox_inches='tight')
    # plt.savefig(f"/Users/morgan.schneider/Documents/conferences/perilsmeeting/iop2_vol{vi}_{filetime}_az{azimuth}_PPI.png", dpi=300)


#%% Save vortex positions to pickle

dat = dict()

for vi in np.arange(7,18):
    if vi == 7: # 082000 UTC
        filetime = vol[vi]['scan_time'][1]
        # vortex 1
        az_rot1 = np.array([290, 291, 292, 293, 294, 295,
                            296, 297, 298, 299, 300,
                            301, 302, 303, 304, 305,
                            306, 307, 308, 309, 310,
                            311, 312, 313, 314])
        r_rot1 = np.array([2.02, 2.04, 2.06, 2.08, 2.10, 2.11,
                           2.12, 2.13, 2.15, 2.16, 2.17,
                           2.19, 2.19, 2.20, 2.20, 2.21,
                           2.23, 2.23, 2.23, 2.24, 2.25,
                           2.27, 2.28, 2.30, 2.30])
        z_rot1 = np.array([0.18, 0.20, 0.20, 0.21, 0.23, 0.24,
                           0.26, 0.30, 0.32, 0.35, 0.38,
                           0.41, 0.44, 0.47, 0.48, 0.50,
                           0.56, 0.61, 0.63, 0.66, 0.68,
                           0.7, 0.72, 0.74, 0.76])
        x_rot1 = r_rot1 * np.sin(az_rot1*np.pi/180)
        y_rot1 = r_rot1 * np.cos(az_rot1*np.pi/180)
        datv1 = {'az':az_rot1, 'r':r_rot1, 'z':z_rot1, 'x':x_rot1, 'y':y_rot1}
        datv = {'vortex1':datv1}
        dat.update({f"{filetime}":datv})
    
    
    if vi == 8: # 082030 UTC
        filetime = vol[vi]['scan_time'][1]
        # vortex 1
        az_rot1 = np.array([300, 301, 302, 303, 304, 305,
                            306, 307, 308, 309, 310,
                            311, 312, 313, 314, 315,
                            316, 317, 318, 319, 320,
                            321, 322, 323, 324, 325])
        r_rot1 = np.array([2.12, 2.12, 2.12, 2.12, 2.13, 2.14,
                           2.14, 2.15, 2.15, 2.15, 2.15,
                           2.16, 2.17, 2.18, 2.19, 2.20,
                           2.20, 2.20, 2.21, 2.22, 2.24,
                           2.23, 2.23, 2.23, 2.23, 2.23])
        z_rot1 = np.array([0.04, 0.04, 0.06, 0.08, 0.10, 0.11,
                           0.12, 0.15, 0.18, 0.22, 0.26,
                           0.29, 0.32, 0.35, 0.38, 0.42,
                           0.46, 0.48, 0.51, 0.52, 0.57,
                           0.58, 0.60, 0.63, 0.68, 0.70])
        x_rot1 = r_rot1 * np.sin(az_rot1*np.pi/180)
        y_rot1 = r_rot1 * np.cos(az_rot1*np.pi/180)
        datv1 = {'az':az_rot1, 'r':r_rot1, 'z':z_rot1, 'x':x_rot1, 'y':y_rot1}
        
        # vortex 2
        az_rot2 = np.array([291, 292, 293, 294, 295,
                            296, 297, 298, 299, 300,
                            301, 302, 303, 304, 305,
                            306, 307, 308, 309, 310,
                            311, 312])
        r_rot2 = np.array([1.85, 1.85, 1.83, 1.82, 1.83,
                           1.83, 1.82, 1.81, 1.80, 1.80,
                           1.80, 1.80, 1.75, 1.74, 1.72,
                           1.70, 1.70, 1.71, 1.72, 1.75,
                           1.80, 1.83])
        z_rot2 = np.array([0.03, 0.04, 0.05, 0.07, 0.09,
                           0.12, 0.14, 0.16, 0.17, 0.19,
                           0.21, 0.25, 0.30, 0.32, 0.35,
                           0.38, 0.41, 0.42, 0.43, 0.44,
                           0.45, 0.48])
        x_rot2 = r_rot2 * np.sin(az_rot2*np.pi/180)
        y_rot2 = r_rot2 * np.cos(az_rot2*np.pi/180)
        datv2 = {'az':az_rot2, 'r':r_rot2, 'z':z_rot2, 'x':x_rot2, 'y':y_rot2}
        
        # vortex 3
        az_rot3 = np.array([320, 321, 322, 323, 324, 325, 326, 327])
        r_rot3 = np.array([2.70, 2.72, 2.73, 2.75, 2.79, 2.82, 2.85, 2.88])
        z_rot3 = np.array([0.50, 0.55, 0.62, 0.67, 0.72, 0.76, 0.78, 0.82])
        x_rot3 = r_rot3 * np.sin(az_rot3*np.pi/180)
        y_rot3 = r_rot3 * np.cos(az_rot3*np.pi/180)
        datv3 = {'az':az_rot3, 'r':r_rot3, 'z':z_rot3, 'x':x_rot3, 'y':y_rot3}
        
        # rotor?
        az_rot4 = np.array([281, 282, 283, 284, 285,
                            286, 287, 288, 289, 290,
                            291, 292, 293, 294, 295])
        r_rot4 = np.array([2.24, 2.18, 2.13, 2.08, 2.04,
                           2.01, 1.99, 1.97, 1.96, 1.95,
                           1.94, 1.92, 1.90, 1.90, 1.90])
        z_rot4 = np.array([0.18, 0.18, 0.20, 0.20, 0.16,
                           0.18, 0.20, 0.22, 0.25, 0.28,
                           0.31, 0.34, 0.37, 0.38, 0.40])
        x_rot4 = r_rot4 * np.sin(az_rot4*np.pi/180)
        y_rot4 = r_rot4 * np.cos(az_rot4*np.pi/180)
        datv4 = {'az':az_rot4, 'r':r_rot4, 'z':z_rot4, 'x':x_rot4, 'y':y_rot4}
        
        datv = {'vortex1':datv1, 'vortex2':datv2, 'vortex3':datv3, 'rotor':datv4}
        dat.update({f"{filetime}":datv})
    
    
    if vi == 9: # 082100 UTC
        filetime = vol[vi]['scan_time'][1]
        # vortex 1
        az_rot1 = np.array([321, 322, 323, 324, 325,
                            326, 327, 328, 329, 330,
                            331, 332, 333, 334, 335])
        r_rot1 = np.array([2.12, 2.14, 2.16, 2.18, 2.19,
                           2.20, 2.21, 2.21, 2.19, 2.18,
                           2.18, 2.17, 2.17, 2.16, 2.15])
        z_rot1 = np.array([0.15, 0.18, 0.21, 0.24, 0.28,
                           0.30, 0.32, 0.34, 0.36, 0.38,
                           0.40, 0.42, 0.45, 0.50, 0.54])
        x_rot1 = r_rot1 * np.sin(az_rot1*np.pi/180)
        y_rot1 = r_rot1 * np.cos(az_rot1*np.pi/180)
        datv1 = {'az':az_rot1, 'r':r_rot1, 'z':z_rot1, 'x':x_rot1, 'y':y_rot1}
        
        # vortex 2
        az_rot2 = np.array([305, 306, 307, 308, 309, 310,
                            311, 312, 313, 314, 315,
                            316, 317, 318, 319, 320,
                            321, 322, 323, 324, 325,
                            326, 327, 328, 329, 330])
        r_rot2 = np.array([1.72, 1.71, 1.71, 1.70, 1.70, 1.69,
                           1.68, 1.68, 1.68, 1.68, 1.67,
                           1.68, 1.68, 1.67, 1.68, 1.68,
                           1.68, 1.69, 1.70, 1.70, 1.72,
                           1.72, 1.72, 1.73, 1.74, 1.74])
        z_rot2 = np.array([0.0, 0.0, 0.0, 0.03, 0.04, 0.06,
                           0.08, 0.13, 0.16, 0.17, 0.18,
                           0.21, 0.22, 0.23, 0.26, 0.30,
                           0.33, 0.39, 0.41, 0.44, 0.46,
                           0.49, 0.53, 0.54, 0.56, 0.60])
        x_rot2 = r_rot2 * np.sin(az_rot2*np.pi/180)
        y_rot2 = r_rot2 * np.cos(az_rot2*np.pi/180)
        datv2 = {'az':az_rot2, 'r':r_rot2, 'z':z_rot2, 'x':x_rot2, 'y':y_rot2}
        
        # vortex 3
        az_rot3 = np.array([323, 324, 325, 326, 327, 328, 329, 330,
                            331, 332, 333, 334, 335,
                            336, 337, 338, 339, 340])
        r_rot3 = np.array([2.51, 2.55, 2.59, 2.62, 2.64, 2.66, 2.68, 2.69,
                           2.70, 2.72, 2.73, 2.74, 2.75,
                           2.76, 2.76, 2.77, 2.78, 2.79])
        z_rot3 = np.array([0.13, 0.15, 0.18, 0.22, 0.26, 0.31, 0.36, 0.42,
                           0.48, 0.55, 0.62, 0.68, 0.74,
                           0.79, 0.84, 0.88, 0.91, 0.93])
        x_rot3 = r_rot3 * np.sin(az_rot3*np.pi/180)
        y_rot3 = r_rot3 * np.cos(az_rot3*np.pi/180)
        datv3 = {'az':az_rot3, 'r':r_rot3, 'z':z_rot3, 'x':x_rot3, 'y':y_rot3}
        
        # vortex 4
        az_rot4 = np.array([327, 328, 329, 330, 331, 332, 333, 334, 335])
        r_rot4 = np.array([2.92, 2.95, 2.98, 3.01, 3.05, 3.08, 3.12, 3.16, 3.20])
        z_rot4 = np.array([0.05, 0.12, 0.28, 0.34, 0.40, 0.45, 0.50, 0.55, 0.62])
        x_rot4 = r_rot4 * np.sin(az_rot4*np.pi/180)
        y_rot4 = r_rot4 * np.cos(az_rot4*np.pi/180)
        datv4 = {'az':az_rot4, 'r':r_rot4, 'z':z_rot4, 'x':x_rot4, 'y':y_rot4}
        
        datv = {'vortex1':datv1, 'vortex2':datv2, 'vortex3':datv3, 'vortex4':datv4}
        dat.update({f"{filetime}":datv})
    
    
    if vi == 10: # 082130 UTC
        filetime = vol[vi]['scan_time'][1]
        # vortex 2
        az_rot2 = np.array([323, 324, 325, 326, 327, 328, 329, 330, 331, 332])
        r_rot2 = np.array([1.62, 1.63, 1.65, 1.66, 1.67, 1.66, 1.66, 1.66, 1.66, 1.66])
        z_rot2 = np.array([0.05, 0.06, 0.08, 0.10, 0.12, 0.13, 0.14, 0.15, 0.17, 0.18])
        x_rot2 = r_rot2 * np.sin(az_rot2*np.pi/180)
        y_rot2 = r_rot2 * np.cos(az_rot2*np.pi/180)
        datv2 = {'az':az_rot2, 'r':r_rot2, 'z':z_rot2, 'x':x_rot2, 'y':y_rot2}
        
        # vortex 3
        az_rot3 = np.array([330, 331, 332, 333, 334, 335, 336, 337])
        r_rot3 = np.array([2.59, 2.60, 2.61, 2.62, 2.63, 2.64, 2.65, 2.67])
        z_rot3 = np.array([0.04, 0.06, 0.08, 0.12, 0.16, 0.20, 0.23, 0.27])
        x_rot3 = r_rot3 * np.sin(az_rot3*np.pi/180)
        y_rot3 = r_rot3 * np.cos(az_rot3*np.pi/180)
        datv3 = {'az':az_rot3, 'r':r_rot3, 'z':z_rot3, 'x':x_rot3, 'y':y_rot3}
        
        # vortex 4
        az_rot4 = np.array([333, 334, 335, 336, 337, 338, 339])
        r_rot4 = np.array([2.93, 2.96, 2.99, 3.02, 3.05, 3.08, 3.11])
        z_rot4 = np.array([0.03, 0.07, 0.11, 0.18, 0.25, 0.27, 0.31])
        x_rot4 = r_rot4 * np.sin(az_rot4*np.pi/180)
        y_rot4 = r_rot4 * np.cos(az_rot4*np.pi/180)
        datv4 = {'az':az_rot4, 'r':r_rot4, 'z':z_rot4, 'x':x_rot4, 'y':y_rot4}
        
        datv = {'vortex2':datv2, 'vortex3':datv3, 'vortex4':datv4}
        dat.update({f"{filetime}":datv})
    
    
    if vi == 12: # 082230 UTC
        filetime = vol[vi]['scan_time'][1]
        # vortex 2
        az_rot2 = np.array([357, 358, 359, 0, 1, 2, 3, 4, 5, 6, 7, 8])
        r_rot2 = np.array([1.87, 1.91, 1.95, 2.00, 2.05, 2.10, 2.15, 2.20, 2.25, 2.32, 2.38, 2.45])
        z_rot2 = np.array([0.09, 0.12, 0.14, 0.16, 0.18, 0.22, 0.27, 0.33, 0.41, 0.48, 0.56, 0.65])
        x_rot2 = r_rot2 * np.sin(az_rot2*np.pi/180)
        y_rot2 = r_rot2 * np.cos(az_rot2*np.pi/180)
        datv2 = {'az':az_rot2, 'r':r_rot2, 'z':z_rot2, 'x':x_rot2, 'y':y_rot2}
        
        # vortex 3
        az_rot3 = np.array([352, 353,  354,  355,  356,  357])
        r_rot3 = np.array([2.95, 3.01, 3.07, 3.12, 3.18, 3.22])
        z_rot3 = np.array([0.16, 0.18, 0.21, 0.24, 0.27, 0.32])
        x_rot3 = r_rot3 * np.sin(az_rot3*np.pi/180)
        y_rot3 = r_rot3 * np.cos(az_rot3*np.pi/180)
        datv3 = {'az':az_rot3, 'r':r_rot3, 'z':z_rot3, 'x':x_rot3, 'y':y_rot3}
        
        # vortex 4
        az_rot4 = np.array([350, 351,  352,  353,  354,  355,  356,  357,  358,  359])
        r_rot4 = np.array([3.12, 3.14, 3.16, 3.18, 3.20, 3.25, 3.30, 3.38, 3.46, 3.55])
        z_rot4 = np.array([0.15, 0.18, 0.23, 0.26, 0.31, 0.38, 0.42, 0.47, 0.53, 0.60])
        x_rot4 = r_rot4 * np.sin(az_rot4*np.pi/180)
        y_rot4 = r_rot4 * np.cos(az_rot4*np.pi/180)
        datv4 = {'az':az_rot4, 'r':r_rot4, 'z':z_rot4, 'x':x_rot4, 'y':y_rot4}
        
        datv = {'vortex2':datv2, 'vortex3':datv3, 'vortex4':datv4}
        dat.update({f"{filetime}":datv})
    
    
    if vi == 13: # 082300 UTC
        filetime = vol[vi]['scan_time'][1]
        # vortex 2
        az_rot2 = np.array([356, 357, 358, 359, 0, 1, 2, 3, 4, 5])
        r_rot2 = np.array([1.77, 1.82, 1.85, 1.89, 1.93, 1.97, 2.01, 2.05, 2.08, 2.10])
        z_rot2 = np.array([0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.12, 0.15, 0.18])
        x_rot2 = r_rot2 * np.sin(az_rot2*np.pi/180)
        y_rot2 = r_rot2 * np.cos(az_rot2*np.pi/180)
        datv2 = {'az':az_rot2, 'r':r_rot2, 'z':z_rot2, 'x':x_rot2, 'y':y_rot2}
        
        # vortex 3
        az_rot3 = np.array([356, 357, 358, 359, 0, 1, 2, 3, 4, 5])
        r_rot3 = np.array([3.35, 3.40, 3.45, 3.52, 3.58, 3.67, 3.78, 3.90, 4.02, 4.16])
        z_rot3 = np.array([0.08, 0.10, 0.13, 0.20, 0.28, 0.38, 0.50, 0.58, 0.66, 0.75])
        x_rot3 = r_rot3 * np.sin(az_rot3*np.pi/180)
        y_rot3 = r_rot3 * np.cos(az_rot3*np.pi/180)
        datv3 = {'az':az_rot3, 'r':r_rot3, 'z':z_rot3, 'x':x_rot3, 'y':y_rot3}
        
        datv = {'vortex2':datv2, 'vortex3':datv3}
        dat.update({f"{filetime}":datv})
    
    
    if vi == 14: # 082330 UTC
        filetime = vol[vi]['scan_time'][1]
        # vortex 3
        az_rot3 = np.array([2, 3, 4, 5, 6, 7, 8])
        r_rot3 = np.array([3.79, 3.85, 3.92, 4.00, 4.10, 4.20, 4.30]) # 7-8 could be 4.50, 4.55
        z_rot3 = np.array([0.10, 0.12, 0.20, 0.30, 0.38, 0.46, 0.56]) # 7-8 could be 0.70, 0.75
        x_rot3 = r_rot3 * np.sin(az_rot3*np.pi/180)
        y_rot3 = r_rot3 * np.cos(az_rot3*np.pi/180)
        datv3 = {'az':az_rot3, 'r':r_rot3, 'z':z_rot3, 'x':x_rot3, 'y':y_rot3}
        
        datv = {'vortex3':datv3}
        dat.update({f"{filetime}":datv})
    
    
    if vi == 15: # 082400 UTC
        filetime = vol[vi]['scan_time'][1]
        # vortex 3
        az_rot3 = np.array([7, 8, 9, 10, 11])
        r_rot3 = np.array([4.20, 4.28, 4.40, 4.55, 4.70])
        z_rot3 = np.array([0.10, 0.13, 0.25, 0.35, 0.50])
        x_rot3 = r_rot3 * np.sin(az_rot3*np.pi/180)
        y_rot3 = r_rot3 * np.cos(az_rot3*np.pi/180)
        datv3 = {'az':az_rot3, 'r':r_rot3, 'z':z_rot3, 'x':x_rot3, 'y':y_rot3}
        
        datv = {'vortex3':datv3}
        dat.update({f"{filetime}":datv})
    
    
    if vi == 16: # 082430 UTC
        filetime = vol[vi]['scan_time'][1]
        az_rot3 = np.array([10, 11, 12, 13, 14])
        r_rot3 = np.array([4.80, 4.85, 4.90, 4.98, 5.05])
        z_rot3 = np.array([0.10, 0.15, 0.20, 0.30, 0.40])
        x_rot3 = r_rot3 * np.sin(az_rot3*np.pi/180)
        y_rot3 = r_rot3 * np.cos(az_rot3*np.pi/180)
        datv3 = {'az':az_rot3, 'r':r_rot3, 'z':z_rot3, 'x':x_rot3, 'y':y_rot3}
        
        datv = {'vortex3':datv3}
        dat.update({f"{filetime}":datv})
    
    
    if vi == 17: # 082500 UTC
        filetime = vol[vi]['scan_time'][1]
        # vortex 3
        az_rot3 = np.array([7, 8, 9, 10, 11, 12])
        r_rot3 = np.array([5.81, 5.86, 5.91, 5.97, 6.04, 6.11])
        z_rot3 = np.array([0.15, 0.20, 0.30, 0.48, 0.65, 0.85])
        x_rot3 = r_rot3 * np.sin(az_rot3*np.pi/180)
        y_rot3 = r_rot3 * np.cos(az_rot3*np.pi/180)
        datv3 = {'az':az_rot3, 'r':r_rot3, 'z':z_rot3, 'x':x_rot3, 'y':y_rot3}
        
        datv = {'vortex3':datv3}
        dat.update({f"{filetime}":datv})


save_to_pickle(dat, '/Users/morgan.schneider/Documents/perils2023/iop2/raxpol_vortex_locs.pkl')



#%% PERiLS 2023 Z, ZDR calibration constants

# IOP2
# KLZK beam height 1.2 km -- Rax range 3.7 km
fn_rax = '/Users/morgan.schneider/Documents/perils2023/iop2/raxpol/CFradial/cfrad.20230303_080850_RaXPol_v38_s13.nc' # 17.5 deg
fn_88d = '/Users/morgan.schneider/Documents/perils2023/iop2/klzk/KLZK20230303_080235_V06.ar2v'
# fn_rax = 'cfrad.20230303_081922_RaXPol_v97_s1.nc' # 19 deg
# fn_88d = '/Users/morgan.schneider/Documents/perils2023/iop2/klzk/KLZK20230303_081934_V06.ar2v'

raxpol = pyart.io.read(fn_rax)
klzk = pyart.io.read(fn_88d)

rax = read_raxpol(fn_rax)
r_lzk = klzk.range['data']/1000
az_lzk = klzk.azimuth['data']
el_lzk = np.median(klzk.elevation['data'][0:719])

r_mat = np.tile(r_lzk, (len(az_lzk), 1))
az_mat = np.transpose(np.tile(az_lzk, (len(r_lzk), 1)))

xx_lzk = r_mat * np.sin(az_mat*np.pi/180) * np.cos(el_lzk*np.pi/180)
yy_lzk = r_mat * np.cos(az_mat*np.pi/180) * np.cos(el_lzk*np.pi/180)
zz_lzk = r_mat * np.sin(el_lzk*np.pi/180)



xr,yr = latlon2xy(raxpol.latitude['data'][0], raxpol.longitude['data'][0], klzk.latitude['data'][0], klzk.longitude['data'][0])


c1 = pyart.graph.RadarDisplay(raxpol)
c2 = pyart.graph.RadarDisplay(klzk)



#%% ZDR correction stuff

fig = plt.figure(figsize=(14,5))

ax1 = fig.add_subplot(121)
c2.plot('reflectivity', 0, vmin=0, vmax=70, cmap='pyart_NWSRef')
c2.plot_range_rings([20,40,60,80,100])
ax2.set_xlim([-100,100])
ax2.set_ylim([-100,100])
# ax1.scatter(0, 0, s=50, c='k')
# ax1.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')

ax2 = fig.add_subplot(122)
c2.plot('velocity', vmin=-34, vmax=34, cmap='pyart_Carbone42')
c2.plot_range_rings([20,40,60,80,100])
ax2.set_xlim([-100,100])
ax2.set_ylim([-100,100])
ax2.scatter(xr, yr, s=25, c='k')

plt.show()

#%% ZDR correction stuff

fig = plt.figure(figsize=(14,5))

ax1 = fig.add_subplot(121)
c1.plot('DBZ', 0, vmin=0, vmax=70, cmap='pyart_NWSRef')
c1.plot_range_rings([5, 10])
ax1.set_xlim([-10,10])
ax1.set_ylim([-10,10])
# ax1.scatter(0, 0, s=50, c='k')
# ax1.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')

ax2 = fig.add_subplot(122)
c2.plot('reflectivity', vmin=0, vmax=70, cmap='pyart_NWSRef')
c2.plot_range_rings([20,40,60,80,100])
ax2.set_xlim([-100,100])
ax2.set_ylim([-100,100])
ax2.scatter(xr, yr, s=25, c='k')

plt.show()



fig = plt.figure(figsize=(14,5))

ax1 = fig.add_subplot(121)
c1.plot('DBZ', 0, vmin=0, vmax=70, cmap='pyart_NWSRef')
c1.plot_range_rings([5, 10])
ax1.set_xlim([-10,10])
ax1.set_ylim([-10,10])
# ax1.scatter(0, 0, s=50, c='k')
# ax1.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')

ax2 = fig.add_subplot(122)
c2.plot('reflectivity', vmin=0, vmax=70, cmap='pyart_NWSRef')
c2.plot_range_rings([20,40,60,80,100])
ax2.set_xlim([xr-10,xr+10])
ax2.set_ylim([yr-10,yr+10])
ax2.scatter(xr, yr, s=25, c='k')

plt.show()



fig = plt.figure(figsize=(14,5))

ax1 = fig.add_subplot(121)
c1.plot('ZDR', 0, vmin=-5, vmax=5, cmap='pyart_NWSRef')
c1.plot_range_rings([5, 10])
ax1.set_xlim([-10,10])
ax1.set_ylim([-10,10])
# ax1.scatter(0, 0, s=50, c='k')
# ax1.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')

ax2 = fig.add_subplot(122)
c2.plot('differential_reflectivity', vmin=-5, vmax=5, cmap='pyart_NWSRef')
c2.plot_range_rings([20,40,60,80,100])
ax2.set_xlim([-100,100])
ax2.set_ylim([-100,100])
ax2.scatter(xr, yr, s=25, c='k')

plt.show()



fig = plt.figure(figsize=(14,5))

ax1 = fig.add_subplot(121)
c1.plot('ZDR', 0, vmin=-5, vmax=5, cmap='pyart_NWSRef')
c1.plot_range_rings([5, 10])
ax1.set_xlim([-10,10])
ax1.set_ylim([-10,10])
# ax1.scatter(0, 0, s=50, c='k')
# ax1.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')

ax2 = fig.add_subplot(122)
c2.plot('differential_reflectivity', vmin=-5, vmax=5, cmap='pyart_NWSRef')
c2.plot_range_rings([20,40,60,80,100])
ax2.set_xlim([xr-10,xr+10])
ax2.set_ylim([yr-10,yr+10])
# ax2.set_xlim([-100,100])
# ax2.set_ylim([-100,100])
ax2.scatter(xr, yr, s=25, c='k')

plt.show()



#%% ZDR correction stuff

r_rax = (xr**2 + yr**2)**0.5
az_rax = np.arctan2(xr,yr)*180/np.pi
ir = np.where(r_lzk >= r_rax-4)[0][0]
# bh_lzk = np.median(klzk.gate_z['data'][0:359,ir])

# j1 = np.where(r_lzk >= r_rax-10)[0][-1]
# j2 = np.where(r_lzk <= r_rax+10)[0][0]

# i1 = np.where(az_lzk[0:719] >= 120)[0]
# i2 = np.where(az_lzk[0:719] <= 140)[0]

ii = slice(359,399)
jj = slice(318,358)

jj1 = slice(66,166)

# Z offset
zmean_lzk = np.min(klzk.fields['reflectivity']['data'][ii,jj])
zmean_rax = np.min(raxpol.fields['DBZ']['data'][:,jj1])

z_offset = zmean_lzk - zmean_rax

# ZDR offset
zdrmean_lzk = np.min(klzk.fields['differential_reflectivity']['data'][ii,jj])
zdrmean_rax = np.min(raxpol.fields['ZDR']['data'][:,jj1])

zdr_offset = zdrmean_lzk - zdrmean_rax




#%% ZDR correction stuff

raxpol2 = raxpol
raxpol2.fields['DBZ']['data'] = raxpol2.fields['DBZ']['data'] + z_offset
raxpol2.fields['ZDR']['data'] = raxpol2.fields['ZDR']['data'] + zdr_offset

c3 = pyart.graph.RadarDisplay(raxpol2)



# fig = plt.figure(figsize=(14,5))

# ax1 = fig.add_subplot(121)
# c1.plot('DBZ', 0, vmin=0, vmax=70, cmap='pyart_NWSRef')
# c1.plot_range_rings([5, 10])
# ax1.set_xlim([-10,10])
# ax1.set_ylim([-10,10])
# # ax1.scatter(0, 0, s=50, c='k')
# # ax1.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')

# ax2 = fig.add_subplot(122)
# c2.plot('reflectivity', vmin=0, vmax=70, cmap='pyart_NWSRef')
# c2.plot_range_rings([20,40,60,80,100])
# ax2.set_xlim([-100,100])
# ax2.set_ylim([-100,100])
# ax2.scatter(xr, yr, s=25, c='k')

# plt.show()



fig = plt.figure(figsize=(14,5))

ax1 = fig.add_subplot(121)
c3.plot('DBZ', vmin=0, vmax=70, cmap='pyart_NWSRef')
c3.plot_range_rings([5, 10])
ax1.set_xlim([-10,10])
ax1.set_ylim([-10,10])
# ax1.scatter(0, 0, s=50, c='k')
# ax1.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')

ax2 = fig.add_subplot(122)
c2.plot('reflectivity', vmin=0, vmax=70, cmap='pyart_NWSRef')
c2.plot_range_rings([20,40,60,80,100])
ax2.set_xlim([xr-10,xr+10])
ax2.set_ylim([yr-10,yr+10])
ax2.scatter(xr, yr, s=25, c='k')

plt.show()



# fig = plt.figure(figsize=(14,5))

# ax1 = fig.add_subplot(121)
# c1.plot('ZDR', 0, vmin=-5, vmax=5, cmap='pyart_NWSRef')
# c1.plot_range_rings([5, 10])
# ax1.set_xlim([-10,10])
# ax1.set_ylim([-10,10])
# # ax1.scatter(0, 0, s=50, c='k')
# # ax1.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')

# ax2 = fig.add_subplot(122)
# c2.plot('differential_reflectivity', vmin=-5, vmax=5, cmap='pyart_NWSRef')
# c2.plot_range_rings([20,40,60,80,100])
# ax2.set_xlim([-100,100])
# ax2.set_ylim([-100,100])
# ax2.scatter(xr, yr, s=25, c='k')

# plt.show()



fig = plt.figure(figsize=(14,5))

ax1 = fig.add_subplot(121)
c3.plot('ZDR', vmin=-5, vmax=5, cmap='pyart_NWSRef')
c3.plot_range_rings([5, 10])
ax1.set_xlim([-10,10])
ax1.set_ylim([-10,10])
# ax1.scatter(0, 0, s=50, c='k')
# ax1.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')

ax2 = fig.add_subplot(122)
c2.plot('differential_reflectivity', vmin=-5, vmax=5, cmap='pyart_NWSRef')
c2.plot_range_rings([20,40,60,80,100])
ax2.set_xlim([xr-10,xr+10])
ax2.set_ylim([yr-10,yr+10])
# ax2.set_xlim([-100,100])
# ax2.set_ylim([-100,100])
ax2.scatter(xr, yr, s=25, c='k')

plt.show()


#%% ZDR correction stuff
# IOP5
# KNQA beam height 1.2-1.5 km -- Rax range 2.5-3.1 km (29 deg), 3.1-4 km (23 deg), 4.6-5.8 km (15 deg)
fp = '/Users/morgan.schneider/Documents/perils2023/iop5/raxpol/CFradial/'
fn_rax = fp+'cfrad.20230405_163236_RaXPol_v136_s3.nc'
fn_88d = '/Users/morgan.schneider/Documents/perils2023/iop5/knqa/KNQA20230405_163216_V06.ar2v'
# knqa 162222, rax 162212 v123 s3
# knqa 162702, rax 162700 v129 s3
# knqa 163216, rax 163236 v136 s3/163148 v135 s3
raxpol = pyart.io.read(fn_rax)
knqa = pyart.io.read(fn_88d)
c1 = pyart.graph.RadarDisplay(raxpol)
c2 = pyart.graph.RadarDisplay(knqa)


xr,yr = latlon2xy(raxpol.latitude['data'][0], raxpol.longitude['data'][0], knqa.latitude['data'][0], knqa.longitude['data'][0])


fig = plt.figure(figsize=(14,5))

ax1 = fig.add_subplot(121)
c2.plot('reflectivity', sweep=0, vmin=0, vmax=70, cmap='pyart_NWSRef')
c2.plot_range_rings([25,50,75,100])
# ax1.set_xlim([-150,150])
# ax1.set_ylim([-150,150])
ax1.set_xlim([xr-15,xr+15])
ax1.set_ylim([yr-15,yr+15])
ax1.scatter(xr, yr, s=25, c='k')
ax1.scatter(xr-15, yr, s=10, c='k')
ax1.scatter(xr, yr+15, s=10, c='k')
# ax1.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')

ax2 = fig.add_subplot(122)
c2.plot('differential_reflectivity', sweep=0, vmin=-6, vmax=6, cmap='pyart_NWSRef')
c2.plot_range_rings([25,50,75,100])
# ax2.set_xlim([-150,150])
# ax2.set_ylim([-150,150])
ax2.set_xlim([xr-15,xr+15])
ax2.set_ylim([yr-15,yr+15])
ax2.scatter(xr, yr, s=25, c='k')

plt.show()



fig = plt.figure(figsize=(14,5))

ax1 = fig.add_subplot(121)
c1.plot('DBZ', vmin=0, vmax=70, cmap='pyart_NWSRef')
c1.plot_range_rings([5,10,15])
ax1.set_xlim([-15,15])
ax1.set_ylim([-15,15])
# ax1.scatter(xr, yr, s=25, c='k')
# ax1.scatter(xr-10, yr, s=10, c='k')
# ax1.scatter(xr, yr+10, s=10, c='k')
# ax1.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')

ax2 = fig.add_subplot(122)
c1.plot('ZDR', vmin=-6, vmax=6, cmap='pyart_NWSRef')
c1.plot_range_rings([5,10,15])
ax2.set_xlim([-15,15])
ax2.set_ylim([-15,15])
# ax2.scatter(xr, yr, s=25, c='k')

plt.show()

#%% ZDR correction stuff

r_nqa = knqa.range['data']/1000
az_nqa = knqa.azimuth['data']
el_nqa = np.median(knqa.elevation['data'][0:719])

r_mat = np.tile(r_nqa, (len(az_nqa), 1))
az_mat = np.transpose(np.tile(az_nqa, (len(r_nqa), 1)))

xx_nqa = r_mat * np.sin(az_mat*np.pi/180) * np.cos(el_nqa*np.pi/180)
yy_nqa = r_mat * np.cos(az_mat*np.pi/180) * np.cos(el_nqa*np.pi/180)
zz_nqa = r_mat * np.sin(el_nqa*np.pi/180)



r_rax = (xr**2 + yr**2)**0.5
az_rax = np.arctan2(xr,yr)*180/np.pi + 360
ir = np.where(r_nqa >= r_rax-4)[0][0]
# bh_lzk = np.median(klzk.gate_z['data'][0:359,ir])

# j1 = np.where(r_lzk >= r_rax-10)[0][-1]
# j2 = np.where(r_lzk <= r_rax+10)[0][0]

# i1 = np.where(az_lzk[0:719] >= 120)[0]
# i2 = np.where(az_lzk[0:719] <= 140)[0]

ii = slice(52,88)
jj = slice(370,430)

ii1 = slice(88,268)
ii2 = slice(448,628)
jj1 = slice(333,400)

# Z offset
zmean_nqa = np.mean(knqa.fields['reflectivity']['data'][ii,jj])
zmean_rax = np.mean(raxpol.fields['DBZ']['data'][ii1,jj1])

z_offset = zmean_nqa - zmean_rax

# ZDR offset
zdrmean_nqa = np.mean(knqa.fields['differential_reflectivity']['data'][ii,jj])
zdrmean_rax = np.mean(raxpol.fields['ZDR']['data'][ii1,jj1])

zdr_offset = zdrmean_nqa - zdrmean_rax

#%% ZDR correction stuff

raxpol2 = pyart.io.read(fn_rax)
raxpol2.fields['DBZ']['data'] = raxpol2.fields['DBZ']['data'] + z_offset
raxpol2.fields['ZDR']['data'] = raxpol2.fields['ZDR']['data'] + zdr_offset
c3 = pyart.graph.RadarDisplay(raxpol2)


fig = plt.figure(figsize=(14,5))

ax1 = fig.add_subplot(121)
c2.plot('reflectivity', sweep=0, vmin=0, vmax=70, cmap='pyart_NWSRef')
c2.plot_range_rings([25,50,75,100])
# ax1.set_xlim([-150,150])
# ax1.set_ylim([-150,150])
ax1.set_xlim([xr-15,xr+15])
ax1.set_ylim([yr-15,yr+15])
ax1.scatter(xr, yr, s=25, c='k')
ax1.scatter(xr-15, yr, s=10, c='k')
ax1.scatter(xr, yr+15, s=10, c='k')
# ax1.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')

ax2 = fig.add_subplot(122)
c2.plot('differential_reflectivity', sweep=0, vmin=-6, vmax=6, cmap='pyart_NWSRef')
c2.plot_range_rings([25,50,75,100])
# ax2.set_xlim([-150,150])
# ax2.set_ylim([-150,150])
ax2.set_xlim([xr-15,xr+15])
ax2.set_ylim([yr-15,yr+15])
ax2.scatter(xr, yr, s=25, c='k')

plt.show()



fig = plt.figure(figsize=(14,5))

ax1 = fig.add_subplot(121)
c3.plot('DBZ', vmin=0, vmax=70, cmap='pyart_NWSRef')
c3.plot_range_rings([5,10,15])
ax1.set_xlim([-15,15])
ax1.set_ylim([-15,15])
# ax1.scatter(xr, yr, s=25, c='k')
# ax1.scatter(xr-10, yr, s=10, c='k')
# ax1.scatter(xr, yr+10, s=10, c='k')
# ax1.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')

ax2 = fig.add_subplot(122)
c3.plot('ZDR', vmin=-6, vmax=6, cmap='pyart_NWSRef')
c3.plot_range_rings([5,10,15])
ax2.set_xlim([-15,15])
ax2.set_ylim([-15,15])
# ax2.scatter(xr, yr, s=25, c='k')

plt.show()







#%% PERiLS 2023 IOP5 heading correction

# if False:
#     fp = '/Users/morgan.schneider/Documents/perils2023/iop5/raxpol/CFradial-corrected/'
#     files = sorted(glob(fp+'*.nc'))
    
#     for file in files:
#         ds = nc.Dataset(file, 'r+')
        
#         az = ds.variables['azimuth'][:] + 90
#         az_new = np.where(az > 360, az - 360, az)
#         ds.variables['azimuth'][:] = az_new
        
#         # az = ds.variables['azimuth'][:] - 90
#         # az_new = np.where(az < 0, az + 360, az)
#         # ds.variables['azimuth'][:] = az_new
        
#         ds.close()


#%% PERiLS 2023 IOP2 and IOP5 Z/ZDR calibration


z_offset = -6
zdr_offset = -3.5


fp = '/Users/morgan.schneider/Documents/perils2023/iop2/raxpol/CFradial_cal/'
files = sorted(glob(fp+'*.nc'))

for file in files:
    ds = nc.Dataset(file, 'r+')
    
    z_new = ds.variables['DBZ'][:] + z_offset
    zdr_new = ds.variables['ZDR'][:] + zdr_offset
    ds.variables['DBZ'][:] = z_new
    ds.variables['ZDR'][:] = zdr_new
    
    ds.close()


fp = '/Users/morgan.schneider/Documents/perils2023/iop5/raxpol/CFradial_cal/'
files = sorted(glob(fp+'*.nc'))

for file in files:
    ds = nc.Dataset(file, 'r+')
    
    z_new = ds.variables['DBZ'][:] + z_offset
    zdr_new = ds.variables['ZDR'][:] + zdr_offset
    ds.variables['DBZ'][:] = z_new
    ds.variables['ZDR'][:] = zdr_new
    
    ds.close()















