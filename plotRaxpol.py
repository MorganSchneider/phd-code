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
filetime = '082206'
fn = glob(fp+f"raxpol/CFradial_cal/*_{filetime}_*.nc")[0]

# rax = read_raxpol(fn)
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
c.plot('DBZ', 0, vmin=0, vmax=70, cmap='NWSRef')
# c.plot('VEL', 0, vmin=-30, vmax=30, cmap='balance')
# c.plot('RHOHV', 0, vmin=0.8, vmax=1, cmap='LangRainbow12')
# c.plot('ZDR', 0, vmin=-5, vmax=5, cmap=nwsdmap)
# c.plot_range_rings([2.5, 5, 7.5, 10])
ax.set_xlim([-5,5])
ax.set_ylim([-5,5])
plt.text(0, 0.2, 'RaXPol', fontsize=12, fontweight='bold')
# b1 = ax.barbs(P1['xx'][i1], P1['yy'][i1], P1['uwind'][i1], P1['vwind'][i1], barbcolor='k', length=9)
# b2 = ax.barbs(P2['xx'][i2], P2['yy'][i2], P2['uwind'][i2], P2['vwind'][i2], barbcolor='k', length=9)
# s1 = ax.scatter(P1['xx'][i1], P1['yy'][i1], s=150, c=P1['Tfast'][i1], cmap=cmaps['temp']['cm'],
#                 vmin=datalims[0], vmax=datalims[1], marker='s', edgecolors='k')
# s2 = ax.scatter(P2['xx'][i2], P2['yy'][i2], s=175, c=P2['Tfast'][i2], cmap=cmaps['temp']['cm'],
#                 vmin=datalims[0], vmax=datalims[1], marker='^', edgecolors='k')
ax.scatter(0, 0, s=30, c='k')
# ax.scatter(-2.6, 0.6, s=20, c='k', marker='x')
# plt.colorbar(s1,label='Sfc temperature (C)')
# plt.legend([s1,s2], ['Probe 1', 'Probe 2'], loc='upper right')
plt.show()

if figsave:
    plt.savefig(fp+f"figs/vel+T_{filetime}.png", dpi=400)


#%% COW data IOP2

cow_times = ['074052', '074952', '080022', '080922', '081952', '082422', '083022', '083922']
cow_filetime = '083922'
cow_fn = glob(fp+f"cow/*_{cow_filetime}.*.nc")[0]
cow = pyart.io.read(cow_fn)


xrax,yrax = latlon2xy(raxpol.latitude['data'][0], raxpol.longitude['data'][0],
                      cow.latitude['data'][0], cow.longitude['data'][0])


c2 = pyart.graph.RadarDisplay(cow)

fig = plt.figure(figsize=(7.5,6))
ax = fig.add_subplot(111)
c2.plot('DBZHCC', 0, vmin=0, vmax=70, cmap='NWSRef')
# c2.plot('VEL', 0, vmin=-30, vmax=30, cmap='Carbone42')
# c2.plot('RHOHV', 0, vmin=0.8, vmax=1, cmap='LangRainbow12')
# c2.plot('ZDRC', 0, vmin=-5, vmax=5, cmap=nwsdmap)
# c2.plot_range_rings([20, 40, 60, 80])
ax.set_xlim([-70,10])
ax.set_ylim([-10,70])
ax.scatter(0, 0, s=40, c='k')
plt.text(-5, 2, 'COW', fontsize=16, fontweight='bold')
ax.scatter(xrax, yrax, s=40, c='k')
plt.text(xrax-8, yrax+2, 'RaXPol', fontsize=16, fontweight='bold')


if figsave:
    plt.savefig(fp+f"figs/cow-imgs/dbz_{cow_filetime}.png", dpi=300)


fig = plt.figure(figsize=(7.5,6))
ax = fig.add_subplot(111)
c2.plot('VEL', 0, vmin=-30, vmax=30, cmap='balance')
# c2.plot_range_rings([20, 40, 60, 80])
ax.set_xlim([-70,10])
ax.set_ylim([-10,70])
ax.scatter(0, 0, s=40, c='k')
plt.text(-5, 2, 'COW', fontsize=16, fontweight='bold')
ax.scatter(xrax, yrax, s=40, c='k')
plt.text(xrax-8, yrax+2, 'RaXPol', fontsize=16, fontweight='bold')

if figsave:
    plt.savefig(fp+f"figs/cow-imgs/vel_{cow_filetime}.png", dpi=300)



#%% Make COW overview figures??? idk what i want from this

fp = '/Users/morgan.schneider/Documents/perils2023/iop2/'
filetime = '082028'
fn = glob(fp+f"raxpol/CFradial_cal/*_{filetime}_*.nc")[0]

raxpol = pyart.io.read(fn)
rax_lat = raxpol.latitude['data'][0]
rax_lon = raxpol.longitude['data'][0]


cow_times = ['074952', '080022', '081952', '082422', '083022']
cow_fn = [glob(fp+f"cow/*_{cowt}.*.nc")[0] for cowt in cow_times]


cow = {}
for i in range(len(cow_fn)):
    tmp = pyart.io.read(cow_fn[i])
    
    cow_lat = tmp.latitude['data'][0]
    cow_lon = tmp.longitude['data'][0]
    
    dat = {'dbz':tmp.fields['DBZHCC'][:].data, 'vel':tmp.fields['VEL'][:].data, 'r':cow.range['data']/1000,
           'az':cow.azimuth['data'], 'elev':cow.elevation['data']}
    cow.update({f"{cow_times[i]}":dat})


xrax,yrax = latlon2xy(rax_lat, rax_lon, cow_lat, cow_lon)



fig,ax = plt.subplots(figsize=(10,10))






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
c.plot('DBZ', 0, vmin=0, vmax=70, cmap='NWSRef')
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
c.plot('VEL', 0, vmin=-30, vmax=30, cmap='Carbone42')
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
    zvort_tmp = np.zeros(shape=(len(f),360,1246,), dtype=float) # vertical pseudovorticity
    hvort_tmp = np.zeros(shape=(len(f),360,1246,), dtype=float) # cross-radial horizontal pseudovorticity
    vort3d_tmp = np.zeros(shape=(len(f),360,1246,), dtype=float) # "3D" (really 2D) pseudovorticity magnitude
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
                zvort_tmp[ii,int(ix),:] = 1/(r*1000) * (v2-v1)/(np.pi/180) # 2/r * dVr/dphi
            else:
                az_tmp[ii,int(ix)] = ix
                xx_tmp[ii,int(ix),:] = r * np.sin(ix*np.pi/180) * np.cos(d['elev']*np.pi/180)
                yy_tmp[ii,int(ix),:] = r * np.cos(ix*np.pi/180) * np.cos(d['elev']*np.pi/180)
                zz_tmp[ii,int(ix),:] = r * np.sin(d['elev']*np.pi/180)
    
    hvort_tmp = -2/(r*1000) * np.gradient(vel_tmp, el_tmp*np.pi/180, axis=0) # -2/r * dVr/dtheta
    vort3d_tmp = np.sqrt(zvort_tmp**2 + hvort_tmp**2)
    div_tmp = np.gradient(vel_tmp, r*1000, axis=2)
    
    vol[vn].update({'dbz':dbz_tmp, 'vel':vel_tmp, 'sw':sw_tmp, 'zdr':zdr_tmp, 'rhohv':rhohv_tmp,
                    'xx':xx_tmp, 'yy':yy_tmp, 'zz':zz_tmp, 'az':az_tmp, 'elev':el_tmp, 'r':r,
                    'zvort':zvort_tmp, 'hvort':hvort_tmp, 'vort3d':vort3d_tmp, 'div':div_tmp,
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

if True:
    print("...Saving to pickle...")
    dbfile = open('/Users/morgan.schneider/Documents/perils2023/iop2/circuit_data.pkl', 'wb')
    data = {'Rax':vol, 'P1':P1, 'P2':P2}
    pickle.dump(data, dbfile)
    dbfile.close()
    del data

del az_tmp,dbz_tmp,vel_tmp,sw_tmp,zdr_tmp,rhohv_tmp,xx_tmp,yy_tmp,zz_tmp,zvort_tmp,hvort_tmp,vort3d_tmp,div_tmp,el_tmp,fname_tmp,time_tmp,files,xx1,yy1,xx2,yy2,P1,P2,d

#%% IOP2 plot reconstructed RHIs and azimuth-height cross sections

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

if 'locs' not in locals():
    dbfile = open('/Users/morgan.schneider/Documents/perils2023/iop2/raxpol_vortex_locs.pkl', 'rb')
    locs = pickle.load(dbfile)
    dbfile.close()

# i1 = np.where(P1['time'] == int(vol[vi]['scan_time'][eli]))[0][0]
# i2 = np.where(P2['time'] == int(vol[vi]['scan_time'][eli]))[0][0]

T_lims = [18,21] # Temperature
Td_lims = [17,19] # Dewpoint
datalims = T_lims
vort_lim = 0.12
div_lim = 0.1

va = vol[0]['va']


vi = 9
eli = 1
filetime = vol[vi]['scan_time'][1]
vortex_num = 1
# vortex 1: vi = 7-9
# vortex 2: vi = 8-13
# vortex 3: vi = 8-17
# vortex 4: vi = 9-12
# rotor: vi = 8

az_rot = locs[filetime][f"vortex{vortex_num}"]['az']
r_rot = locs[filetime][f"vortex{vortex_num}"]['r']
# z_rot = locs[filetime][f"vortex{vortex_num}"]['z']
# x_rot = locs[filetime][f"vortex{vortex_num}"]['x']
# y_rot = locs[filetime][f"vortex{vortex_num}"]['y']

x_rot = np.array([])
y_rot = np.array([])

for key in list(locs[filetime].keys()):
    if 'vortex' in key:
        x_rot = np.append(x_rot, locs[filetime][key]['x'])
        y_rot = np.append(y_rot, locs[filetime][key]['y'])

# vortex_num = '-all'

# # azimuth = round(np.mean(az_rot))
azimuth = az_rot[round(len(az_rot)/2)]
irot = np.where(np.isclose(az_rot, azimuth))[0][0]

# # if azimuth not in az_rot:
# #     print(f"Invalid azimuth, must be between {az_rot[0]}-{az_rot[-1]} degrees")
# #     azimuth = az_rot[0]


azi = np.where(vol[vi]['az'][eli,:].round(0) == azimuth)[0][0]
rr = (vol[vi]['xx'][:,azi,:]**2 + vol[vi]['yy'][:,azi,:]**2)**0.5



if False:
    if vi == 7: # 082000
        xl = [-2.5, -1.0]
        yl = [-0.5, 2.5]
        
    elif vi == 8: # 082030
        xl = [-2.25, -0.75]
        yl = [0.0, 3.0]
        
    elif vi == 9: # 082100
        xl = [-2.0, -0.5]
        yl = [0.5, 3.5]
        
    elif vi == 10: # 082130
        xl = [-1.75, -0.25]
        yl = [1.0, 4.0]
        
    elif vi == 12: # 082230
        xl = [-1.0, 0.5]
        yl = [1.25, 4.25]
        
    elif vi == 13: # 082300
        xl = [-0.75, 0.75]
        yl = [1.5, 4.5]
        
    elif vi == 14: # 082330
        xl = [-0.25, 1.25]
        yl = [2.5, 5.5]
        
    elif vi == 15: # 082400
        xl = [0.0, 1.5]
        yl = [3.0, 6.0]
        
    elif vi == 16: # 082430
        xl = [0.25, 1.75]
        yl = [3.5, 6.5]
        
    elif vi == 17: # 082500
        xl = [0.25, 1.75]
        yl = [4.0, 7.0]


if vi <= 15:
    rlim = 6
    zlim = 2
else:
    rlim = 8
    zlim = 2.5



#%% Actual dissertation/paper figures

xls = [[-2.75, -1.25], [-2.5, -1.0], [-2.25, -0.75], [-2.0, -0.5], [-1.25, 0.25], 
       [-0.75, 0.75], [-0.5, 1.0], [0.0, 1.5], [0.25, 1.75], [0.25, 1.75]]
yls = [[0.5, 3.5], [0.5, 3.5], [1.0, 4.0], [1.0, 4.0], [1.5, 4.5], 
       [2.0, 5.0], [2.5, 5.5], [3.0, 6.0], [3.5, 6.5], [4.5, 7.5]]
sharex = False; sharey = False

xls = [[-5.0, 0.0], [-4.5, 0.5], [-4.0, 1.0], [-3.5, 1.5], [-3.0, 2.0], 
       [-2.5, 2.5], [-2.0, 3.0], [-1.5, 3.5], [-1.0, 4.0], [-0.5, 4.5]]
yls = [[-5.0, 5.0], [-4.5, 5.5], [-4.0, 6.0], [-3.5, 6.5], [-3.0, 7.0], 
       [-2.5, 7.5], [-2.0, 8.0], [-1.5, 8.5], [-1.0, 9.0], [-0.5, 9.5]]
sharex = False; sharey = False

# xls = [[-5.0, 5.0], [-5.0, 5.0], [-5.0, 5.0], [-5.0, 5.0], [-5.0, 5.0],
#        [-5.0, 5.0], [-5.0, 5.0], [-5.0, 5.0], [-5.0, 5.0], [-5.0, 5.0]]
# yls = [[-10.0, 10.0], [-10.0, 10.0], [-10.0, 10.0], [-10.0, 10.0], [-10.0, 10.0],
#        [-10.0, 10.0], [-10.0, 10.0], [-10.0, 10.0], [-10.0, 10.0], [-10.0, 10.0]]
# sharex = True; sharey = True



# dbz overview, vel overview
overview_flag = [0, 0]


figsave = False


# DBZ overview figure - 2x5 panels
if overview_flag[0]:
    # figsize (13,9) for the 5x10 km and 10x20 km, (13.5,9) for the 1.5x3 km
    fig,ax = plt.subplots(2, 5, figsize=(13,9), sharex=sharex, sharey=sharey, subplot_kw=dict(box_aspect=2), layout='constrained')
    
    vis = [7, 8, 9, 10, 12, 13, 14, 15, 16, 17]
    panels = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)', 'g)', 'h)', 'i)', 'j)']
    
    for i in range(5):
        if i == 4:
            cbflag = True
        else:
            cbflag = False
        
        plot_cfill(vol[vis[i]]['xx'][eli,:,:], vol[vis[i]]['yy'][eli,:,:], np.ma.masked_array(vol[vis[i]]['dbz'][eli,:,:], vol[vis[i]]['dbz'][eli,:,:]<1), 'dbz', ax[0,i],
                   datalims=[0,70], xlims=xls[i], ylims=yls[i], cbar=cbflag, cbfs=16, cbts=14)
        ax[0,i].scatter(0, 0, s=200, c='k', marker='.')
        if xls[0][1] - xls[0][0] == 5:
            ax[0,i].text(xls[i][0]+0.1, yls[i][-1]-1, f"{panels[i]}", fontsize=28, fontweight='bold')
        elif xls[0][1] - xls[0][0] == 10:
            ax[0,i].text(xls[i][0]+0.3, yls[i][-1]-2, f"{panels[i]}", fontsize=28, fontweight='bold')
        else:
            ax[0,i].text(xls[i][0]+0.05, yls[i][-1]-0.3, f"{panels[i]}", fontsize=28, fontweight='bold')
        ax[0,i].set_title(f"{vol[vis[i]]['scan_time'][1]} UTC", fontsize=18)
    
    for i in range(5):
        if i == 4:
            cbflag = True
        else:
            cbflag = False
        
        plot_cfill(vol[vis[i+5]]['xx'][eli,:,:], vol[vis[i+5]]['yy'][eli,:,:], np.ma.masked_array(vol[vis[i+5]]['dbz'][eli,:,:], vol[vis[i+5]]['dbz'][eli,:,:]<1), 'dbz', ax[1,i],
                   datalims=[0,70], xlims=xls[i+5], ylims=yls[i+5], cbar=cbflag, cbfs=16, cbts=14)
        ax[1,i].scatter(0, 0, s=200, c='k', marker='.')
        if xls[0][1] - xls[0][0] == 5:
            ax[1,i].text(xls[i+5][0]+0.2, yls[i+5][-1]-1, f"{panels[i+5]}", fontsize=28, fontweight='bold')
        elif xls[0][1] - xls[0][0] == 10:
            ax[1,i].text(xls[i+5][0]+0.4, yls[i+5][-1]-2, f"{panels[i+5]}", fontsize=28, fontweight='bold')
        else:
            ax[1,i].text(xls[i+5][0]+0.05, yls[i+5][-1]-0.3, f"{panels[i+5]}", fontsize=28, fontweight='bold')
        ax[1,i].set_title(f"{vol[vis[i+5]]['scan_time'][1]} UTC", fontsize=18)
        ax[1,i].set_xlabel('E-W distance (km)', fontsize=16)
    
    ax[0,0].set_ylabel('N-S distance (km)', fontsize=16)
    ax[1,0].set_ylabel('N-S distance (km)', fontsize=16)
    
    if figsave:
        plt.savefig(ip+'overview_dbz_PPI_10x20km.png', dpi=300)



# VEL overview figure
if overview_flag[1]:
    # figsize (13,9) for the 5x10 km and 10x20 km, (13.5,9) for the 1.5x3 km
    fig,ax = plt.subplots(2, 5, figsize=(13,9), sharex=sharex, sharey=sharey, subplot_kw=dict(box_aspect=2), layout='constrained')
    
    vis = [7, 8, 9, 10, 12, 13, 14, 15, 16, 17]
    panels = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)', 'g)', 'h)', 'i)', 'j)']
    
    for i in range(5):
        # xrot = np.array([])
        # yrot = np.array([])
        # filetime = vol[vis[i]]['scan_time'][1]
        # for key in list(locs[filetime].keys()):
        #     if 'vortex' in key:
        #         xrot = np.append(xrot, locs[filetime][key]['x'])
        #         yrot = np.append(yrot, locs[filetime][key]['y'])
        # if 'vortex3' in locs[filetime].keys():
        #     xrot3 = locs[filetime]['vortex3']['x']
        #     yrot3 = locs[filetime]['vortex3']['y']
        # if 'vortex4' in locs[filetime].keys():
        #     xrot4 = locs[filetime]['vortex4']['x']
        #     yrot4 = locs[filetime]['vortex4']['y']
                
        if i == 4:
            cbflag = True
        else:
            cbflag = False
        
        plot_cfill(vol[vis[i]]['xx'][eli,:,:], vol[vis[i]]['yy'][eli,:,:], np.ma.masked_array(vol[vis[i]]['vel'][eli,:,:], vol[vis[i]]['vel'][eli,:,:]<-100), 'vel', ax[0,i],
                   datalims=[-30,30], xlims=xls[i], ylims=yls[i], cbar=cbflag, cbfs=16, cbts=14, cmap='balance')
        ax[0,i].scatter(0, 0, s=200, c='k', marker='.')
        # ax[0,i].scatter(xrot, yrot, s=40, c='k', marker='.')
        # if 'vortex4' in locs[filetime].keys():
        #     ax[0,i].scatter(xrot4, yrot4, s=30, fc='w', ec='k', marker='o', linewidth=1)
        # elif 'vortex3' in locs[filetime].keys():
        #     ax[0,i].scatter(xrot3, yrot3, s=30, fc='w', ec='k', marker='o', linewidth=1)
        if xls[0][1] - xls[0][0] == 5:
            ax[0,i].text(xls[i][0]+0.1, yls[i][-1]-1, f"{panels[i]}", color='k', fontsize=28, fontweight='bold')
        elif xls[0][1] - xls[0][0] == 10:
            ax[0,i].text(xls[i][0]+0.3, yls[i][-1]-2, f"{panels[i]}", color='k', fontsize=28, fontweight='bold')
        else:
            ax[0,i].text(xls[i][0]+0.05, yls[i][-1]-0.3, f"{panels[i]}", color='k', fontsize=28, fontweight='bold')
        ax[0,i].set_title(f"{vol[vis[i]]['scan_time'][1]} UTC", fontsize=18)
    
    for i in range(5):
        # xrot = np.array([])
        # yrot = np.array([])
        # filetime = vol[vis[i+5]]['scan_time'][1]
        # for key in list(locs[filetime].keys()):
        #     if 'vortex' in key:
        #         xrot = np.append(xrot, locs[filetime][key]['x'])
        #         yrot = np.append(yrot, locs[filetime][key]['y'])
        # if 'vortex3' in locs[filetime].keys():
        #     xrot3 = locs[filetime]['vortex3']['x']
        #     yrot3 = locs[filetime]['vortex3']['y']
        # if 'vortex4' in locs[filetime].keys():
        #     xrot4 = locs[filetime]['vortex4']['x']
        #     yrot4 = locs[filetime]['vortex4']['y']
        
        if i == 4:
            cbflag = True
        else:
            cbflag = False
        
        plot_cfill(vol[vis[i+5]]['xx'][eli,:,:], vol[vis[i+5]]['yy'][eli,:,:], np.ma.masked_array(vol[vis[i+5]]['vel'][eli,:,:], vol[vis[i+5]]['vel'][eli,:,:]<-100), 'vel', ax[1,i],
                   datalims=[-30,30], xlims=xls[i+5], ylims=yls[i+5], cbar=cbflag, cbfs=16, cbts=14, cmap='balance')
        ax[1,i].scatter(0, 0, s=200, c='k', marker='.')
        # ax[1,i].scatter(xrot, yrot, s=40, c='k', marker='.')
        # if 'vortex3' in locs[filetime].keys():
        #     ax[1,i].scatter(xrot3, yrot3, s=30, fc='w', ec='k', marker='o', linewidth=1)
        if xls[0][1] - xls[0][0] == 5:
            ax[1,i].text(xls[i+5][0]+0.2, yls[i+5][-1]-1, f"{panels[i+5]}", color='k', fontsize=28, fontweight='bold')
        elif xls[0][1] - xls[0][0] == 10:
            ax[1,i].text(xls[i+5][0]+0.4, yls[i+5][-1]-2, f"{panels[i+5]}", color='k', fontsize=28, fontweight='bold')
        else:
            ax[1,i].text(xls[i+5][0]+0.05, yls[i+5][-1]-0.3, f"{panels[i+5]}", color='k', fontsize=28, fontweight='bold')
        ax[1,i].set_title(f"{vol[vis[i+5]]['scan_time'][1]} UTC", fontsize=18)
        ax[1,i].set_xlabel('E-W distance (km)', fontsize=16)
    
    ax[0,0].set_ylabel('N-S distance (km)', fontsize=16)
    ax[1,0].set_ylabel('N-S distance (km)', fontsize=16)
    
    if figsave:
        plt.savefig(ip+f"overview_vel_PPI_5x10km.png", dpi=300)




# weird little 3 panel cross section figures
if False:
    ip = "/Users/morgan.schneider/Documents/perils2023/iop2/figs/cross-section-panels/"
    
    xl = [-1, 2]
    yl = [3, 6]
    
    # Vertical pseudovorticity PPI
    fig,ax = plt.subplots(1, 1, figsize=(7.5,6), subplot_kw=dict(box_aspect=1), layout='constrained')
    c = plot_cfill(vol[vi]['xx'][eli,:,:], vol[vi]['yy'][eli,:,:], np.max(vol[vi]['zvort'][1:-2,:,:],axis=0), 'vort', ax,
               datalims=[0,0.25], xlims=xl, ylims=yl, cbfs=12, cbts=12, cmap='LangRainbow12', cbar=False)
    ax.scatter(x_rot, y_rot, s=30, facecolor='w', edgecolor='k', marker='o', linewidth=1)
    # ax.set_title(f"{filetime} UTC Maximum vertical pseudovorticity", fontsize=14)
    ax.set_xlabel('E-W distance from radar (km)', fontsize=14)
    ax.set_ylabel('N-S distance from radar (km)', fontsize=14)
    ax.plot([0,x_rot[0]], [0,y_rot[0]], '-k', linewidth=3)
    ax.plot([0,x_rot[-1]], [0,y_rot[-1]], '-k', linewidth=3)
    ax.plot([0,np.median(x_rot)], [0,np.median(y_rot)], '--k', linewidth=3)
    cb = plt.colorbar(c, ax=ax, extend='max', aspect=30)
    cb.set_label("Pseudovorticity (s$^{-1}$)", fontsize=14)
    # cb.set_ticks(np.arange(0, 0.22, 0.02))
    cb.ax.tick_params(labelsize=12)
    
    if figsave:
        plt.savefig(ip+f"{filetime}_vortPPI.png", dpi=300)
    
    
    # Velocity reconstructed RHI
    fig,ax = plt.subplots(1, 1, figsize=(8,3), layout='constrained')
    plot_cfill(rr, vol[vi]['zz'][:,azi,:], np.ma.masked_array(vol[vi]['vel'][:,azi,:], vol[vi]['dbz'][:,azi,:]<1), 'vel', ax, datalims=[-va,va], xlims=[0,rlim], ylims=[0,zlim], cmap='balance')
    ax.scatter(r_rot[irot], z_rot[irot]+0.03, s=30, c='k', marker='x', linewidth=2)
    ax.set_xlabel('Range from radar (km)', fontsize=12)
    ax.set_ylabel('Height ARL (km)', fontsize=12)
    # ax.set_title(f"{filetime} UTC radial velocity reconstructed RHI", fontsize=12)
    ax.invert_xaxis()
    
    if figsave:
        plt.savefig(ip+f"{filetime}_velRHI.png", dpi=300)
    
    
    # Composite 2D vorticity cross section
    ir1 = np.where(np.isclose(r_rot[0], rr[eli,:], rtol=0.01))[0][0]
    ir2 = np.where(np.isclose(r_rot[-1], rr[eli,:], rtol=0.01))[0][0]
    ia1 = np.where(np.isclose(az_rot[0], vol[vi]['az'][eli,:], atol=0.1))[0][0]
    ia2 = np.where(np.isclose(az_rot[-1], vol[vi]['az'][eli,:], atol=0.1))[0][0]
    
    
    # pulling the actual columns from the data instead of WRF linear cross sections
    vort3d_cross = np.zeros(shape=(len(vol[vi]['zz'][:,0,0]),len(az_rot)), dtype=float)
    for i in range(len(az_rot)):
        ir = np.where(np.isclose(r_rot[i], rr[eli,:], rtol=0.01))[0][0]
        vort3d_cross[:,i] = vol[vi]['vort3d'][:,az_rot[i],ir]
    
    # az_rot2 = az_rot
    az_rot2 = np.array([az_rot[i]+360 if az_rot[i]<az_rot[0] else az_rot[i] for i in range(len(az_rot))])
    if ir2 > ir1:
        zz_rot = np.nanmean(vol[vi]['zz'][:,ia1,slice(ir1,ir2+1)], axis=1)
    else:
        zz_rot = np.nanmean(vol[vi]['zz'][:,ia1,slice(ir2,ir1+1)], axis=1)
    
    rm = np.median(r_rot)
    dl = rm * np.pi/180
    l = rm * (az_rot2[-1] - az_rot2[0])*np.pi/180
    z_max = np.max(z_rot)
    
    az_lims = [az_rot2[0], az_rot2[-1]]
    z_lims = [0, np.round(np.max(z_rot), decimals=1)]
    xtick_pos = np.arange(az_lims[0], az_lims[1]+1, 2)
    xtick_lab = [f"{int(pos % 360)}" for pos in xtick_pos]
    
    
    # fig,ax = plt.subplots(1, 1, figsize=(6,6), subplot_kw=dict(box_aspect=z_max/l), layout='constrained')
    fig,ax = plt.subplots(1, 1, figsize=(6,6), layout='constrained')
    c = plot_cfill(az_rot2, zz_rot, vort3d_cross, 'vort', ax, datalims=[0,0.25], xlims=az_lims, ylims=z_lims,
                   cmap='LangRainbow12', cbfs=12, cbts=12, cbar=False)
    # ax.set_title(f"{filetime} UTC 3-D pseudovorticity cross-section", fontsize=12)
    ax.set_ylabel('Height ARL (km)', fontsize=14)
    ax.set_xlabel('Azimuth (deg)', fontsize=14)
    ax.set_xticks(xtick_pos)
    ax.set_xticklabels(xtick_lab)
    cb = plt.colorbar(c, ax=ax, extend='max', aspect=30)
    cb.set_label("Pseudovorticity (s$^{-1}$)", fontsize=14)
    # cb.set_ticks(np.arange(0, 0.22, 0.02))
    cb.ax.tick_params(labelsize=12)
    
    if figsave:
        plt.savefig(ip+f"{filetime}_vortCS.png", dpi=300)
    
    

# need to add more azimuths to vol 10 , maybe 12, maybe 8
# recheck r_rot for vol 9? and vol 7. and vol 14. 
# use vortex 4 for vol 12
# redo vol 17
# i think vol 13 is okay?
# fix some of the z_rots just for the RHI plots like vol 12 but just do it manually


# vorticity time heights !!!
if False:
    vis = [8, 9, 10, 12, 13, 14, 15, 16]
    r1 = [2.4, 2.4, 2.4, 2.8, 3.3, 3.7, 4.0, 4.5, 5.7]
    r2 = [3.5, 3.5, 3.5, 4.0, 4.3, 4.5, 4.8, 5.2, 6.3]
    az1 = [310, 320, 325, 345, 350, 0, 0, 5, 5]
    az2 = [330, 345, 345, 5, 10, 10, 15, 20, 20]
    zvort_max = np.zeros(shape=(len(vis),13), dtype=float)
    z_max = np.zeros(shape=(len(vis),13), dtype=float)
    vrot_max = np.zeros(shape=(len(vis),13), dtype=float)
    times = []
    
    for i in range(len(vis)):
        filetime = vol[vis[i]]['scan_time'][1]
        r_rot = locs[filetime][f"vortex3"]['r']
        az_rot = locs[filetime][f"vortex3"]['az']
        x_rot = locs[filetime][f"vortex3"]['x']
        y_rot = locs[filetime][f"vortex3"]['y']
        
        ir1 = np.where(np.isclose(r1[i], rr[eli,:], rtol=0.01))[0][0]
        ir2 = np.where(np.isclose(r2[i], rr[eli,:], rtol=0.01))[0][0]
        ia1 = np.where(np.isclose(az1[i], vol[vis[i]]['az'][eli,:], atol=0.1))[0][0]
        ia2 = np.where(np.isclose(az2[i], vol[vis[i]]['az'][eli,:], atol=0.1))[0][0]
        
        ir3 = np.where(np.isclose(np.min(r_rot)-0.07, rr[eli,:], rtol=0.01))[0][0]
        ir4 = np.where(np.isclose(np.max(r_rot)+0.07, rr[eli,:], rtol=0.01))[0][0]
        ia3 = np.where(np.isclose((az_rot[0]-3) % 360, vol[vis[i]]['az'][eli,:], atol=0.1))[0][0]
        ia4 = np.where(np.isclose((az_rot[1]+3) % 360, vol[vis[i]]['az'][eli,:], atol=0.1))[0][0]
        
        if ia2 < ia1:
            zvort_tmp = np.append(vol[vis[i]]['zvort'][:,ia1:,ir1:ir2+1], vol[vis[i]]['zvort'][:,:ia2+1,ir1:ir2+1], axis=1)
            z_tmp = np.append(vol[vis[i]]['zz'][:,ia1:,ir1:ir2+1], vol[vis[i]]['zz'][:,:ia2+1,ir1:ir2+1], axis=1)
        else:
            zvort_tmp = vol[vis[i]]['zvort'][:,ia1:ia2+1,ir1:ir2+1]
            z_tmp = vol[vis[i]]['zz'][:,ia1:ia2+1,ir1:ir2+1]
        
        if ia4 < ia3:
            vel_tmp = np.append(vol[vis[i]]['vel'][:,ia3:,ir3:ir4+1], vol[vis[i]]['vel'][:,:ia4+1,ir3:ir4+1], axis=1)
            xx_tmp = np.append(vol[vis[i]]['xx'][:,ia3:,ir3:ir4+1], vol[vis[i]]['xx'][:,:ia4+1,ir3:ir4+1], axis=1)
            yy_tmp = np.append(vol[vis[i]]['yy'][:,ia3:,ir3:ir4+1], vol[vis[i]]['yy'][:,:ia4+1,ir3:ir4+1], axis=1)
        else:
            vel_tmp = vol[vis[i]]['vel'][:,ia3:ia4+1,ir3:ir4+1]
            xx_tmp = vol[vis[i]]['xx'][:,ia3:ia4+1,ir3:ir4+1]
            yy_tmp = vol[vis[i]]['yy'][:,ia3:ia4+1,ir3:ir4+1]
        
        zvort_tmp = np.ma.masked_array(zvort_tmp, zvort_tmp>0.7)
        
        zvort_max[i,:zvort_tmp.shape[0]] = np.max(zvort_tmp, axis=(1,2))
        # zvort_zmax[i,:z_tmp.shape[0]] = z_tmp[(zvort_tmp == zvort_max)]
        vrot_max[i,:vel_tmp.shape[0]] = (np.max(vel_tmp, axis=(1,2)) - np.min(vel_tmp, axis=(1,2))) / 2
        for j in range(zvort_tmp.shape[0]):
            m = int(vol[vis[i]]['scan_time'][j][2:4])
            s = int(vol[vis[i]]['scan_time'][j][4:6])
            t = datetime(2023,3,3,8,m,s)
            times.append(t)
            iz = np.where(zvort_tmp[j,:,:] == zvort_max[i,j])
            iaz = iz[0][0]
            irz = iz[1][0]
            z_max[i,j] = z_tmp[j,iaz,irz]
            # vrot_max[i,j] = (np.max(vel_tmp[j,:,:]) - np.min(vel_tmp[j,:,:])) / 2
        
        if False:
            x1 = r1[i] * np.sin(np.append(np.arange(270,360),np.arange(0,90))*np.pi/180)
            y1 = r1[i] * np.cos(np.append(np.arange(270,360),np.arange(0,90))*np.pi/180)
            x2 = r2[i] * np.sin(np.append(np.arange(270,360),np.arange(0,90))*np.pi/180)
            y2 = r2[i] * np.cos(np.append(np.arange(270,360),np.arange(0,90))*np.pi/180)
            
            xl = [-4,3]
            yl = [0,7]
            
            # Vertical pseudovorticity PPI
            fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
            plot_cfill(vol[vis[i]]['xx'][eli,:,:], vol[vis[i]]['yy'][eli,:,:], np.max(vol[vis[i]]['zvort'][1:,:,:],axis=0), 'vort', ax,
                       datalims=[0,0.2], xlims=xl, ylims=yl, cbfs=12, cmap='LangRainbow12')
            ax.scatter(x_rot, y_rot, s=30, facecolor='w', edgecolor='k', marker='o', linewidth=1)
            ax.scatter(x1, y1, s=10, color='k', marker='.')
            ax.scatter(x2, y2, s=10, color='k', marker='.')
            ax.set_title(f"{filetime} UTC Maximum vertical pseudovorticity", fontsize=14)
            ax.set_xlabel('E-W distance from radar (km)', fontsize=12)
            ax.set_ylabel('N-S distance from radar (km)', fontsize=12)
            ax.plot([0,x_rot[0]], [0,y_rot[0]], '-k', linewidth=3)
            ax.plot([0,x_rot[-1]], [0,y_rot[-1]], '-k', linewidth=3)
            ax.plot([0,np.median(x_rot)], [0,np.median(y_rot)], '--k', linewidth=3)
            
            plt.show()
        
        if False:
            xl = [-4,3]
            yl = [0,7]
            eli2 = 1
            
            imax = np.where(vel_tmp[eli2,:,:] == np.max(vel_tmp[eli2,:,:]))
            xmax = xx_tmp[eli2,imax[0][0],imax[1][0]]
            ymax = yy_tmp[eli2,imax[0][0],imax[1][0]]
            imin = np.where(vel_tmp[eli2,:,:] == np.min(vel_tmp[eli2,:,:]))
            xmin = xx_tmp[eli2,imin[0][0],imin[1][0]]
            ymin = yy_tmp[eli2,imin[0][0],imin[1][0]]
            
            fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
            plot_cfill(vol[vis[i]]['xx'][eli2,:,:], vol[vis[i]]['yy'][eli2,:,:], vol[vis[i]]['vel'][eli2,:,:], 'vel', ax,
                       datalims=[-30,30], xlims=xl, ylims=yl, cbfs=12, cmap='balance')
            # ax.scatter(x_rot, y_rot, s=30, facecolor='w', edgecolor='k', marker='o', linewidth=1)
            ax.scatter(xmin, ymin, s=30, facecolor='w', edgecolor='k', marker='o', linewidth=1)
            ax.scatter(xmax, ymax, s=30, facecolor='w', edgecolor='k', marker='o', linewidth=1)
            ax.set_title(f"{filetime} UTC {vol[vis[i]]['elev'][eli2]:.1f}$^{{\circ}}$ radial velocity", fontsize=14)
            ax.set_xlabel('E-W distance from radar (km)', fontsize=12)
            ax.set_ylabel('N-S distance from radar (km)', fontsize=12)
            ax.plot([0,x_rot[0]], [0,y_rot[0]], '-k', linewidth=3)
            ax.plot([0,x_rot[-1]], [0,y_rot[-1]], '-k', linewidth=3)
            ax.plot([0,np.median(x_rot)], [0,np.median(y_rot)], '--k', linewidth=3)
            
            plt.show()
    
    
    zvort_max = zvort_max.flatten()
    z_max = z_max.flatten()
    vrot_max = vrot_max.flatten()
    
    zvort_max = zvort_max[(z_max > 0)]
    vrot_max = vrot_max[(z_max > 0)]
    z_max = z_max[(z_max > 0)]
    
    
    fig,ax = plt.subplots(1, 1, figsize=(7,4), layout='constrained')
    c = ax.scatter(times, z_max, c=zvort_max, cmap='HomeyerRainbow', marker='.', s=150, vmin=0, vmax=0.4)
    ax.set_xlabel('Time (UTC)', fontsize=12)
    ax.set_ylabel('Height ARL (km)', fontsize=12)
    ax.set_title('Maximum vertical pseudovorticity', fontsize=14)
    cb = plt.colorbar(c, ax=ax)
    cb.set_label("\u03B6' (s$^{-1}$)", fontsize=12)
    ax.set_ylim([0,1.6])
    ax.set_xlim([datetime(2023,3,3,8,20,0), datetime(2023,3,3,8,25,0)])
    ax.xaxis.set_minor_locator(mdates.SecondLocator(interval=30))
    ax.xaxis.set_minor_locator(mdates.MinuteLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    
    if figsave:
        plt.savefig(ip+f"timeheight_zvort.png", dpi=300)
    
    
    fig,ax = plt.subplots(1, 1, figsize=(7,4), layout='constrained')
    c = ax.scatter(times, z_max, c=vrot_max, cmap='HomeyerRainbow', marker='.', s=150, vmin=2, vmax=10)
    ax.set_xlabel('Time (UTC)', fontsize=12)
    ax.set_ylabel('Height ARL (km)', fontsize=12)
    ax.set_title('Maximum $V_{rot}$', fontsize=14)
    cb = plt.colorbar(c, ax=ax)
    cb.set_label("Velocity (m s$^{-1}$)", fontsize=12)
    ax.set_ylim([0,1.6])
    ax.set_xlim([datetime(2023,3,3,8,20,0), datetime(2023,3,3,8,25,0)])
    ax.xaxis.set_minor_locator(mdates.SecondLocator(interval=30))
    ax.xaxis.set_minor_locator(mdates.MinuteLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    
    if figsave:
        plt.savefig(ip+f"timeheight_vrot.png", dpi=300)
    
    
    
    
#%% Regrid, no advection correction


fp = '/Users/morgan.schneider/Documents/perils2023/iop2/raxpol/CFradial_cal/'
vols = [102, 105, 108, 111, 115, 118, 121] # just the actual analysis volumes - volumes 0820--0823:30
# vols = [102]

### Regrid velocity fields to Cartesian grid ###
meso_x = [-1950.0, -1600.0, -1300.0, -1050.0, -150.0, 300.0, 600.0]
meso_y = [1200.0, 1600.0, 2000.0, 2500.0, 3600.0, 4100.0, 4700.0]

zmax = 2000
delta = 40 #grid spacing in m - 30? 40? 50? i can't tell what i need to do

nx = int(np.ceil(10000/delta))
ny = int(np.ceil(10000/delta))
nz = int(np.ceil(zmax/delta))
# nz = 26

all_files = sorted(glob(fp+'*.nc'))

for vn in range(len(vols)):
    print(f"Volume {vols[vn]}")
    
    # files = sorted(glob(fp+f"cfrad.20230303_*_v{vols[vn]}_*.nc"))
    inds = [i for i,s in enumerate(all_files) if f"v{vols[vn]}" in s]
    if vols[vn] == 111:
        a = slice(inds[1], inds[-1]+1)
    else:
        a = slice(inds[1], inds[-1]+2)
    files = all_files[a]
    
    dbz = np.zeros(shape=(len(files)*360,1246), dtype=float)
    vel = np.zeros(shape=(len(files)*360,1246), dtype=float)
    # vel = np.zeros(shape=(len(files),360,1246), dtype=float)
    elev = np.zeros(shape=(len(files)*360,), dtype=float)
    azim = np.zeros(shape=(len(files)*360,), dtype=float)
    gate_x = np.zeros(shape=(len(files)*360,1246), dtype=float)
    gate_y = np.zeros(shape=(len(files)*360,1246), dtype=float)
    gate_z = np.zeros(shape=(len(files)*360,1246), dtype=float)
    sweep_times = ['' for i in range(len(files))]
    zvort = np.zeros(shape=(len(files)*360,1246,), dtype=float) # vertical pseudovorticity
    # hvort = np.zeros(shape=(len(files)*360,1246,), dtype=float) # cross-radial horizontal pseudovorticity
    # vort3d_tmp = np.zeros(shape=(len(files)*360,1246,), dtype=float) # "3D" (really 2D) pseudovorticity magnitude
    
    xmax,xmin = meso_x[vn]+(delta*nx/2), meso_x[vn]-(delta*nx/2)
    ymax,ymin = meso_y[vn]+(delta*ny/2), meso_y[vn]-(delta*ny/2)
    
    for i in range(len(files)):
        rax = pyart.io.read(files[i])
        dbz_tmp = rax.fields['DBZ']['data']
        vel_tmp = rax.fields['VEL']['data']
        el_tmp = rax.elevation['data']
        x_tmp = rax.gate_x['data']
        y_tmp = rax.gate_y['data']
        z_tmp = rax.gate_altitude['data']
        r = rax.range['data']
        az_tmp = rax.azimuth['data']
        nrays = rax.nrays * len(files)
        time_datetime = datetime.strptime(rax.metadata['time_coverage_start'], "%Y-%m-%dT%H:%M:%SZ")
        sweep_times[i] = time_datetime.strftime("%H%M%S")
        
        for ix in np.linspace(0,359,360):
            if ix in np.floor(rax.azimuth['data']):
                ind = np.where(np.floor(rax.azimuth['data']) == ix)[0][0]
                
                dbz[int(i*360+ix),:] = dbz_tmp[ind,:]
                vel[int(i*360+ix),:] = vel_tmp[ind,:]
                # vel[i,int(ix),:] = vel_tmp[ind,:]
                elev[int(i*360+ix)] = el_tmp[ind]
                azim[int(i*360+ix)] = az_tmp[ind]
                gate_x[int(i*360+ix),:] = x_tmp[ind,:]
                gate_y[int(i*360+ix),:] = y_tmp[ind,:]
                gate_z[int(i*360+ix),:] = z_tmp[ind,:]
                
                if ind == 0:
                    v1 = vel_tmp[-1,:]
                    v2 = vel_tmp[ind+1,:]
                elif ind == len(az_tmp)-1:
                    v1 = vel_tmp[ind-1,:]
                    v2 = vel_tmp[0,:]
                else:
                    v1 = vel_tmp[ind-1,:]
                    v2 = vel_tmp[ind+1,:]
                zvort[int(i*360+ix),:] = 1/r * (v2-v1)/(np.pi/180) # 2/r * dVr/dphi
    
    zvort_metadata = {'long_name':'vertical pseudovorticity', 'standard_name':'inferred vertical pseudovorticity',
                      'units':'1/s', '_FillValue':-32768, 'grid_mapping':'grid_mapping', 'coordinates':'time range',
                      'data':zvort}
    
    dbz = np.ma.masked_array(dbz, dbz < 1)
    dbz[:, :10] = np.ma.masked
    
    center_ind = int(np.round(len(files)/2))
    radar = pyart.io.read(files[center_ind])
    
    radar.gate_x['data'] = gate_x
    radar.gate_y['data'] = gate_y
    radar.gate_altitude['data'] = gate_z
    radar.nrays = 360 * len(files)
    radar.fields['DBZ']['data'] = dbz
    radar.fields['VEL']['data'] = vel
    radar.add_field('ZVORT', zvort_metadata)
    radar.elevation['data'] = elev
    radar.azimuth['data'] = azim
    center_time = sweep_times[center_ind]
    volume_time = sweep_times[1]
    
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_masked('DBZ')
    
    grid = pyart.map.grid_from_radars((radar,), grid_shape=(nz,ny,nx), gatefilters=(gatefilter,),
                                      grid_limits=((0,zmax), (ymin,ymax), (xmin,xmax)), fields=['DBZ','VEL','ZVORT'])
    
    data = {'dbz_grid':grid.fields['DBZ']['data'], 'vel_grid':grid.fields['VEL']['data'], 'zvort_grid':grid.fields['ZVORT']['data'],
            'volume_time':volume_time, 'center_time':center_time, 'sweep_times':sweep_times,
            'x_grid':grid.x['data'], 'y_grid':grid.y['data'], 'z_grid':grid.z['data']}
    
    save_to_pickle(data, f"/Users/morgan.schneider/Documents/perils2023/iop2/raxpol/vel_grid_{volume_time}.pkl", new_pkl=True)




#%% Advection correction and regrid


fp = '/Users/morgan.schneider/Documents/perils2023/iop2/raxpol/CFradial_cal/'
vols = [102, 105, 108, 111, 115, 118, 121] # just the actual analysis volumes - volumes 0820--0823:30
# vols = [102] # for testing/debugging

wspd_storm = [16.0, 17.7, 18.4, 19.8, 21.9, 21.5, 22.0] # storm motion
u_storm = [10.8, 10.6, 10.6, 11.0, 11.7, 11.8, 12.2]
v_storm = [11.7, 13.2, 14.9, 16.4, 18.4, 17.9, 18.2]
meso_x = [-1950.0, -1600.0, -1300.0, -1050.0, -150.0, 300.0, 600.0] #couplet positions at volume center times
meso_y = [1200.0, 1600.0, 2000.0, 2500.0, 3600.0, 4100.0, 4700.0]


zmax = 2000 #upper z limit
delta = 40 #grid spacing in m

# number grid points in cartesian domain
nx = int(np.ceil(10000/delta))
ny = int(np.ceil(10000/delta))
nz = int(np.ceil(zmax/delta))


all_files = sorted(glob(fp+'*.nc'))

# iterate over volumes
for vn in range(len(vols)):
    print(f"Volume {vols[vn]}")
    
    inds = [i for i,s in enumerate(all_files) if f"v{vols[vn]}" in s]
    if vols[vn] == 111: # dumb volume
        a = slice(inds[1], inds[-1]+1)
    else:
        a = slice(inds[1], inds[-1]+2)
    files = all_files[a]
    
    center_ind = np.round(len(files)/2) # index of center sweep in volume
    
    u = u_storm[vn]
    v = v_storm[vn]
    
    x_advected = np.zeros(shape=(len(files)*360,1246), dtype=float)
    y_advected = np.zeros(shape=(len(files)*360,1246), dtype=float)
    
    gate_z = np.zeros(shape=(len(files)*360,1246), dtype=float)
    dbz = np.zeros(shape=(len(files)*360,1246), dtype=float)
    vel = np.zeros(shape=(len(files)*360,1246), dtype=float)
    elev = np.zeros(shape=(len(files)*360,), dtype=float)
    azim = np.zeros(shape=(len(files)*360,), dtype=float)
    sweep_times = ['' for i in range(len(files))]
    zvort = np.zeros(shape=(len(files)*360,1246,), dtype=float) # vertical pseudovorticity
    # hvort = np.zeros(shape=(len(files)*360,1246,), dtype=float) # cross-radial horizontal pseudovorticity
    # vort3d_tmp = np.zeros(shape=(len(files)*360,1246,), dtype=float) # "3D" (really 2D) pseudovorticity magnitude
    
    # cartesian domain limits
    xmax,xmin = meso_x[vn]+(delta*nx/2), meso_x[vn]-(delta*nx/2)
    ymax,ymin = meso_y[vn]+(delta*ny/2), meso_y[vn]-(delta*ny/2)
    
    # iterate over sweeps
    for i in range(len(files)):
        rax = pyart.io.read(files[i])
        gate_x = rax.gate_x['data']
        gate_y = rax.gate_y['data']
        z_tmp = rax.gate_altitude['data']
        dbz_tmp = rax.fields['DBZ']['data']
        vel_tmp = rax.fields['VEL']['data']
        el_tmp = rax.elevation['data']
        az_tmp = rax.azimuth['data']
        r = rax.range['data']
        nrays = rax.nrays * len(files)
        time_datetime = datetime.strptime(rax.metadata['time_coverage_start'], "%Y-%m-%dT%H:%M:%SZ")
        sweep_times[i] = time_datetime.strftime("%H%M%S")
        
        # Advection correction
        print(f"...Advecting sweep {i+1}...")
        dt = 2 * (center_ind - i)
        dx = u * dt # linear translation
        dy = v * dt
        x_new = gate_x + dx
        y_new = gate_y + dy
        
        # Combine all sweeps in the volume into 3d arrays - iterate over azimuth raxpol sucks
        for ix in np.linspace(0,359,360):
            if ix in np.floor(rax.azimuth['data']):
                ind = np.where(np.floor(rax.azimuth['data']) == ix)[0][0]
                
                x_advected[int(i*360+ix),:] = x_new[ind,:]
                y_advected[int(i*360+ix),:] = y_new[ind,:]
                gate_z[int(i*360+ix),:] = z_tmp[ind,:]
                dbz[int(i*360+ix),:] = dbz_tmp[ind,:]
                vel[int(i*360+ix),:] = vel_tmp[ind,:]
                elev[int(i*360+ix)] = el_tmp[ind]
                azim[int(i*360+ix)] = az_tmp[ind]
                
                if ind == 0:
                    v1 = vel_tmp[-1,:]
                    v2 = vel_tmp[ind+1,:]
                elif ind == len(az_tmp)-1:
                    v1 = vel_tmp[ind-1,:]
                    v2 = vel_tmp[0,:]
                else:
                    v1 = vel_tmp[ind-1,:]
                    v2 = vel_tmp[ind+1,:]
                zvort[int(i*360+ix),:] = 1/r * (v2-v1)/(np.pi/180) # 2/r * dVr/dphi
    
    zvort_metadata = {'long_name':'vertical pseudovorticity', 'standard_name':'inferred vertical pseudovorticity',
                      'units':'1/s', '_FillValue':-32768, 'grid_mapping':'grid_mapping', 'coordinates':'time range',
                      'data':zvort}
    
    dbz = np.ma.masked_array(dbz, dbz<1) # mask low SNR
    dbz[:, :10] = np.ma.masked # mask first 10 range gates (300 m)
    
    # Put 3d fields back into radar object for regridding
    radar = pyart.io.read(files[int(center_ind)])
    radar.nrays = 360 * len(files) # have to do this because raxpol is stupid and writes each sweep separately
    radar.gate_x['data'] = x_advected # replace original gate x/y with advection-corrected positions
    radar.gate_y['data'] = y_advected
    radar.gate_altitude['data'] = gate_z
    radar.fields['DBZ']['data'] = dbz # overwrite sweep dbz/vel with the 3d fields
    radar.fields['VEL']['data'] = vel
    radar.add_field('ZVORT', zvort_metadata) # add zvort field to radar object
    radar.elevation['data'] = elev
    radar.azimuth['data'] = azim
    center_time = sweep_times[int(center_ind)]
    volume_time = sweep_times[1] # actual volume start time
    
    # pyart told me to do this idk
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_masked('DBZ')
    
    # Regrid with advection-corrected fields
    print('...Regridding...')
    grid = pyart.map.grid_from_radars((radar,), grid_shape=(nz,ny,nx), gatefilters=(gatefilter,),
                                      grid_limits=((0,zmax), (ymin,ymax), (xmin,xmax)), fields=['DBZ','VEL','ZVORT'])
    
    # save to pickle
    data = {'dbz_grid':grid.fields['DBZ']['data'], 'vel_grid':grid.fields['VEL']['data'], 'zvort_grid':grid.fields['ZVORT']['data'],
            'volume_time':volume_time, 'center_time':center_time, 'sweep_times':sweep_times,
            'x_grid':grid.x['data'], 'y_grid':grid.y['data'], 'z_grid':grid.z['data']}
    
    print('...Saving to pickle...')
    save_to_pickle(data, f"/Users/morgan.schneider/Documents/perils2023/iop2/raxpol/rax_grid_advected_{volume_time}.pkl", new_pkl=True)
    



#%% Reconstructed RHIs and PPIs from advection-corrected gridded data (includes actual dissertation/paper figures)

vol_time = '082130'

figsave = False

fp = '/Users/morgan.schneider/Documents/perils2023/iop2/raxpol/'
ip = '/Users/morgan.schneider/Documents/perils2023/iop2/figs/cross-section-panels/'

# xl = [np.min(xx), np.max(xx)]
# yl = [np.min(yy), np.max(yy)]

dbfile = open(f"/Users/morgan.schneider/Documents/perils2023/iop2/raxpol/vel_grid_{vol_time}.pkl", 'rb')
data = pickle.load(dbfile)
vel = data['vel_grid']
zvort = data['zvort_grid']
xx = data['x_grid']
yy = data['y_grid']
zz = data['z_grid']
sweep_times = data['sweep_times']
center_time = data['center_time']
dbfile.close()


if vol_time == '082030':
    xl = [-3.0, 0.0]
    yl = [0.5, 3.5]
    xlims = [0, 1.0]
    zlims = [0, 1.25]
    r_rot = 2.65
    az_rot = np.arange(314, 332)
    azimuth = 321
elif vol_time == '082100':
    xl = [-2.5, 0.5]
    yl = [0.5, 3.5]
    xlims = [0, 1]
    zlims = [0, 1.25]
    r_rot = 2.7
    az_rot = np.arange(323, 342)
elif vol_time == '082130':
    xl = [-2.5, 0.5]
    yl = [1.0, 4.0]
    xlims = [0, 1]
    zlims = [0, 1.25]
    r_rot = 2.72
    az_rot = np.arange(329, 343)
elif vol_time == '082230':
    xl = [-2.0, 1.0]
    yl = [2.0, 5.0]
    xlims = [0, 1]
    zlims = [0, 1.25]
    r_rot = 3.45
    az_rot = np.arange(351, 363)
    # r_rot = np.linspace(2.95, 3.85, 7)
    # az_rot = np.arange(353, 360)
elif vol_time == '082300':
    xl = [-1.5, 1.5]
    yl = [2.5, 5.5]
    xlims = [0, 1]
    zlims = [0, 1.25]
    r_rot = 3.8
    az_rot = np.array([356, 357, 358, 359, 0, 1, 2, 3, 4, 5, 6, 7])
    # r_rot = np.linspace(3.4, 4.2, 9)
    # az_rot = np.array([357, 358, 359, 0, 1, 2, 3, 4, 5])
elif vol_time == '082330':
    xl = [-1.0, 2.0]
    yl = [2.5, 5.5]
    xlims = [0, 1]
    zlims = [0, 1.25]
    r_rot = 4.2
    az_rot = np.arange(2, 11)
    # r_rot = np.linspace(3.85, 4.5, 5)
    # az_rot = np.arange(3, 8)
elif vol_time == '082000':
    xl = [-3.5, -0.5]
    yl = [0.0, 3.0]
    xlims = [0, 1]
    ylims = [0, 1.25]
    r_rot = 2.65
    az_rot = np.arange(314, 332)
else:
    xl = [np.min(xx), np.max(xx)]
    yl = [np.min(yy), np.max(yy)]

x_rot = r_rot * np.sin(az_rot*np.pi/180)
y_rot = r_rot * np.cos(az_rot*np.pi/180)



dbfile = open(f"/Users/morgan.schneider/Documents/perils2023/iop2/raxpol/rax_grid_advected_{vol_time}.pkl", 'rb')
data = pickle.load(dbfile)
dbz_advected = data['dbz_grid']
vel_advected = data['vel_grid']
zvort_advected = data['zvort_grid']
dbfile.close()


if 'locs' not in locals():
    dbfile = open('/Users/morgan.schneider/Documents/perils2023/iop2/raxpol_vortex_locs.pkl', 'rb')
    locs = pickle.load(dbfile)
    dbfile.close()

# vortex_num = 3
# az_rot = locs[vol_time][f"vortex{vortex_num}"]['az']
# r_rot = locs[vol_time][f"vortex{vortex_num}"]['r']
# z_rot = locs[vol_time][f"vortex{vortex_num}"]['z']
# x_rot = locs[vol_time][f"vortex{vortex_num}"]['x']
# y_rot = locs[vol_time][f"vortex{vortex_num}"]['y']



# for key in list(locs[vol_time].keys()):
#     if 'vortex' in key:
#         x_rot = np.append(x_rot, locs[vol_time][key]['x'])
#         y_rot = np.append(y_rot, locs[vol_time][key]['y'])

# vortex_num = '-all'

# # azimuth = round(np.mean(az_rot))
azimuth = az_rot[round(len(az_rot)/2)]
# azimuth = 322
irot = np.where(np.isclose(az_rot, azimuth))[0][0]


# azi = np.where(vol[vi]['az'][eli,:].round(0) == azimuth)[0][0]
# rr = (vol[vi]['xx'][:,azi,:]**2 + vol[vi]['yy'][:,azi,:]**2)**0.5





xx = xx/1000; yy = yy/1000; zz = zz/1000

az_ray = azimuth
r_ray = 315
x_ray = r_ray * np.sin(az_ray*np.pi/180)
y_ray = r_ray * np.cos(az_ray*np.pi/180)

# RHI cross section of uncorrected gridded velocity
xi1 = np.where(np.abs(xx) == np.min(np.abs(xx)))[0][0]
yi1 = np.where(np.abs(yy) == np.min(np.abs(yy)))[0][0]
xi2 = np.where(np.abs(xx-x_ray) == np.min(np.abs(xx-x_ray)))[0][0]
yi2 = np.where(np.abs(yy-y_ray) == np.min(np.abs(yy-y_ray)))[0][0]

xy = wrf.xy(vel, start_point=(xi1,yi1), end_point=(xi2,yi2))
vel_cs = wrf.interp2dxy(vel, xy) # original velocity field
vel_rhi = vel_cs.data

xy = wrf.xy(vel_advected, start_point=(xi1,yi1), end_point=(xi2,yi2))
vel_adv_cs = wrf.interp2dxy(vel_advected, xy) # advection-corrected velocity field
vel_adv_rhi = vel_adv_cs.data

xy = wrf.xy(zvort_advected, start_point=(xi1,yi1), end_point=(xi2,yi2))
zvort_adv_cs = wrf.interp2dxy(zvort_advected, xy) # advection-corrected velocity field
zvort_adv_rhi = zvort_adv_cs.data

rr = np.linspace(0, r_ray, vel_rhi.shape[1])
R = np.sqrt(xx**2 + yy**2)


# Composite vertical pseudovorticity cross sections
zvort_cross = np.zeros(shape=(len(zz),len(x_rot)), dtype=float)
zvort_adv_cross = np.zeros(shape=(len(zz),len(x_rot)), dtype=float)
x_cross = np.zeros(shape=(len(x_rot),), dtype=float)
for i in range(len(x_rot)):
    ix = np.where(np.abs(xx - x_rot[i]) == np.min(np.abs(xx - x_rot[i])))[0][0]
    iy = np.where(np.abs(yy - y_rot[i]) == np.min(np.abs(yy - y_rot[i])))[0][0]
    zvort_cross[:,i] = zvort[:,iy,ix]
    zvort_adv_cross[:,i] = zvort_advected[:,iy,ix]
    x_cross[i] = np.sqrt((x_rot[i]-x_rot[0])**2 + (y_rot[i]-y_rot[0])**2)

i_rhi = np.where(az_rot == azimuth)[0][0]
x_rhi = x_cross[i_rhi]


# Advection-corrected everything
if False:
    # Max vertical vorticity PPI
    fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
    plot_cfill(xx, yy, np.max(zvort_advected, axis=0), 'vort', ax, datalims=[0,0.15],
               xlims=xl, ylims=yl, cbfs=12, cmap='LangRainbow12')
    ax.scatter(x_rot, y_rot, s=30, marker='o', facecolor='w', edgecolor='k')
    ax.set_xlabel('E-W distance from radar (km)', fontsize=14)
    ax.set_ylabel('N-S distance from radar (km)', fontsize=14)
    ax.plot([0,x_rot[0]], [0,y_rot[0]], '-k', linewidth=3)
    ax.plot([0,x_rot[-1]], [0,y_rot[-1]], '-k', linewidth=3)
    ax.plot([0,x_rot[irot]], [0,y_rot[irot]], '--k', linewidth=3)
    ax.set_xlabel('x distance from radar (km)', fontsize=12)
    ax.set_ylabel('y distance from radar (km)', fontsize=12)
    # ax.set_title(f"{vol_time}z maximum vertical pseudovorticity (advection-corrected)", fontsize=12)
    if figsave:
        plt.savefig(ip+f"{vol_time}_vortPPI_advected.png", dpi=300)
    
    
    # Velocity RHI
    fig,ax = plt.subplots(1, 1, figsize=(9,3), subplot_kw=dict(box_aspect=1/3), layout='constrained')
    plot_cfill(rr, zz, vel_adv_rhi, 'vel', ax, datalims=[-20,20], xlims=[0,r_ray], ylims=[0,2], cbfs=12)
    ax.invert_xaxis()
    ax.set_xlabel('Range from radar (km)', fontsize=12)
    ax.set_ylabel('Height ARL (km)', fontsize=12)
    # ax.set_title(f"{vol_time}z gridded velocity RHI (advection-corrected)", fontsize=12)
    if figsave:
        plt.savefig(ip+f"{vol_time}_velRHI_advected.png", dpi=300)
    
    
    # Vertical vorticity cross section
    fig,ax = plt.subplots(1, 1, figsize=(6.5,6), subplot_kw=dict(box_aspect=1.25), layout='constrained')
    c = plot_cfill(x_cross, zz, zvort_adv_cross, 'vort', ax, datalims=[0,0.15], xlims=xlims, ylims=zlims,
                   cmap='LangRainbow12', cbfs=12, cbts=12, cbar=False)
    # ax.plot([x_rhi, x_rhi], [0, 1.25], '-k', linewidth=1)
    # ax.set_title(f"{vol_time}z vertical pseudovorticity (advection-corrected)", fontsize=12)
    ax.set_ylabel('Height ARL (km)', fontsize=14)
    ax.set_xlabel('Arc length (km)', fontsize=14)
    cb = plt.colorbar(c, ax=ax, extend='max', aspect=30)
    cb.set_label("Vertical pseudovorticity (s$^{-1}$)", fontsize=14)
    cb.ax.tick_params(labelsize=12)
    if figsave:
        plt.savefig(ip+f"{vol_time}_vortCS_advected.png", dpi=300)




# fig,ax = plt.subplots(1, 1, figsize=(9,3), subplot_kw=dict(box_aspect=1/3), layout='constrained')
# plot_cfill(rr, zz, np.ma.masked_array(zvort_adv_rhi, zvort_adv_rhi==0), 'vort', ax, datalims=[0, 0.15], xlims=[0,r_ray], ylims=[0,2], cbfs=12, cmap='LangRainbow12')
# ax.invert_xaxis()
# ax.set_xlabel('Range from radar (km)', fontsize=12)
# ax.set_ylabel('Height ARL (km)', fontsize=12)
# ax.plot([r_rot, r_rot], [0, 1.25], '-k', linewidth=1.25)
# plt.show()

figsave = False


# Uncorrected everything
if False:
    # Max vertical vorticity PPI
    fig,ax = plt.subplots(1, 1, figsize=(8,6), layout='constrained')
    plot_cfill(xx, yy, np.max(zvort, axis=0), 'vort', ax, datalims=[0,0.12],
               xlims=xl, ylims=yl, cbfs=12, cmap='LangRainbow12')
    ax.scatter(x_rot, y_rot, s=30, marker='o', facecolor='w', edgecolor='k')
    ax.set_xlabel('E-W distance from radar (km)', fontsize=14)
    ax.set_ylabel('N-S distance from radar (km)', fontsize=14)
    ax.plot([0,x_rot[0]], [0,y_rot[0]], '-k', linewidth=3)
    ax.plot([0,x_rot[-1]], [0,y_rot[-1]], '-k', linewidth=3)
    ax.plot([0,np.median(x_rot)], [0,np.median(y_rot)], '--k', linewidth=3)
    ax.set_xlabel('x distance from radar (km)', fontsize=12)
    ax.set_ylabel('y distance from radar (km)', fontsize=12)
    ax.set_title(f"{vol_time}z maximum vertical pseudovorticity (no correction)", fontsize=12)
    if figsave:
        plt.savefig(ip+f"{vol_time}_vortPPI_uncorrected.png", dpi=300)
    
    
    # Velocity RHI
    fig,ax = plt.subplots(1, 1, figsize=(10,4), layout='constrained')
    plot_cfill(rr, zz, vel_rhi, 'vel', ax, datalims=[-20,20], xlims=[0,r_ray], ylims=[0,2], cbfs=12)
    ax.invert_xaxis()
    ax.set_xlabel('Range from radar (km)', fontsize=12)
    ax.set_ylabel('Height ARL (km)', fontsize=12)
    ax.set_title(f"{vol_time}z gridded velocity RHI (no correction)", fontsize=12)
    if figsave:
        plt.savefig(ip+f"{vol_time}_velRHI_uncorrected.png", dpi=300)
    
    
    # Vertical vorticity cross section
    xlims = [0, 1]
    zlims = [0, 1.25]
    
    fig,ax = plt.subplots(1, 1, figsize=(7,6), layout='constrained')
    c = plot_cfill(x_cross, zz, zvort_cross, 'vort', ax, datalims=[0,0.12], xlims=xlims, ylims=zlims,
                   cmap='LangRainbow12', cbfs=12, cbts=12, cbar=False)
    ax.set_title(f"{vol_time}z vertical pseudovorticity (no correction)", fontsize=12)
    ax.set_ylabel('Height ARL (km)', fontsize=14)
    ax.set_xlabel('Arc length (km)', fontsize=14)
    # ax.set_xticks(xtick_pos)
    # ax.set_xticklabels(xtick_lab)
    cb = plt.colorbar(c, ax=ax, extend='max', aspect=30)
    cb.set_label("Vertical pseudovorticity (s$^{-1}$)", fontsize=14)
    # cb.set_ticks(np.arange(0, 0.22, 0.02))
    cb.ax.tick_params(labelsize=12)
    if figsave:
        plt.savefig(ip+f"{vol_time}_vortCS_uncorrected.png", dpi=300)



# Gridded velocity reconstructed RHIs (uncorrected and corrected)
if False:
    fig,ax = plt.subplots(1, 1, figsize=(10,4), layout='constrained')
    plot_cfill(rr, zz, vel_rhi, 'vel', ax, datalims=[-20,20], xlims=[0,r_ray], ylims=[0,2], cbfs=12)
    ax.invert_xaxis()
    ax.set_xlabel('Range from radar (km)', fontsize=12)
    ax.set_ylabel('Height ARL (km)', fontsize=12)
    ax.set_title(f"{vol_time}z uncorrected velocity RHI", fontsize=12)
    if figsave:
        plt.savefig(ip+f"cross-section-panels/{vol_time}_velRHI_uncorrected.png", dpi=300)
    
    
    fig,ax = plt.subplots(1, 1, figsize=(10,4), layout='constrained')
    plot_cfill(rr, zz, vel_adv_rhi, 'vel', ax, datalims=[-20,20], xlims=[0,r_ray], ylims=[0,2], cbfs=12)
    ax.invert_xaxis()
    ax.set_xlabel('Range from radar (km)', fontsize=12)
    ax.set_ylabel('Height ARL (km)', fontsize=12)
    ax.set_title(f"{vol_time}z advection-corrected velocity RHI", fontsize=12)
    if figsave:
        plt.savefig(ip+f"cross-section-panels/{vol_time}_velRHI_advected.png", dpi=300)



# Gridded velocity PPIs (uncorrected and corrected)
if False:
    fig,ax = plt.subplots(1, 1, figsize=(8,6), layout='constrained')
    plot_cfill(xx, yy, vel[10,:,:], 'vel', ax, datalims=[-30,30], 
               xlims=[np.min(xx),np.max(xx)], ylims=[np.min(yy),np.max(yy)], cbfs=12)
    ax.set_xlabel('x distance from radar (km)', fontsize=12)
    ax.set_ylabel('y distance from radar (km)', fontsize=12)
    ax.set_title(f"Gridded velocity PPI at {zz[10]*1000:.1f} km ARL (no correction)", fontsize=12)
    
    
    fig,ax = plt.subplots(1, 1, figsize=(8,6), layout='constrained')
    plot_cfill(xx, yy, vel_advected[10,:,:], 'vel', ax, datalims=[-30,30], 
               xlims=[np.min(xx),np.max(xx)], ylims=[np.min(yy),np.max(yy)], cbfs=12)
    ax.set_xlabel('x distance from radar (km)', fontsize=12)
    ax.set_ylabel('y distance from radar (km)', fontsize=12)
    ax.set_title(f"Gridded velocity PPI at {zz[10]*1000:.1f} km ARL (advection-corrected)", fontsize=12)
    
    plt.show()



# Gridded maximum vertical pseudovorticity PPIs (uncorrected and corrected)
if False:
    fig,ax = plt.subplots(1, 1, figsize=(8,6), layout='constrained')
    plot_cfill(xx, yy, np.max(zvort, axis=0), 'vort', ax, datalims=[0,0.15],
               xlims=[np.min(xx),np.max(xx)], ylims=[np.min(yy),np.max(yy)], cbfs=12, cmap='LangRainbow12')
    ax.scatter(x_rot, y_rot, s=30, marker='o', facecolor='w', edgecolor='k')
    ax.set_xlabel('x distance from radar (km)', fontsize=12)
    ax.set_ylabel('y distance from radar (km)', fontsize=12)
    ax.set_title(f"Gridded max zvort PPI (no correction)", fontsize=12)
    
    
    fig,ax = plt.subplots(1, 1, figsize=(8,6), layout='constrained')
    plot_cfill(xx, yy, np.max(zvort_advected, axis=0), 'vort', ax, datalims=[0,0.15],
               xlims=xl, ylims=yl, cbfs=12, cmap='LangRainbow12')
    ax.scatter(x_rot, y_rot, s=30, marker='o', facecolor='w', edgecolor='k')
    ax.set_xlabel('x distance from radar (km)', fontsize=12)
    ax.set_ylabel('y distance from radar (km)', fontsize=12)
    ax.set_title(f"Gridded max zvort PPI (advection-corrected)", fontsize=12)
    
    plt.show()



# Composite vertical pseudovorticity cross section (uncorrected and advection corrected)
if False:
    # Uncorrected
    # fig,ax = plt.subplots(1, 1, figsize=(6,6), subplot_kw=dict(box_aspect=z_max/l), layout='constrained')
    fig,ax = plt.subplots(1, 1, figsize=(7,6), layout='constrained')
    c = plot_cfill(x_cross, zz, zvort_cross, 'vort', ax, datalims=[0,0.12], xlims=xlims, ylims=zlims,
                   cmap='LangRainbow12', cbfs=12, cbts=12, cbar=False)
    ax.set_title(f"{vol_time}z uncorrected vertical pseudovorticity", fontsize=12)
    ax.set_ylabel('Height ARL (km)', fontsize=14)
    ax.set_xlabel('Arc length (km)', fontsize=14)
    # ax.set_xticks(xtick_pos)
    # ax.set_xticklabels(xtick_lab)
    cb = plt.colorbar(c, ax=ax, extend='max', aspect=30)
    cb.set_label("Vertical pseudovorticity (s$^{-1}$)", fontsize=14)
    # cb.set_ticks(np.arange(0, 0.22, 0.02))
    cb.ax.tick_params(labelsize=12)
    # if figsave:
    #     plt.savefig(ip+f"cross-section-panels/{vol_time}_vortCS_uncorrected.png", dpi=300)
    
    
    # Advection corrected
    fig,ax = plt.subplots(1, 1, figsize=(7,6), layout='constrained')
    c = plot_cfill(x_cross, zz, zvort_adv_cross, 'vort', ax, datalims=[0,0.12], xlims=xlims, ylims=zlims,
                   cmap='LangRainbow12', cbfs=12, cbts=12, cbar=False)
    ax.set_title(f"{vol_time}z advection-corrected vertical pseudovorticity", fontsize=12)
    ax.set_ylabel('Height ARL (km)', fontsize=14)
    ax.set_xlabel('Arc length (km)', fontsize=14)
    # ax.set_xticks(xtick_pos)
    # ax.set_xticklabels(xtick_lab)
    cb = plt.colorbar(c, ax=ax, extend='max', aspect=30)
    cb.set_label("Vertical pseudovorticity (s$^{-1}$)", fontsize=14)
    # cb.set_ticks(np.arange(0, 0.22, 0.02))
    cb.ax.tick_params(labelsize=12)
    if figsave:
        plt.savefig(ip+f"cross-section-panels/{vol_time}_vortCS_advected.png", dpi=300)


#%% Worms

fp = '/Users/morgan.schneider/Documents/perils2023/iop2/raxpol/'
ip = '/Users/morgan.schneider/Documents/perils2023/iop2/figs/'

if 'vol' not in locals():
    dbfile = open('/Users/morgan.schneider/Documents/perils2023/iop2/circuit_data.pkl', 'rb')
    tmp = pickle.load(dbfile)
    vol = tmp['Rax']
    dbfile.close()

vi = 10
vol_time = vol[vi]['scan_time'][1]


dbfile = open(f"/Users/morgan.schneider/Documents/perils2023/iop2/raxpol/vel_grid_{vol_time}.pkl", 'rb')
data = pickle.load(dbfile)
vel_grid = data['vel_grid']
zvort_grid = data['zvort_grid']
xx = data['x_grid']/1000
yy = data['y_grid']/1000
zz = data['z_grid']/1000
sweep_times = data['sweep_times']
center_time = data['center_time']
dbfile.close()

dbfile = open(f"/Users/morgan.schneider/Documents/perils2023/iop2/raxpol/rax_grid_advected_{vol_time}.pkl", 'rb')
data = pickle.load(dbfile)
dbz_advected = data['dbz_grid']
vel_advected = data['vel_grid']
zvort_advected = data['zvort_grid']
dbfile.close()

cref = np.max(dbz_advected, axis=0)



figsave = False

# Worms stuff
if True:
    cref = np.max(dbz_advected, axis=0)
    
    xl = [-4, 1]
    yl = [0, 5]
    
    # raw raxpol 0.9 deg zvort
    fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
    c = plot_cfill(vol[vi]['xx'][0,:,:], vol[vi]['yy'][0,:,:], vol[vi]['zvort'][0,:,:], 'vort', ax, datalims=[0,0.2],
               xlims=[-4,1], ylims=[0,5], cbfs=12, cmap=cmocean.cm.dense)
    ax.contour(xx, yy, cref, levels=[40], colors=['k'], linewidths=[1.25])
    ax.contour(xx, yy, np.max(zvort_advected, axis=0), levels=[0.08, 0.10, 0.12],
               colors=['darkorange', 'r', 'maroon'], linewidths=[2])
    ax.set_xlabel('x distance from radar (km)', fontsize=12)
    ax.set_ylabel('y distance from radar (km)', fontsize=12)
    ax.set_title(f"{vol_time} UTC 0.9$^{{\circ}}$ vertical pseudovorticity", fontsize=14)
    
    if figsave:
        plt.savefig(ip+f"worms-{vol_time}_v2.png", dpi=300)
        
        
    # gridded max zvort
    fig,ax = plt.subplots(1, 1, figsize=(8,6), subplot_kw=dict(box_aspect=1), layout='constrained')
    plot_cfill(xx, yy, np.max(zvort_advected, axis=0), 'vort', ax, datalims=[0,0.15],
               xlims=[-4,1], ylims=[0,5], cbfs=12, cmap='LangRainbow12')
    ax.contour(xx, yy, cref, levels=[40], colors=['k'], linewidth=1.25)
    # ax.scatter(x_rot, y_rot, s=30, marker='o', facecolor='w', edgecolor='k')
    ax.set_xlabel('x distance from radar (km)', fontsize=12)
    ax.set_ylabel('y distance from radar (km)', fontsize=12)
    # ax.set_title(f"Gridded max zvort PPI (advection-corrected)", fontsize=12)
    ax.set_title(f"{vol_time} UTC maximum vertical pseudovorticity", fontsize=14)
    
    if figsave:
        plt.savefig(ip+f"worms-max-{vol_time}.png", dpi=300)




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
        az_rot3 = np.array([314, 315, 316, 317, 318, 319,
                           320, 321, 322, 323, 324, 325,
                           326, 327, 328])
        r_rot3 = np.array([2.50, 2.54, 2.58, 2.62, 2.65, 2.68,
                          2.70, 2.72, 2.73, 2.75, 2.76, 2.78,
                          2.79, 2.80, 2.81])
        z_rot3 = np.array([0.12, 0.16, 0.22, 0.28, 0.38, 0.45,
                          0.50, 0.55, 0.62, 0.67, 0.72, 0.76,
                          0.78, 0.82, 0.85])
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
        # az_rot3 = np.array([330, 331, 332, 333, 334, 335, 336, 337])
        # r_rot3 = np.array([2.59, 2.60, 2.61, 2.62, 2.63, 2.64, 2.65, 2.67])
        # z_rot3 = np.array([0.04, 0.06, 0.08, 0.12, 0.16, 0.20, 0.23, 0.27])
        az_rot3 = np.array([332, 333, 334, 335, 336, 337, 338])
        r_rot3 = np.array([2.61, 2.62, 2.63, 2.64, 2.65, 2.67, 2.68])
        z_rot3 = np.array([0.08, 0.12, 0.16, 0.20, 0.23, 0.27, 0.30])
        x_rot3 = r_rot3 * np.sin(az_rot3*np.pi/180)
        y_rot3 = r_rot3 * np.cos(az_rot3*np.pi/180)
        datv3 = {'az':az_rot3, 'r':r_rot3, 'z':z_rot3, 'x':x_rot3, 'y':y_rot3}
        
        # vortex 4
        # az_rot4 = np.array([333, 334, 335, 336, 337, 338, 339])
        # r_rot4 = np.array([2.93, 2.96, 2.99, 3.02, 3.05, 3.08, 3.11])
        # z_rot4 = np.array([0.03, 0.07, 0.11, 0.18, 0.25, 0.27, 0.31])
        az_rot4 = np.array([334, 335, 336, 337, 338])
        r_rot4 = np.array([2.82, 2.84, 2.85, 2.87, 2.89])
        z_rot4 = np.array([0.16, 0.20, 0.24, 0.28, 0.30])
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
        
        datv = {'vortex2':datv2, 'vortex3':datv4, 'vortex4':datv3}
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
c2.plot('reflectivity', 0, vmin=0, vmax=70, cmap='NWSRef')
c2.plot_range_rings([20,40,60,80,100])
ax2.set_xlim([-100,100])
ax2.set_ylim([-100,100])
# ax1.scatter(0, 0, s=50, c='k')
# ax1.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')

ax2 = fig.add_subplot(122)
c2.plot('velocity', vmin=-34, vmax=34, cmap='Carbone42')
c2.plot_range_rings([20,40,60,80,100])
ax2.set_xlim([-100,100])
ax2.set_ylim([-100,100])
ax2.scatter(xr, yr, s=25, c='k')

plt.show()

#%% ZDR correction stuff

fig = plt.figure(figsize=(14,5))

ax1 = fig.add_subplot(121)
c1.plot('DBZ', 0, vmin=0, vmax=70, cmap='NWSRef')
c1.plot_range_rings([5, 10])
ax1.set_xlim([-10,10])
ax1.set_ylim([-10,10])
# ax1.scatter(0, 0, s=50, c='k')
# ax1.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')

ax2 = fig.add_subplot(122)
c2.plot('reflectivity', vmin=0, vmax=70, cmap='NWSRef')
c2.plot_range_rings([20,40,60,80,100])
ax2.set_xlim([-100,100])
ax2.set_ylim([-100,100])
ax2.scatter(xr, yr, s=25, c='k')

plt.show()



fig = plt.figure(figsize=(14,5))

ax1 = fig.add_subplot(121)
c1.plot('DBZ', 0, vmin=0, vmax=70, cmap='NWSRef')
c1.plot_range_rings([5, 10])
ax1.set_xlim([-10,10])
ax1.set_ylim([-10,10])
# ax1.scatter(0, 0, s=50, c='k')
# ax1.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')

ax2 = fig.add_subplot(122)
c2.plot('reflectivity', vmin=0, vmax=70, cmap='NWSRef')
c2.plot_range_rings([20,40,60,80,100])
ax2.set_xlim([xr-10,xr+10])
ax2.set_ylim([yr-10,yr+10])
ax2.scatter(xr, yr, s=25, c='k')

plt.show()



fig = plt.figure(figsize=(14,5))

ax1 = fig.add_subplot(121)
c1.plot('ZDR', 0, vmin=-5, vmax=5, cmap='NWSRef')
c1.plot_range_rings([5, 10])
ax1.set_xlim([-10,10])
ax1.set_ylim([-10,10])
# ax1.scatter(0, 0, s=50, c='k')
# ax1.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')

ax2 = fig.add_subplot(122)
c2.plot('differential_reflectivity', vmin=-5, vmax=5, cmap='NWSRef')
c2.plot_range_rings([20,40,60,80,100])
ax2.set_xlim([-100,100])
ax2.set_ylim([-100,100])
ax2.scatter(xr, yr, s=25, c='k')

plt.show()



fig = plt.figure(figsize=(14,5))

ax1 = fig.add_subplot(121)
c1.plot('ZDR', 0, vmin=-5, vmax=5, cmap='NWSRef')
c1.plot_range_rings([5, 10])
ax1.set_xlim([-10,10])
ax1.set_ylim([-10,10])
# ax1.scatter(0, 0, s=50, c='k')
# ax1.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')

ax2 = fig.add_subplot(122)
c2.plot('differential_reflectivity', vmin=-5, vmax=5, cmap='NWSRef')
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
# c1.plot('DBZ', 0, vmin=0, vmax=70, cmap='NWSRef')
# c1.plot_range_rings([5, 10])
# ax1.set_xlim([-10,10])
# ax1.set_ylim([-10,10])
# # ax1.scatter(0, 0, s=50, c='k')
# # ax1.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')

# ax2 = fig.add_subplot(122)
# c2.plot('reflectivity', vmin=0, vmax=70, cmap='NWSRef')
# c2.plot_range_rings([20,40,60,80,100])
# ax2.set_xlim([-100,100])
# ax2.set_ylim([-100,100])
# ax2.scatter(xr, yr, s=25, c='k')

# plt.show()



fig = plt.figure(figsize=(14,5))

ax1 = fig.add_subplot(121)
c3.plot('DBZ', vmin=0, vmax=70, cmap='NWSRef')
c3.plot_range_rings([5, 10])
ax1.set_xlim([-10,10])
ax1.set_ylim([-10,10])
# ax1.scatter(0, 0, s=50, c='k')
# ax1.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')

ax2 = fig.add_subplot(122)
c2.plot('reflectivity', vmin=0, vmax=70, cmap='NWSRef')
c2.plot_range_rings([20,40,60,80,100])
ax2.set_xlim([xr-10,xr+10])
ax2.set_ylim([yr-10,yr+10])
ax2.scatter(xr, yr, s=25, c='k')

plt.show()



# fig = plt.figure(figsize=(14,5))

# ax1 = fig.add_subplot(121)
# c1.plot('ZDR', 0, vmin=-5, vmax=5, cmap='NWSRef')
# c1.plot_range_rings([5, 10])
# ax1.set_xlim([-10,10])
# ax1.set_ylim([-10,10])
# # ax1.scatter(0, 0, s=50, c='k')
# # ax1.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')

# ax2 = fig.add_subplot(122)
# c2.plot('differential_reflectivity', vmin=-5, vmax=5, cmap='NWSRef')
# c2.plot_range_rings([20,40,60,80,100])
# ax2.set_xlim([-100,100])
# ax2.set_ylim([-100,100])
# ax2.scatter(xr, yr, s=25, c='k')

# plt.show()



fig = plt.figure(figsize=(14,5))

ax1 = fig.add_subplot(121)
c3.plot('ZDR', vmin=-5, vmax=5, cmap='NWSRef')
c3.plot_range_rings([5, 10])
ax1.set_xlim([-10,10])
ax1.set_ylim([-10,10])
# ax1.scatter(0, 0, s=50, c='k')
# ax1.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')

ax2 = fig.add_subplot(122)
c2.plot('differential_reflectivity', vmin=-5, vmax=5, cmap='NWSRef')
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
c2.plot('reflectivity', sweep=0, vmin=0, vmax=70, cmap='NWSRef')
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
c2.plot('differential_reflectivity', sweep=0, vmin=-6, vmax=6, cmap='NWSRef')
c2.plot_range_rings([25,50,75,100])
# ax2.set_xlim([-150,150])
# ax2.set_ylim([-150,150])
ax2.set_xlim([xr-15,xr+15])
ax2.set_ylim([yr-15,yr+15])
ax2.scatter(xr, yr, s=25, c='k')

plt.show()



fig = plt.figure(figsize=(14,5))

ax1 = fig.add_subplot(121)
c1.plot('DBZ', vmin=0, vmax=70, cmap='NWSRef')
c1.plot_range_rings([5,10,15])
ax1.set_xlim([-15,15])
ax1.set_ylim([-15,15])
# ax1.scatter(xr, yr, s=25, c='k')
# ax1.scatter(xr-10, yr, s=10, c='k')
# ax1.scatter(xr, yr+10, s=10, c='k')
# ax1.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')

ax2 = fig.add_subplot(122)
c1.plot('ZDR', vmin=-6, vmax=6, cmap='NWSRef')
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
c2.plot('reflectivity', sweep=0, vmin=0, vmax=70, cmap='NWSRef')
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
c2.plot('differential_reflectivity', sweep=0, vmin=-6, vmax=6, cmap='NWSRef')
c2.plot_range_rings([25,50,75,100])
# ax2.set_xlim([-150,150])
# ax2.set_ylim([-150,150])
ax2.set_xlim([xr-15,xr+15])
ax2.set_ylim([yr-15,yr+15])
ax2.scatter(xr, yr, s=25, c='k')

plt.show()



fig = plt.figure(figsize=(14,5))

ax1 = fig.add_subplot(121)
c3.plot('DBZ', vmin=0, vmax=70, cmap='NWSRef')
c3.plot_range_rings([5,10,15])
ax1.set_xlim([-15,15])
ax1.set_ylim([-15,15])
# ax1.scatter(xr, yr, s=25, c='k')
# ax1.scatter(xr-10, yr, s=10, c='k')
# ax1.scatter(xr, yr+10, s=10, c='k')
# ax1.text(-1, 0.4, 'RaXPol', fontsize=10, fontweight='bold')

ax2 = fig.add_subplot(122)
c3.plot('ZDR', vmin=-6, vmax=6, cmap='NWSRef')
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















