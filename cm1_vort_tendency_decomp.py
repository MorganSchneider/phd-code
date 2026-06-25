# -*- coding: utf-8 -*-
"""
Created on Wed Jan 14 12:16:31 2026
Python version 3.11
@author: mschne28

Calculate parcel vorticity tilting plan views -- These are parcel-centered composites as in the plan view plots in Figs. 14-16 (as of R1).

This code calculates
1) Decomposition of zeta tilting term into streamwise and crosswise tilting components
2) Exchange of omega_H between streamwise and crosswise directions due to parcel heading changes

--FILES NEEDED--
1) cm1out_000033-000053.nc -- these are on bigbang under /merger/merger-125m
2) (maybe) cm1out_pdata.nc -- also on bigbang under /merger/merger-125m
3) traj_MV1.pkl -- Trajectory data of all parcels in the MV every 5 min. The pickles might all be on bigbang under /merger/merger-125m/pickles?
4) traj_clusters_210min_v2.pkl + traj_clusters_220min_v2.pkl -- Indices of parcels in the MV at 210/220 min clustered by source region
5) storm_motion.pkl -- Estimated storm motion at every model output time (every 1 min) from my MV boxes

Streamwise/crosswise exchange is formulated from Adlerman et al. 1999 Eq. 2-3 / Schenkman et al. 2014 Eq. 1-2
This encompasses the calculations in lines 102-117 and 212-215

******Please let me know if I'm calculating anything wrong, especially since I can't test it myself******

"""

#%% Load modules

import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc
import pickle
from scipy.interpolate import RegularGridInterpolator
from matplotlib.ticker import MultipleLocator

#%% User inputs

####################################
### SET THESE VARIABLES YOURSELF ###
####################################

mvtime = 210 #Analysis time - either 210 or 220 min

fp = '/Volumes/Promise_Pegasus_70TB/merger/merger-125m/' #path to wherever the files are - currently assuming they're all in the same folder


#%% Everything else below this line should be automated unless the files are not all saved in the same place,
#   in which case just change 'fp' in nc.Dataset() and open() to the correct directory as needed

if mvtime == 210:
    fnum = 43 #CM1 output file number
    # Need to have CM1 output files 33-43 -- pull from bigbang
elif mvtime == 220:
    fnum = 53
    # Need to have CM1 output files 43-53 -- pull from bigbang


# Load model grid from output file
# ds = nc.Dataset(fp + f"cm1out_{fnum:06d}.nc")
# xh = ds.variables['xh'][:].data #Load coordinates
# yh = ds.variables['yh'][:].data
# zh = ds.variables['z'][:].data
# ds.close()

dbfile = open(fp + 'coords.pkl', 'rb')
co = pickle.load(dbfile)
xh = co['xh']
yh = co['yh']
zh = co['zh']
dbfile.close()



# Estimated storm motion at every output time (every 1 min) from my MV boxes
dbfile = open(fp + "storm_motion.pkl", 'rb')
sm = pickle.load(dbfile)
u_storm = sm['u_storm']
v_storm = sm['v_storm']
dbfile.close()


# Trajectory data of all parcels in the MV every 5 min
dbfile = open(fp + "traj_MV1.pkl", 'rb')
traj = pickle.load(dbfile)
# Pickle structure
#   traj.keys() = dict_keys(['195min', '200min', '205min', '210min', '215min', '220min', '225min', '230min', '235min', '239min'])
#   traj['210min'].keys() = dict_keys(['pids', 'x', 'y', 'z', 'w', 'zvort', 'b', 'vpg', 'u', 'v'])
dbfile.close()


# Indices of parcels in the MV at a single time by source region
dbfile = open(fp + f"traj_clusters_{mvtime}min_v2.pkl", 'rb')
ccs = pickle.load(dbfile)
cc = ccs['mv1']
# cc = 0 -> Low-level environmental inflow
# cc = 1 -> Mid-level environmental inflow
# cc = 2 -> Supercell outflow
# cc = 3 -> QLCS outflow
dbfile.close()


# Get mid-level parcels (cc = 1)
pids_ml = traj[f"{mvtime}min"]['pids'][(cc == 1)] #indices of ML parcels in the MV at mvtime
x_ml = traj[f"{mvtime}min"]['x'][:,(cc == 1)]/1000 #x position of ML parcels in the MV at mvtime (in km)
y_ml = traj[f"{mvtime}min"]['y'][:,(cc == 1)]/1000
z_ml = traj[f"{mvtime}min"]['z'][:,(cc == 1)]/1000



# Interpolate storm motion from model output freq (60 s) to parcel output freq (15 s)
ptime = np.linspace(10800, 14400, 241) # Parcel output times
mtimes = np.linspace(10800, 14400, 61) # CM1 output times

u_storm_interp = RegularGridInterpolator((mtimes,), u_storm)
v_storm_interp = RegularGridInterpolator((mtimes,), v_storm)
u_storm_prcl = u_storm_interp((ptime,))
v_storm_prcl = v_storm_interp((ptime,))

# Calculate SR parcel direction for SW/CW exchange terms
u_ml = traj[f"{mvtime}min"]['u'][:,(cc == 1)]
v_ml = traj[f"{mvtime}min"]['v'][:,(cc == 1)]
us_ml = u_ml - np.tile(u_storm_prcl, [len(cc[(cc==1)]), 1]).transpose() #SR parcel velocity
vs_ml = v_ml - np.tile(v_storm_prcl, [len(cc[(cc==1)]), 1]).transpose()
psi = np.arctan2(vs_ml, us_ml) #SR parcel direction (from Adlerman and Schenkman papers)
psi2 = np.unwrap(psi, axis=0)
dpsi_dt = np.gradient(psi, ptime, axis=0) #time ROC of SR parcel direction
dpsi_dt = np.gradient(psi2, ptime, axis=0)




times = np.zeros(shape=(11,), dtype=float) #use the 10 minutes leading up to mvtime

# Initialize arrays for tilting terms - specific tilting/exchange components for streamwise and crosswise
e_sw_cw = np.zeros(shape=(len(pids_ml), 11, 33, 33), dtype=float) #exchange from crosswise to streamwise (for d_sw/dt)
e_cw_sw = np.zeros(shape=(len(pids_ml), 11, 33, 33), dtype=float) #exchange from streamwise to crosswise (for d_cw/dt)
t_z_sw = np.zeros(shape=(len(pids_ml), 11, 33, 33), dtype=float) #tilting from streamwise to vertical
t_z_cw = np.zeros(shape=(len(pids_ml), 11, 33, 33), dtype=float) #tilting from crosswise to vertical


# Loop through CM1 output files
m = 0
for fn in np.arange(fnum-10, fnum+1):
    print(f"cm1out_{fn:06d}")
    
    # Open CM1 output files
    ds = nc.Dataset(fp + f"cm1out_{fn:06d}.nc")
    mtime = ds.variables['time'][:].data[0] #model time
    it = np.where(ptime == mtime)[0][0] #where parcel time = model time
    
    # Get max/min x, y, and z positions of all ML parcels at model time + indices of closest grid point
    xmin = np.min(x_ml[it,:]);  ix1 = np.abs(xh - xmin).argmin()
    xmax = np.max(x_ml[it,:]);  ix2 = np.abs(xh - xmax).argmin()
    ymin = np.min(y_ml[it,:]);  iy1 = np.abs(yh - ymin).argmin()
    ymax = np.max(y_ml[it,:]);  iy2 = np.abs(yh - ymax).argmin()
    zmax = np.max(z_ml[it,:]);  iz1 = np.abs(zh - zmax).argmin()
    ix = slice(ix1-20, ix2+21) #2.5 km buffer surrounding all ML parcels at each file time
    iy = slice(iy1-20, iy2+21)
    iz = slice(0, iz1+8) #upper buffer above all ML parcels of ??? km idk it's stretched
    
    # Load model fields only within the indices above (to reduce data size and script runtime)
    u = ds.variables['uinterp'][:].data[0,iz,iy,ix]
    v = ds.variables['vinterp'][:].data[0,iz,iy,ix]
    w = ds.variables['winterp'][:].data[0,iz,iy,ix]
    xvort = ds.variables['xvort'][:].data[0,iz,iy,ix]
    yvort = ds.variables['yvort'][:].data[0,iz,iy,ix]
    
    # calculate SR winds + wind speed on model grid
    u_sr = u - u_storm[fn-13]
    v_sr = v - v_storm[fn-13]
    ws_sr = np.sqrt(u_sr**2 + v_sr**2)
    # calculate streamwise and crosswise vorticity fields
    svort = ((u_sr/ws_sr) * xvort) + ((v_sr/ws_sr) * yvort)
    cvort = ((-v_sr/ws_sr) * xvort) + ((u_sr/ws_sr) * yvort)
    
    # Calculate velocity gradients from model grid
    dwdx = np.gradient(w, xh[ix]*1000, axis=2)
    dwdy = np.gradient(w, yh[iy]*1000, axis=1)
    # dwds = (u_sr/ws_sr)*dwdx + (v_sr/ws_sr)*dwdy #gradient in streamwise direction
    # dwdn = (-v_sr/ws_sr)*dwdx + (u_sr/ws_sr)*dwdy #gradient in crosswise direction
    
    ds.close()
    
    del u, v, w, xvort, yvort #clear variables to free up memory
    
    
    
    # Loop through ML parcels
    for p in range(len(pids_ml)):
        xp = x_ml[it,p] #x position of individual ML parcel at ptime=mtime
        yp = y_ml[it,p]
        zp = z_ml[it,p]
        
        ixp = np.abs(xh[ix]-xp).argmin() #indices of model grid point closest to parcel position
        iyp = np.abs(yh[iy]-yp).argmin()
        k = np.abs(zh[iz]-zp).argmin()
        i = slice(ixp-16,ixp+17) #2 km buffer surrounding parcel
        j = slice(iyp-16,iyp+17)
        
        # Indices-
        # t_* variables are shape (parcel ID, time, y, x) for x-y plan view - indexed [p,m,:,:]
        # all other variables from model grid are shape (z, y, x) - indexed [k,j,i]
        
        # Vertical vorticity tilting components
        t_z_sw[p,m,:,:] = svort[k,j,i] * ((u_sr[k,j,i]/ws_sr[k,j,i]) * dwdx[k,j,i] + 
                                          (v_sr[k,j,i]/ws_sr[k,j,i]) * dwdy[k,j,i])
        t_z_cw[p,m,:,:] = cvort[k,j,i] * ((-v_sr[k,j,i]/ws_sr[k,j,i]) * dwdx[k,j,i] +
                                           (u_sr[k,j,i]/ws_sr[k,j,i]) * dwdy[k,j,i])
        # Exchange from crosswise to streamwise
        e_sw_cw[p,m,:,:] = cvort[k,j,i] * dpsi_dt[it,p]
        # Exchange from streamwise to crosswise
        e_cw_sw[p,m,:,:] = -svort[k,j,i] * dpsi_dt[it,p]
        
        
    times[m] = mtime/60 #save time array as minutes because i can do what i want
    m = m + 1
    
    del u_sr, v_sr, ws_sr, svort, cvort, dwdx, dwdy #clear variables to free up memory


# Save data to pickle file
if False:
    data = {'time':times, 'tilt_z_sw':t_z_sw, 'tilt_z_cw':t_z_cw,
            'exch_sw_cw':e_sw_cw, 'exch_cw_sw':e_cw_sw} #put all save data into a dict
    dbfile = open(fp + f"vten_tilt_{mvtime}min_parcels.pkl", 'wb') #open new file to save data
    pickle.dump(data, dbfile)
    dbfile.close()

    
del t_z_sw, t_z_cw, e_sw_cw, e_cw_sw


#%% Load data

from CM1utils import *



fp = 'C:/Users/mschne28/Documents/merger/merger-125m/'

time = 220

### Choose averaging times ###
if time == 210:
    times = np.arange(206,209)
    # times = np.arange(203,207) #crosswise
    times_cw = np.arange(203,207) #crosswise
elif time == 220:
    times = np.arange(217,220)
    # times = np.arange(213,218) #crosswise
    times_cw = np.arange(213,218) #crosswise


figsave = False


# Load coordinate data
dbfile = open(fp+'coords.pkl', 'rb')
coords = pickle.load(dbfile)
xh = coords['xh']
yh = coords['yh']
zh = coords['zh']
ptime = coords['ptime']
dbfile.close()


# Load filtered parcel data
dbfile = open(fp+'traj_MV1.pkl', 'rb')
traj = pickle.load(dbfile)
dbfile.close()

# Load source regions and pull mid-level source
dbfile = open(fp+f"traj_clusters_{time}min_v2.pkl", 'rb')
tmp = pickle.load(dbfile)
cc = tmp['mv1']
dbfile.close()

z_ml = traj[f"{time}min"]['z'][:,(cc==1)]
z_median = np.median(z_ml, axis=1)
u_ml = traj[f"{time}min"]['u'][:,(cc==1)]
v_ml = traj[f"{time}min"]['v'][:,(cc==1)]


# Load vorticity tendencies
dbfile = open(fp+f"vten_traj_{time}min_parcels.pkl", 'rb')
vten = pickle.load(dbfile)
stime = vten['time']
dbfile.close()

dbfile = open(fp+f"vten_tilt_{time}min_parcels.pkl", 'rb')
vten2 = pickle.load(dbfile)
dbfile.close()

dbfile = open(fp+f"vten_tilt_{time}min_parcels-new.pkl", 'rb')
vten3 = pickle.load(dbfile)
dbfile.close()

# Choose parcel ending layer
dz = 1000
itmv = np.where(ptime == time*60)[0][0]
iz = ((z_ml[itmv,:] >= z_median[itmv]-dz) & (z_ml[itmv,:] <= z_median[itmv]+dz))





it0 = np.where(stime == times[0])[0][0]
itf = np.where(stime == times[-1])[0][0]
it = slice(it0,itf+1)

itc = slice(np.argmin(abs(stime-times_cw[0])), np.argmin(abs(stime-times_cw[-1]))+1)


stretch_x = np.mean(vten['stretch_x'][iz,it,:,:], axis=(0,1))
stretch_y = np.mean(vten['stretch_y'][iz,it,:,:], axis=(0,1))
stretch_z = np.mean(vten['stretch_z'][iz,it,:,:], axis=(0,1))
stretch_h = np.mean(vten['stretch_h'][iz,it,:,:], axis=(0,1))
stretch_sw = np.mean(vten['stretch_sw'][iz,it,:,:], axis=(0,1))
stretch_cw = np.mean(vten['stretch_cw'][iz,itc,:,:], axis=(0,1))

tilt_x = np.mean(vten['tilt_x'][iz,it,:,:], axis=(0,1))
tilt_y = np.mean(vten['tilt_y'][iz,it,:,:], axis=(0,1))
tilt_z = np.mean(vten['tilt_z'][iz,it,:,:], axis=(0,1))
tilt_z_sw = np.mean(vten2['tilt_z_sw'][iz,it,:,:], axis=(0,1))
tilt_z_cw = np.mean(vten2['tilt_z_cw'][iz,it,:,:], axis=(0,1))
tilt_h = np.mean(vten['tilt_h'][iz,it,:,:], axis=(0,1))
tilt_sw = np.mean(vten['tilt_sw'][iz,it,:,:], axis=(0,1))
tilt_cw = np.mean(vten['tilt_cw'][iz,it,:,:], axis=(0,1))
exch_sw_cw = np.mean(vten2['exch_sw_cw'][iz,it,:,:], axis=(0,1))
exch_cw_sw = np.mean(vten2['exch_cw_sw'][iz,itc,:,:], axis=(0,1))
exch3_sw_cw = np.mean(vten3['exch_sw_cw'][iz,it,:,:], axis=(0,1))
exch3_cw_sw = np.mean(vten3['exch_cw_sw'][iz,itc,:,:], axis=(0,1))

bcl_x = np.mean(vten['bcl_x'][iz,it,:,:], axis=(0,1))
bcl_y = np.mean(vten['bcl_y'][iz,it,:,:], axis=(0,1))
bcl_h = np.mean(vten['bcl_h'][iz,it,:,:], axis=(0,1))
bcl_sw = np.mean(vten['bcl_sw'][iz,it,:,:], axis=(0,1))
bcl_cw = np.mean(vten['bcl_cw'][iz,itc,:,:], axis=(0,1))




#% Attempt to correct the exchange terms

esw_old = vten2['exch_sw_cw']
ecw_old = vten2['exch_cw_sw']
esw_new = np.zeros(shape=esw_old.shape)
ecw_new = np.zeros(shape=ecw_old.shape)

# Estimated storm motion at every output time (every 1 min) from my MV boxes
dbfile = open(fp + "storm_motion.pkl", 'rb')
sm = pickle.load(dbfile)
u_storm = sm['u_storm']
v_storm = sm['v_storm']
dbfile.close()


# Interpolate storm motion from model output freq (60 s) to parcel output freq (15 s)
mtimes = np.linspace(10800, 14400, 61) # CM1 output times

u_storm_interp = RegularGridInterpolator((mtimes,), u_storm)
v_storm_interp = RegularGridInterpolator((mtimes,), v_storm)
u_storm_prcl = u_storm_interp((ptime,))
v_storm_prcl = v_storm_interp((ptime,))

# Calculate SR parcel direction for SW/CW exchange terms
us_ml = u_ml - np.tile(u_storm_prcl, [len(cc[(cc==1)]), 1]).transpose() #SR parcel velocity
vs_ml = v_ml - np.tile(v_storm_prcl, [len(cc[(cc==1)]), 1]).transpose()
psi = np.arctan2(vs_ml, us_ml) #SR parcel direction (from Adlerman and Schenkman papers)
psi2 = np.unwrap(psi, axis=0)
dpsi_dt = np.gradient(psi, ptime, axis=0) #time ROC of SR parcel direction
dpsi2_dt = np.gradient(psi2, ptime, axis=0)

for p in range(esw_old.shape[0]):
    for m in range(esw_old.shape[1]):
        itm = np.argmin(abs(ptime/60 - stime[m])) #where parcel time = model time
        esw_new[p,m,:,:] = esw_old[p,m,:,:] / dpsi_dt[itm,p] * dpsi2_dt[itm,p]
        ecw_new[p,m,:,:] = ecw_old[p,m,:,:] / dpsi_dt[itm,p] * dpsi2_dt[itm,p]

exch2_sw_cw = np.mean(esw_new[iz,it,:,:], axis=(0,1))
exch2_cw_sw = np.mean(ecw_new[iz,itc,:,:], axis=(0,1))



#%


# does each term contribute positively or negatively to crosswise vorticity in its current direction (AKA magnitude)
dbfile = open(fp+f"hvort_traj_{time}min.pkl", 'rb')
vort_traj = pickle.load(dbfile)
vort_x = vort_traj['xvort_ml']
vort_y = vort_traj['yvort_ml']
vort_sw = vort_traj['vort_sw_ml']
vort_cw = vort_traj['vort_cw_ml']
dbfile.close()

itp0 = np.where(ptime == times[0]*60)[0][0]
itpf = np.where(ptime == times[-1]*60)[0][0]
itp = slice(itp0,itpf+1)
itpc = slice(np.argmin(abs(ptime-times_cw[0]*60)), np.argmin(abs(ptime-times_cw[-1]*60))+1)
vort_sw = vort_sw[itp,iz]
vort_cw = vort_cw[itpc,iz]


cw_sign = vort_cw / np.abs(vort_cw)

scw = vten['stretch_cw'][iz,itc,:,:]
tcw = vten['tilt_cw'][iz,itc,:,:]
bcw = vten['bcl_cw'][iz,itc,:,:]
# ecw = vten2['exch_cw_sw'][iz,itc,:,:]
ecw = ecw_new[iz,itc,:,:]

for i in range(len(times_cw)):
    for j in range(len(iz[(iz)])):
        scw[j,i,:,:] = scw[j,i,:,:] * cw_sign[i,j]
        tcw[j,i,:,:] = tcw[j,i,:,:] * cw_sign[i,j]
        bcw[j,i,:,:] = bcw[j,i,:,:] * cw_sign[i,j]
        ecw[j,i,:,:] = ecw[j,i,:,:] * cw_sign[i,j]

stretch_cw2 = np.mean(scw, axis=(0,1))
tilt_cw2 = np.mean(tcw, axis=(0,1))
bcl_cw2 = np.mean(bcw, axis=(0,1))
exch_cw2 = np.mean(ecw, axis=(0,1))

xvort = np.mean(np.median(vort_x[itp,iz], axis=0))
yvort = np.mean(np.median(vort_y[itp,iz], axis=0))
xvortc = np.mean(np.median(vort_x[itpc,iz], axis=0))
yvortc = np.mean(np.median(vort_y[itpc,iz], axis=0))

# Load storm motion
dbfile = open(fp+'storm_motion.pkl', 'rb')
sm = pickle.load(dbfile)
u_storm = sm['u_storm']
v_storm = sm['v_storm']
ws_storm = np.sqrt(u_storm**2 + v_storm**2)
dbfile.close()


u_sr = np.mean(np.median(u_ml[itp,iz], axis=0)) - np.mean(u_storm[times[0]-180:times[-1]-179])
v_sr = np.mean(np.median(v_ml[itp,iz], axis=0)) - np.mean(v_storm[times[0]-180:times[-1]-179])
ws_sr = np.mean(ws_storm[times[0]-180:times[-1]-179])

u_src = np.mean(np.median(u_ml[itpc,iz], axis=0)) - np.mean(u_storm[times_cw[0]-180:times_cw[-1]-179])
v_src = np.mean(np.median(v_ml[itpc,iz], axis=0)) - np.mean(v_storm[times_cw[0]-180:times_cw[-1]-179])
ws_src = np.mean(ws_storm[times_cw[0]-180:times_cw[-1]-179])
# swvort = np.mean(np.median(vort_sw, axis=0))
# cwvort = np.mean(np.median(vort_cw, axis=0))
# swvort_x = u_sr/ws_sr * swvort
# swvort_y = v_sr/ws_sr * swvort
# cwvort_x = -v_sr/ws_sr * cwvort
# cwvort_y = u_sr/ws_sr * cwvort

#% Make the plots

from matplotlib.ticker import MultipleLocator

from matplotlib.patches import FancyArrowPatch

def add_vectors(ax, x, y, dx, dy, lengthscale=10, *args, **kwargs):
    # SMALLER lengthscale for SHORTER arrows -- this is REVERSED from scale kwarg in quiver
    for i in range(len(x)):
        x1 = x[i]
        y1 = y[i]
        dx1 = dx[i]*lengthscale
        dy1 = dy[i]*lengthscale
        x2 = x1 + dx1
        y2 = y1 + dy1
        arrow = FancyArrowPatch((x1, y1), (x2, y2), *args, **kwargs)
        ax.add_patch(arrow)
        
    return arrow


# figsave = False


xx = np.linspace(-2, 2, 33)
yy = np.linspace(-2, 2, 33)
lims = [-2.5e-4, 2.5e-4]
# lims = [-5e-4, 5e-4]
levs = np.linspace(lims[0], lims[1], 51)
levs = np.append(np.append([-0.01], levs), [0.01])
xl = [-1, 1]
yl = [-1, 1]

cm = 'balance'



#%% Plot zvort tendency

figsave = False

lims = [-2.5e-4, 2.5e-4]
# lims = [-5e-4, 5e-4]
levs = np.linspace(lims[0], lims[1], 51)
levs = np.append(np.append([-0.01], levs), [0.01])



### Vertical ###
fig,ax = plt.subplots(1, 2, figsize=(7.5,3.5), sharex=True, sharey=True, layout='constrained', subplot_kw=dict(box_aspect=1))

plot_contourf(xx, yy, tilt_z, 'zvort', ax[0], levels=levs, datalims=lims, xlims=xl, ylims=yl, cmap=cm, cbar=False)
ax[0].scatter(0, 0, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[0].set_xlabel('x (km)', fontsize=14)
ax[0].set_ylabel('y (km)', fontsize=14)
# ax[0].text(0.2, 0.8, 'Tilting', color='k', fontsize=18, fontweight='bold')
a_wind = add_vectors(ax[0], [0], [0], [u_sr], [v_sr],
            lengthscale=0.03, arrowstyle='simple', mutation_scale=15, ec='k', fc='k', lw=1)
a_vort = add_vectors(ax[0], [0], [0], [xvort], [yvort],
            lengthscale=20, arrowstyle='simple', mutation_scale=15, ec='k', fc='lightgray', lw=1)
# ax[0].legend(handles=[a_vort,a_wind], labels=["\u03c9$_H$","V$_{SR}$"], loc='lower left', fontsize=10)

c = plot_contourf(xx, yy, stretch_z, 'zvort', ax[1], levels=levs, datalims=lims, xlims=xl, ylims=yl, cmap=cm, cbar=False)
cb = plt.colorbar(c, ax=ax[1], extend='both')
cb.set_label("d\u03B6/dt (s$^{-2}$)", fontsize=12)
cb.set_ticks(np.linspace(-2.5e-4, 2.5e-4, 11))
cb.formatter.set_powerlimits((0,0))
ax[1].scatter(0, 0, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[1].set_xlabel('x (km)', fontsize=14)
# ax[1].text(-0.25, 0.8, 'Stretching', color='k', fontsize=18, fontweight='bold')
a_wind = add_vectors(ax[1], [0], [0], [u_sr], [v_sr],
            lengthscale=0.03, arrowstyle='simple', mutation_scale=15, ec='k', fc='k', lw=1)
a_vort = add_vectors(ax[1], [0], [0], [xvort], [yvort],
            lengthscale=20, arrowstyle='simple', mutation_scale=15, ec='k', fc='lightgray', lw=1)

for i in range(2):
    ax[i].xaxis.set_major_locator(MultipleLocator(0.5))
    ax[i].xaxis.set_minor_locator(MultipleLocator(0.25))
    ax[i].yaxis.set_major_locator(MultipleLocator(0.5))
    ax[i].yaxis.set_minor_locator(MultipleLocator(0.25))
    ax[i].tick_params(axis='both', labelsize=12)

plt.suptitle(f"{times[0]}-{times[-1]} min", fontsize=16)

if figsave:
    plt.savefig(fp+f"figs/new_vort_figs/zvort_composite_{times[0]}-{times[-1]}min.png", dpi=300)





fig,ax = plt.subplots(1, 3, figsize=(9.5,2.98), sharex=True, sharey=True, layout='constrained', subplot_kw=dict(box_aspect=1))

plot_contourf(xx, yy, tilt_z_sw, 'zvort', ax[0], levels=levs, datalims=lims, xlims=xl, ylims=yl, cmap=cm, cbar=False)
ax[0].scatter(0, 0, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[0].set_xlabel('x (km)', fontsize=14)
ax[0].set_ylabel('y (km)', fontsize=14)
# ax[0].text(-0.4, 0.8, "Tilting (SW)", color='k', fontsize=16, fontweight='bold')
a_wind = add_vectors(ax[0], [0], [0], [u_sr], [v_sr],
            lengthscale=0.03, arrowstyle='simple', mutation_scale=15, ec='k', fc='k', lw=1)
a_vort = add_vectors(ax[0], [0], [0], [xvort], [yvort],
            lengthscale=20, arrowstyle='simple', mutation_scale=15, ec='k', fc='lightgray', lw=1)
# ax[0].legend(handles=[a_vort,a_wind], labels=["\u03c9$_H$","V$_{SR}$"], loc='lower left', fontsize=10)

plot_contourf(xx, yy, tilt_z_cw, 'zvort', ax[1], levels=levs, datalims=lims, xlims=xl, ylims=yl, cmap=cm, cbar=False)
ax[1].scatter(0, 0, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[1].set_xlabel('x (km)', fontsize=14)
ax[1].set_title(f"{times[0]}-{times[-1]} min", fontsize=16)
# ax[1].text(-0.4, 0.8, "Tilting (CW)", color='k', fontsize=16, fontweight='bold')
a_wind = add_vectors(ax[1], [0], [0], [u_sr], [v_sr],
            lengthscale=0.03, arrowstyle='simple', mutation_scale=15, ec='k', fc='k', lw=1)
a_vort = add_vectors(ax[1], [0], [0], [xvort], [yvort],
            lengthscale=20, arrowstyle='simple', mutation_scale=15, ec='k', fc='lightgray', lw=1)

c = plot_contourf(xx, yy, stretch_z, 'zvort', ax[2], levels=levs, datalims=lims, xlims=xl, ylims=yl, cmap=cm, cbar=False)
cb = plt.colorbar(c, ax=ax[2], extend='both')
cb.set_label("d\u03B6/dt (s$^{-2}$)", fontsize=12)
cb.set_ticks(np.linspace(-2.5e-4, 2.5e-4, 11))
cb.formatter.set_powerlimits((0,0))
ax[2].scatter(0, 0, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[2].set_xlabel('x (km)', fontsize=14)
# ax[2].text(-0.25, 0.8, "Stretching", color='k', fontsize=16, fontweight='bold')
a_wind = add_vectors(ax[2], [0], [0], [u_sr], [v_sr],
            lengthscale=0.03, arrowstyle='simple', mutation_scale=15, ec='k', fc='k', lw=1)
a_vort = add_vectors(ax[2], [0], [0], [xvort], [yvort],
            lengthscale=20, arrowstyle='simple', mutation_scale=15, ec='k', fc='lightgray', lw=1)

for i in range(3):
    ax[i].xaxis.set_major_locator(MultipleLocator(0.5))
    ax[i].xaxis.set_minor_locator(MultipleLocator(0.25))
    ax[i].yaxis.set_major_locator(MultipleLocator(0.5))
    ax[i].yaxis.set_minor_locator(MultipleLocator(0.25))
    ax[i].tick_params(axis='both', labelsize=12)

if figsave:
    plt.savefig(fp+f"figs/new_vort_figs/zvort_composite_{times[0]}-{times[-1]}min_decomposed.png", dpi=300)


#%% Plot hvort tendency

figsave = False

lims = [-2.5e-4, 2.5e-4]
# lims = [-5e-4, 5e-4]
levs = np.linspace(lims[0], lims[1], 51)
levs = np.append(np.append([-0.01], levs), [0.01])



### Horizontal ###
fig,ax = plt.subplots(1, 3, figsize=(9.5,2.98), sharex=True, sharey=True, layout='constrained', subplot_kw=dict(box_aspect=1))

plot_contourf(xx, yy, tilt_h, 'zvort', ax[0], levels=levs, datalims=lims, xlims=xl, ylims=yl, cmap=cm, cbar=False)
ax[0].scatter(0, 0, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[0].set_xlabel('x (km)', fontsize=14)
ax[0].set_ylabel('y (km)', fontsize=14)
# ax[0].text(0.2, 0.8, 'Tilting', color='k', fontsize=16, fontweight='bold')
a_wind = add_vectors(ax[0], [0], [0], [u_sr], [v_sr],
            lengthscale=0.03, arrowstyle='simple', mutation_scale=15, ec='k', fc='k', lw=1)
a_vort = add_vectors(ax[0], [0], [0], [xvort], [yvort],
            lengthscale=20, arrowstyle='simple', mutation_scale=15, ec='k', fc='lightgray', lw=1)
# ax[0].legend(handles=[a_vort,a_wind], labels=["\u03c9$_H$","V$_{SR}$"], loc='lower left', fontsize=10)

plot_contourf(xx, yy, stretch_h, 'zvort', ax[1], levels=levs, datalims=lims, xlims=xl, ylims=yl, cmap=cm, cbar=False)
ax[1].scatter(0, 0, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[1].set_title(f"{times[0]}-{times[-1]} min", fontsize=16)
ax[1].set_xlabel('x (km)', fontsize=14)
# ax[1].text(-0.25, 0.8, 'Stretching', color='k', fontsize=16, fontweight='bold')
a_wind = add_vectors(ax[1], [0], [0], [u_sr], [v_sr],
            lengthscale=0.03, arrowstyle='simple', mutation_scale=15, ec='k', fc='k', lw=1)
a_vort = add_vectors(ax[1], [0], [0], [xvort], [yvort],
            lengthscale=20, arrowstyle='simple', mutation_scale=15, ec='k', fc='lightgray', lw=1)

c = plot_contourf(xx, yy, bcl_h, 'zvort', ax[2], levels=levs, datalims=lims, xlims=xl, ylims=yl, cmap=cm, cbar=False)
cb = plt.colorbar(c, ax=ax[2], extend='both')
cb.set_label("d\u03c9$_H$/dt (s$^{-2}$)", fontsize=12)
cb.set_ticks(np.linspace(-2.5e-4, 2.5e-4, 11))
cb.formatter.set_powerlimits((0,0))
ax[2].scatter(0, 0, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[2].set_xlabel('x (km)', fontsize=14)
# ax[2].text(-0.2, 0.8, 'Baroclinic', color='k', fontsize=16, fontweight='bold')
a_wind = add_vectors(ax[2], [0], [0], [u_sr], [v_sr],
            lengthscale=0.03, arrowstyle='simple', mutation_scale=15, ec='k', fc='k', lw=1)
a_vort = add_vectors(ax[2], [0], [0], [xvort], [yvort],
            lengthscale=20, arrowstyle='simple', mutation_scale=15, ec='k', fc='lightgray', lw=1)

for i in range(3):
    ax[i].xaxis.set_major_locator(MultipleLocator(0.5))
    ax[i].xaxis.set_minor_locator(MultipleLocator(0.25))
    ax[i].yaxis.set_major_locator(MultipleLocator(0.5))
    ax[i].yaxis.set_minor_locator(MultipleLocator(0.25))
    ax[i].tick_params(axis='both', labelsize=12)

if figsave:
    plt.savefig(fp+f"figs/new_vort_figs/hvort_composite_{times[0]}-{times[-1]}min.png", dpi=300)



#%% Plot streamwise vort tendency

figsave = False

lims = [-1.5e-4, 1.5e-4]
levs = np.linspace(lims[0], lims[1], 41)
levs = np.append(np.append([-0.01], levs), [0.01])



### Streamwise ###

fig,ax = plt.subplots(1, 4, figsize=(12,2.94), sharex=True, sharey=True, layout='constrained', subplot_kw=dict(box_aspect=1))

plot_contourf(xx, yy, exch3_sw_cw, 'zvort', ax[0], levels=levs, datalims=lims, xlims=xl, ylims=yl, cmap=cm, cbar=False)
ax[0].scatter(0, 0, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[0].set_xlabel('x (km)', fontsize=12)
ax[0].set_ylabel('y (km)', fontsize=12)
a_wind = add_vectors(ax[0], [0], [0], [u_sr], [v_sr],
            lengthscale=0.03, arrowstyle='simple', mutation_scale=15, ec='k', fc='k', lw=1)
a_vort = add_vectors(ax[0], [0], [0], [xvort], [yvort],
            lengthscale=20, arrowstyle='simple', mutation_scale=15, ec='k', fc='lightgray', lw=1)
if time == 210:
    ax[0].text(-0.2, 0.8, 'Exchange', color='k', fontsize=16, fontweight='bold',
               bbox=dict(facecolor='w', alpha=0.6, edgecolor='w', boxstyle='Square, pad=0.15'))
    ax[0].legend(handles=[a_vort,a_wind], labels=["\u03c9$_H$","V$_{SR}$"], loc='lower left', fontsize=10)

plot_contourf(xx, yy, tilt_sw, 'zvort', ax[1], levels=levs, datalims=lims, xlims=xl, ylims=yl, cmap=cm, cbar=False)
ax[1].scatter(0, 0, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[1].set_xlabel('x (km)', fontsize=12)
ax[1].set_title(f"                                   {times[0]}-{times[-1]} min", fontsize=18)
a_wind = add_vectors(ax[1], [0], [0], [u_sr], [v_sr],
            lengthscale=0.03, arrowstyle='simple', mutation_scale=15, ec='k', fc='k', lw=1)
a_vort = add_vectors(ax[1], [0], [0], [xvort], [yvort],
            lengthscale=20, arrowstyle='simple', mutation_scale=15, ec='k', fc='lightgray', lw=1)
if time == 210:
    ax[1].text(0.18, 0.8, 'Tilting', color='k', fontsize=16, fontweight='bold',
               bbox=dict(facecolor='w', alpha=0.6, edgecolor='w', boxstyle='Square, pad=0.15'))

plot_contourf(xx, yy, stretch_sw, 'zvort', ax[2], levels=levs, datalims=lims, xlims=xl, ylims=yl, cmap=cm, cbar=False)
ax[2].scatter(0, 0, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
# ax[2].set_title(f"{times[0]}-{times[-1]} min", fontsize=16)
ax[2].set_xlabel('x (km)', fontsize=12)
a_wind = add_vectors(ax[2], [0], [0], [u_sr], [v_sr],
            lengthscale=0.03, arrowstyle='simple', mutation_scale=15, ec='k', fc='k', lw=1)
a_vort = add_vectors(ax[2], [0], [0], [xvort], [yvort],
            lengthscale=20, arrowstyle='simple', mutation_scale=15, ec='k', fc='lightgray', lw=1)
if time == 210:
    ax[2].text(-0.3, 0.8, 'Stretching', color='k', fontsize=16, fontweight='bold',
               bbox=dict(facecolor='w', alpha=0.6, edgecolor='w', boxstyle='Square, pad=0.15'))

c = plot_contourf(xx, yy, bcl_sw, 'zvort', ax[3], levels=levs, datalims=lims, xlims=xl, ylims=yl, cmap=cm, cbar=False)
cb = plt.colorbar(c, ax=ax[3], extend='both')
cb.set_label("d\u03c9$_{{sw}}$/dt (s$^{-2}$)", fontsize=12)
cb.set_ticks(np.linspace(lims[0], lims[1], 11))
cb.formatter.set_powerlimits((0,0))
ax[3].set_xlabel('x (km)', fontsize=12)
ax[3].scatter(0, 0, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
a_wind = add_vectors(ax[3], [0], [0], [u_sr], [v_sr],
            lengthscale=0.03, arrowstyle='simple', mutation_scale=15, ec='k', fc='k', lw=1)
a_vort = add_vectors(ax[3], [0], [0], [xvort], [yvort],
            lengthscale=20, arrowstyle='simple', mutation_scale=15, ec='k', fc='lightgray', lw=1)
if time == 210:
    ax[3].text(-0.22, 0.8, 'Baroclinic', color='k', fontsize=16, fontweight='bold',
               bbox=dict(facecolor='w', alpha=0.6, edgecolor='w', boxstyle='Square, pad=0.15'))

for i in range(4):
    ax[i].xaxis.set_major_locator(MultipleLocator(0.5))
    ax[i].xaxis.set_minor_locator(MultipleLocator(0.25))
    ax[i].yaxis.set_major_locator(MultipleLocator(0.5))
    ax[i].yaxis.set_minor_locator(MultipleLocator(0.25))
    ax[i].tick_params(axis='both', labelsize=12)

if figsave:
    plt.savefig(fp+f"figs/new_vort_figs/swvort_composite_{times[0]}-{times[-1]}min_new.png", dpi=300)




#% Plot crosswise vort tendency

# figsave = False

# lims = [-1e-4, 1e-4]
# levs = np.linspace(lims[0], lims[1], 41)
# levs = np.append(np.append([-0.01], levs), [0.01])



### Crosswise ###

fig,ax = plt.subplots(1, 4, figsize=(12,2.94), sharex=True, sharey=True, layout='constrained', subplot_kw=dict(box_aspect=1))

plot_contourf(xx, yy, exch3_cw_sw, 'zvort', ax[0], levels=levs, datalims=lims, xlims=xl, ylims=yl, cmap=cm, cbar=False)
ax[0].scatter(0, 0, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[0].set_xlabel('x (km)', fontsize=12)
ax[0].set_ylabel('y (km)', fontsize=12)
a_wind = add_vectors(ax[0], [0], [0], [u_src], [v_src],
            lengthscale=0.03, arrowstyle='simple', mutation_scale=15, ec='k', fc='k', lw=1)
a_vort = add_vectors(ax[0], [0], [0], [xvortc], [yvortc],
            lengthscale=20, arrowstyle='simple', mutation_scale=15, ec='k', fc='lightgray', lw=1)
if time == 210:
    ax[0].text(-0.2, 0.8, 'Exchange', color='k', fontsize=16, fontweight='bold',
               bbox=dict(facecolor='w', alpha=0.6, edgecolor='w', boxstyle='Square, pad=0.15'))
    ax[0].legend(handles=[a_vort,a_wind], labels=["\u03c9$_H$","V$_{SR}$"], loc='lower left', fontsize=10)

plot_contourf(xx, yy, tilt_cw, 'zvort', ax[1], levels=levs, datalims=lims, xlims=xl, ylims=yl, cmap=cm, cbar=False)
ax[1].scatter(0, 0, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[1].set_xlabel('x (km)', fontsize=12)
ax[1].set_title(f"                                   {times_cw[0]}-{times_cw[-1]} min", fontsize=18)
a_wind = add_vectors(ax[1], [0], [0], [u_src], [v_src],
            lengthscale=0.03, arrowstyle='simple', mutation_scale=15, ec='k', fc='k', lw=1)
a_vort = add_vectors(ax[1], [0], [0], [xvortc], [yvortc],
            lengthscale=20, arrowstyle='simple', mutation_scale=15, ec='k', fc='lightgray', lw=1)
if time == 210:
    ax[1].text(0.18, 0.8, 'Tilting', color='k', fontsize=16, fontweight='bold',
               bbox=dict(facecolor='w', alpha=0.6, edgecolor='w', boxstyle='Square, pad=0.15'))

plot_contourf(xx, yy, stretch_cw, 'zvort', ax[2], levels=levs, datalims=lims, xlims=xl, ylims=yl, cmap=cm, cbar=False)
ax[2].scatter(0, 0, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
ax[2].set_xlabel('x (km)', fontsize=12)
# ax[2].set_title(f"{times[0]}-{times[-1]} min", fontsize=16)
a_wind = add_vectors(ax[2], [0], [0], [u_src], [v_src],
            lengthscale=0.03, arrowstyle='simple', mutation_scale=15, ec='k', fc='k', lw=1)
a_vort = add_vectors(ax[2], [0], [0], [xvortc], [yvortc],
            lengthscale=20, arrowstyle='simple', mutation_scale=15, ec='k', fc='lightgray', lw=1)
if time == 210:
    ax[2].text(-0.3, 0.8, 'Stretching', color='k', fontsize=16, fontweight='bold',
               bbox=dict(facecolor='w', alpha=0.6, edgecolor='w', boxstyle='Square, pad=0.15'))

c = plot_contourf(xx, yy, bcl_cw, 'zvort', ax[3], levels=levs, datalims=lims, xlims=xl, ylims=yl, cmap=cm, cbar=False)
cb = plt.colorbar(c, ax=ax[3], extend='both')
cb.set_label("d\u03c9$_{{cw}}$/dt (s$^{-2}$)", fontsize=12)
cb.set_ticks(np.linspace(lims[0], lims[1], 11))
cb.formatter.set_powerlimits((0,0))
ax[3].set_xlabel('x (km)', fontsize=12)
ax[3].scatter(0, 0, s=80, edgecolor='k', facecolor='w', marker='o', linewidth=1.5)
a_wind = add_vectors(ax[3], [0], [0], [u_src], [v_src],
            lengthscale=0.03, arrowstyle='simple', mutation_scale=15, ec='k', fc='k', lw=1)
a_vort = add_vectors(ax[3], [0], [0], [xvortc], [yvortc],
            lengthscale=20, arrowstyle='simple', mutation_scale=15, ec='k', fc='lightgray', lw=1)
if time == 210:
    ax[3].text(-0.22, 0.8, 'Baroclinic', color='k', fontsize=16, fontweight='bold',
               bbox=dict(facecolor='w', alpha=0.6, edgecolor='w', boxstyle='Square, pad=0.15'))

for i in range(4):
    ax[i].xaxis.set_major_locator(MultipleLocator(0.5))
    ax[i].xaxis.set_minor_locator(MultipleLocator(0.25))
    ax[i].yaxis.set_major_locator(MultipleLocator(0.5))
    ax[i].yaxis.set_minor_locator(MultipleLocator(0.25))
    ax[i].tick_params(axis='both', labelsize=12)

if figsave:
    plt.savefig(fp+f"figs/new_vort_figs/cwvort_composite_{times_cw[0]}-{times_cw[-1]}min_new.png", dpi=300)





#%% Rewriting all tendency calculations!

from CM1utils import *


mvtime = 210

if mvtime == 210:
    fnum = 43
elif mvtime == 220:
    fnum = 53

calc_zvort = False
calc_xvort = False
calc_yvort = False
calc_hvort = False
calc_swvort = False
calc_cwvort = False
























