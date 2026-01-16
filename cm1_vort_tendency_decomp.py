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
3) traj_MV1.pkl -- Trajectory data of all parcels in the MV every 5 min
4) traj_clusters_210min_v2.pkl + traj_clusters_220min_v2.pkl -- Indices of parcels in the MV at 210/220 min, clustered by source region
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
ds = nc.Dataset(fp + f"cm1out_{fnum:06d}.nc")
xh = ds.variables['xh'][:].data #Load coordinates
yh = ds.variables['yh'][:].data
zh = ds.variables['z'][:].data
ds.close()




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
dpsi_dt = np.gradient(psi, ptime, axis=0) #time ROC of SR parcel direction




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
if True:
    data = {'time':times, 'tilt_z_sw':t_z_sw, 'tilt_z_cw':t_z_cw,
            'exch_sw_cw':e_sw_cw, 'exch_cw_sw':e_cw_sw} #put all save data into a dict
    dbfile = open(fp + f"vten_tilt_{mvtime}min_parcels.pkl", 'wb') #open new file to save data
    pickle.dump(data, dbfile)
    dbfile.close()

    
del t_z_sw, t_z_cw, e_sw_cw, e_cw_sw








