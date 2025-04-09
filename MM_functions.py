#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import scipy
import pandas as pd
import numpy as np
from RaxpolUtils import *

#%%

def compute_wdir_bias_rasmussen(dat,R_earth=6371000,analysis_window_s=30):
    #Read in the file and compute acceleration and bearing
    #d = xr.open_dataset(file_name)
    lat,lon = dat.lat.values,dat.lon.values
    wdir,wspd = dat.wdir.values,dat.wspd.values
    time = dat.time.values
    
    irem = np.where((np.isnan(wdir))|(np.isnan(wspd)))[0]
    ikeep = np.setdiff1d(np.arange(lat.size),irem)
    
    lat,lon,time_ = lat[ikeep],lon[ikeep],time[ikeep]
    wdir,wspd = wdir[ikeep],wspd [ikeep]
    
    #Compute heading and speed in 30s moving averages
    dt = 30 #seconds
    
    speed, bearing = compute_speed_and_bearing(lat, lon, dt=dt)
    
    if speed.size == 0:
        return np.nan

    
    acceleration = np.zeros((speed.size))*np.nan
    acceleration[1:-1] = (speed[2:]-speed[:-2]) / 2 #m s-2
    acceleration[0] = speed[1]-speed[0]
    acceleration[-1] = speed[-1]-speed[-2]

    #Compute WDIR bias (Rasmussen 2024)
    lat_deg_N = lat[30:]
    lon_deg_E = lon[30:]
    dir_deg = wdir[30:]
    spd_mPerS = wspd[30:]
    u_ob = -np.sin(np.radians(dir_deg)) * spd_mPerS
    v_ob = -np.cos(np.radians(dir_deg)) * spd_mPerS
    
    accel = acceleration
    hdg_veh_deg = bearing
    spd_veh_mPerS = speed
    u_veh = np.sin(np.radians(hdg_veh_deg)) * spd_veh_mPerS
    v_veh = np.cos(np.radians(hdg_veh_deg)) * spd_veh_mPerS
    time_epoch = time_
    
    cos_hdg_veh = np.cos(np.radians(hdg_veh_deg))
    hdg_cos_change_one_minute = np.abs(np.roll(cos_hdg_veh, -29) - np.roll(cos_hdg_veh, 29))
    
    #print(hdg_veh_deg[1000:1020], spd_veh_mPerS[1000:1020], u_veh[1000:1020], v_veh[1000:1020])
    # print(time_epoch)
    # print(lat_deg_N)
    # print(lon_deg_E)
    # print(dir_deg)
    # print(spd_mPerS)
    
    
    
    # compute vehicle heading and speed across analysis_window_s (30) sec deltas centered on the ob. Future versions of bias correction
    # will need the actual GPS vehicle heading and speed for better accuracy.
    roll_length = int((analysis_window_s - 2) / 2)
    
    # cos_hdg_veh = np.cos(np.radians(hdg_veh_deg))
    # hdg_cos_change_one_minute = np.abs(np.roll(cos_hdg_veh, -29) - np.roll(cos_hdg_veh, 29))
    # plt.plot(epochtime, hdg_cos_change_one_minute) # Create figure and axis objects
    # plt.show() # Display plot to screen
    
    
    # Create dictionaries of the variables needed for computing bias... vehicle u, v, and the RMY propeller speed. Each dictionary
    # entry represents one cluster of proximate obs. To be proximate the obs must occur within a one-minute window with the vehicle
    # moving at two headings nearly orthogonal to each other during that window. Each dictionary entry is a list of values for
    # that cluster
    
    cluster_num = -1
    time_last = 0.0
    cluster_veh_u_dict = {}
    cluster_veh_v_dict = {}
    cluster_veh_RMY_spd_dict = {}
    cluster_min_id = []
    cluster_max_id = []
    
    #print(time_epoch.size, hdg_cos_change_one_minute.size)
    for i in range(time_epoch.size-60):
        if((hdg_cos_change_one_minute[i] > 0.9) &
          (spd_veh_mPerS[i] > 10.0) &
          ((accel[i]*1.0) < 1.0)):
    
            # if it's been > 60 s since the last time we hit the orthogonal route criterion, we assume we're in a new cluster of
            # useful obs. Also if we're at the last ob, we do the if block to close off the last cluster.
            
            if((time_epoch[i] - time_last) > 60.0 or i == (lat_deg_N.size-1)):
                # start new cluster
                
                cluster_num += 1
                time_last = time_epoch[i]
                
                cluster_RMY_spd_list = []
                cluster_veh_u_list = []
                cluster_veh_v_list = []
                cluster_min_id.append(10000000)
                cluster_max_id.append(0)
                
                cluster_veh_u_dict[cluster_num] = cluster_veh_u_list
                cluster_veh_v_dict[cluster_num] = cluster_veh_v_list
                cluster_veh_RMY_spd_dict[cluster_num] = cluster_RMY_spd_list
    
                # if(cluster_num > 0):
                #     print('cluster veh u:',cluster_veh_u_dict[cluster_num-1])
                #     print('cluster veh v:',cluster_veh_v_dict[cluster_num-1])
                #     print('cluster RMY spd:',cluster_veh_RMY_spd_dict[cluster_num-1])
                         
            # print(time_epoch[i], 'veh hdg,spd',hdg_veh_deg[i], spd_veh_mPerS[i], 'brkt hdg',
            #       hdg_veh_deg[i-29], hdg_veh_deg[i+29],dir_deg[i], spd_mPerS[i])
    
            # At this ith ob, we know the +29th and -29th ob (one minute apart) are at greatly disparate road orientations and
            # sufficient vehicle speed to help in our calculation; hence we collect the -29th and 29th
    
            # -29th:
            u_rel = u_veh[i-29] - u_ob[i-29]
            v_rel = v_veh[i-29] - v_ob[i-29]
            rmy_prop_speed_mPerS = np.sqrt((u_rel * u_rel) + (v_rel * v_rel))
            # print('   obs comps',u_ob,v_ob,'rel comps',u_rel,v_rel,'prop spd',rmy_prop_speed_mPerS)
            
            cluster_veh_u_list.append(u_veh[i-29])
            cluster_veh_v_list.append(v_veh[i-29])
            cluster_RMY_spd_list.append(rmy_prop_speed_mPerS)
    
            # +29th
            u_rel = u_veh[i+29] - u_ob[i+29]
            v_rel = v_veh[i+29] - v_ob[i+29]
            rmy_prop_speed_mPerS = np.sqrt((u_rel * u_rel) + (v_rel * v_rel))
            # print('   obs comps',u_ob,v_ob,'rel comps',u_rel,v_rel,'prop spd',rmy_prop_speed_mPerS)
            
            cluster_veh_u_list.append(u_veh[i+29])
            cluster_veh_v_list.append(v_veh[i+29])
            cluster_RMY_spd_list.append(rmy_prop_speed_mPerS)
    
            if(i < cluster_min_id[cluster_num]): cluster_min_id[cluster_num] = i
            if(i > cluster_max_id[cluster_num]): cluster_max_id[cluster_num] = i;
    
    # Now go through the sufficient-sized ob clusters and compute the RMY pointing bias in each one...
    def SFunc(Vveh, utrue, vtrue):
        return (np.sqrt(((Vveh[0] - utrue) * (Vveh[0] - utrue)) + ((Vveh[1] - vtrue) * (Vveh[1] - vtrue))))
    
    sum_bias = 0.0
    count_bias = 0.0
    
    for i_cluster in range(len(cluster_veh_u_dict)):
        if(len(cluster_veh_u_dict[i_cluster]) > 20):
            # print('will process cluster',i_cluster,'with',len(cluster_veh_u_dict[i_cluster]),'obs')
    
            S= np.array(cluster_veh_RMY_spd_dict[i_cluster])
            uveh= np.array(cluster_veh_u_dict[i_cluster])
            vveh= np.array(cluster_veh_v_dict[i_cluster])
            p0 = 0.0, 10.0  # first guess for utrue, vtrue
            try:
                popt, _ = scipy.optimize.curve_fit(SFunc, (uveh, vveh), S, p0)
            except RuntimeError:
                return np.nan
            # print(popt)
    
            # pretty sure the popt values here, utrue and vtrue in the SFunc, are the solved vehicle-wind, best fit for the whole
            # cluster. So to get the ground relative wind we have to add back in the vehicle motion at each observing point.
    
            u_wind_cluster_area = popt[0]
            v_wind_cluster_area = popt[1]
    
            # Compute the RMY pointing angle bias for this cluster
            pointing_angle_bias = []
            for i_ob in range(cluster_min_id[i_cluster], cluster_max_id[i_cluster]):
                
                # 1. Compute actual RMY pointing angle relative to vehicle dead-ahead for this ob
                u_actual_rel = u_veh[i_ob] - u_ob[i_ob]
                v_actual_rel = v_veh[i_ob] - v_ob[i_ob]
                actual_RMY_point_deg = np.degrees(np.arctan2(v_actual_rel, u_actual_rel))
                if(actual_RMY_point_deg < 0.0): actual_RMY_point_deg += 360.0
    
                #2. Compute the expected RMY pointing angle based on the cluster-area wind
                u_expected_rel = u_veh[i_ob] - u_wind_cluster_area
                v_expected_rel = v_veh[i_ob] - v_wind_cluster_area
                expected_RMY_point_deg = np.degrees(np.arctan2(v_expected_rel, u_expected_rel))
                if(expected_RMY_point_deg < 0.0): expected_RMY_point_deg += 360.0
    
                pointing_angle_bias.append(expected_RMY_point_deg - actual_RMY_point_deg)
                sum_bias += (expected_RMY_point_deg - actual_RMY_point_deg)
                count_bias += 1.0
                # print(cluster_min_id[i_cluster],'-',cluster_max_id[i_cluster],'pointing angles actual',actual_RMY_point_deg,'expected',expected_RMY_point_deg,'bias',pointing_angle_bias[-1])
    
            numpy_bias = np.array(pointing_angle_bias)
            median_bias = np.median(numpy_bias)
            #print("uniform comps",popt,"median bias", median_bias,"in",cluster_max_id[i_cluster]-cluster_min_id[i_cluster],"samples")
    
    #print("whole-day mean bias is", sum_bias/count_bias,"from",count_bias,"samples")
    mean_bias = sum_bias/count_bias if count_bias > 0 else sum_bias

    #Correct the wind speed and direction-------------------------------
    RMY_bias_corr_deg = mean_bias
    u_ob_rel = u_veh - u_ob
    v_ob_rel = v_veh - v_ob
    ob_RMY_speed = np.sqrt((u_ob_rel * u_ob_rel) + (v_ob_rel * v_ob_rel))
    
    ob_RMY_point_deg = np.degrees(np.arctan2(v_ob_rel, u_ob_rel))
    corrected_RMY_point_deg = ob_RMY_point_deg + RMY_bias_corr_deg
    u_corrected_rel = np.cos(np.radians(corrected_RMY_point_deg)) * ob_RMY_speed
    v_corrected_rel = np.sin(np.radians(corrected_RMY_point_deg)) * ob_RMY_speed
    
    u_corrected = u_veh - u_corrected_rel
    v_corrected = v_veh - v_corrected_rel
    corrected_dir_deg = 270.0 - np.degrees(np.arctan2(v_corrected, u_corrected))
    corrected_spd_mPerS = np.sqrt((u_corrected*u_corrected) + (v_corrected*v_corrected))

    #Put back the nans to keep the same size----------------------
    acorrected_dir_deg,acorrected_spd_mPerS = np.zeros((dat.lat.values.size))*np.nan,np.zeros((dat.lat.values.size))*np.nan
    au_corrected,av_corrected = np.zeros((dat.lat.values.size))*np.nan,np.zeros((dat.lat.values.size))*np.nan
    
    acorrected_dir_deg[ikeep[analysis_window_s:]] = corrected_dir_deg
    acorrected_spd_mPerS[ikeep[analysis_window_s:]] = corrected_spd_mPerS
    au_corrected[ikeep[analysis_window_s:]] = u_corrected
    av_corrected[ikeep[analysis_window_s:]] = v_corrected
        
    return acorrected_dir_deg,acorrected_spd_mPerS,au_corrected,av_corrected
    
    
def compute_distances(meso_times, meso_lats, meso_lons, point_times, point_lats, point_lons, max_time_diff=60):
    """
    Compute distances from each point to the mesocyclone location at the closest time step.

    Parameters:
    - meso_times (array-like): Timestamps of mesocyclone locations.
    - meso_lats (array-like): Latitudes of mesocyclone locations.
    - meso_lons (array-like): Longitudes of mesocyclone locations.
    - point_times (array-like): Timestamps of points.
    - point_lats (array-like): Latitudes of points.
    - point_lons (array-like): Longitudes of points.
    - max_time_diff (int): Maximum allowed time difference in seconds.

    Returns:
    - distances (array): Distances (km) between mesocyclone and points.
    - matched_times (array): Matched mesocyclone times for each point.
    """

    distances = []
    matched_times = []

    # Convert times to pandas datetime for easy matching
    meso_times = pd.to_datetime(meso_times)
    point_times = pd.to_datetime(point_times)

    for p_time, p_lat, p_lon in zip(point_times, point_lats, point_lons):
        # Find the mesocyclone time closest to the point time
        time_diffs = np.abs((meso_times - p_time).total_seconds())
        closest_idx = np.argmin(time_diffs)
        
        # Check if the time difference is within the allowed window
        if time_diffs[closest_idx] <= max_time_diff:
            matched_times.append(meso_times[closest_idx])
            
            # Compute distance
            meso_loc = (meso_lats[closest_idx], meso_lons[closest_idx])
            point_loc = (p_lat, p_lon)
            # distance_km = get_dx_dy(meso_lons[closest_idx],meso_lats[closest_idx],p_lon,p_lat)
            distance_km = latlon2xy(p_lat,p_lon,meso_lats[closest_idx],meso_lons[closest_idx])
            distances.append(distance_km)
        else:
            # If no close time match, store NaN
            matched_times.append(np.nan)
            distances.append([np.nan,np.nan])

    return np.array(distances)
