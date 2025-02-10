# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 12:32:34 2023

@author: morgan.schneider
"""

####################
### Load modules ###
####################

from CM1utils import *
from metpy.plots import SkewT, Hodograph
import metpy.calc as mc
from metpy.units import units

#%%

# Effective Shear Algorithm
def effective_layer(p, T, Td, h, height_layer=False):
    '''This function determines the effective inflow layer for a convective sounding.
    
    Input:
      - p: sounding pressure with units
      - T: sounding temperature with units
      - Td: sounding dewpoint temperature with units
      - h: sounding heights with units
      
    Returns:
      - pbot/hbot, ptop/htop: pressure/height of the bottom level, pressure/height of the top level
    '''
    
    pbot = None
    
    for i in range(p.shape[0]):
        prof = mc.parcel_profile(p[i:], T[i], Td[i])
        sbcape, sbcin = mc.cape_cin(p[i:], T[i:], Td[i:], prof)
        if sbcape >= 100 * units('J/kg') and sbcin > -250 * units('J/kg'):
            pbot = p[i]
            hbot = h[i]
            bot_idx = i
            break
    if not pbot:
        return None, None
    
    for i in range(bot_idx+1, p.shape[0]):
        prof = mc.parcel_profile(p[i:], T[i], Td[i])
        sbcape, sbcin = mc.cape_cin(p[i:], T[i:], Td[i:], prof)
        if sbcape < 100 * units('J/kg') or sbcin < -250 * units('J/kg'):
            ptop = p[i]
            htop = h[i]
            break
            
    if height_layer:
        return hbot, htop
    else:
        return pbot, ptop

#%% Load obs data

save_flag = False
iop = 2
snd_time = 1959

snd_fp = f"/Users/morgan.schneider/Documents/PERiLS_LIDAR_Soundings/IOP{iop}/"
snd_fn = glob(snd_fp+f"*_{snd_time}*.csv")[0]

df = pd.read_csv(snd_fn,header=2)

sndst = np.array(list(df['Filtered Temperature (K)']))
sndtd = np.array(list(df['Filtered Dewpoint (K)']))
sndsp = np.array(list(df['Filtered Pressure (mb)']))
sndwd = np.array(list(df['Filtered Wind Dir']))
sndws = np.array(list(df['Filtered Wind Spd (m/s)']))
sndsh = np.array(list(df['Filtered Altitude (m)']))
sndqv = np.array(list(df['Calculated Mixing Ratio']))
sndth = sndst*(1000./sndsp)**0.286

wnd = mc.wind_components(sndws*units('m/s'), sndwd*units.deg)
sndsu = wnd[0].magnitude[:]
sndsv = wnd[1].magnitude[:]

#%% load cm1 data

#fp = f"/Users/morgan.schneider/Documents/cm1/PERiLS/IOP{iop}/{snd_time}/dry/"
#ff = glob(fp+'cm1out_0000*.nc')
#fn = ff[-1]
fp = '/Volumes/Promise_Pegasus_70TB/merger/merger-125m/base/'
fn = fp+'cm1out_000001.nc'
#ds = read_cm1out(fn,['z','th','qv','prs','u0','v0','th0','qv0','prs0'])

ds = nc.Dataset(fn)

z = ds.variables['z'][:].data # km
u0 = ds.variables['u0'][:].data[0,:,0,0] # m/s
v0 = ds.variables['v0'][:].data[0,:,0,0] # m/s
th0 = ds.variables['th0'][:].data[0,:,0,0] # K
qv0 = ds.variables['qv0'][:].data[0,:,0,0] # kg/kg
prs0 = ds.variables['prs0'][:].data[0,:,0,0] # Pa

ds.close()



T0 = th0 * (prs0/100000.)**0.286
e0 = (qv0 * prs0/100) / (0.622+qv0)
Td0 = 243.5 / ((17.67/(np.log(e0/6.112)))-1) + 273.15



#%% obs sounding

# 1 m/s = 1.944 kt

eff_snd = effective_layer(sndsp*units.hPa, sndst*units.K, sndtd*units.K, sndsh*units.m, height_layer=True)
ebot_snd = eff_snd[0].magnitude
etop_snd = eff_snd[1].magnitude

bwnd = mc.bunkers_storm_motion(sndsp*units.hPa, sndsu*units('m/s'), sndsv*units('m/s'), sndsh*units.m)
uBR_snd = bwnd[0].magnitude[0]
vBR_snd = bwnd[0].magnitude[1]
uBL_snd = bwnd[1].magnitude[0]
vBL_snd = bwnd[1].magnitude[1]
u06_snd = bwnd[2].magnitude[0]
v06_snd = bwnd[2].magnitude[1]
smBR_snd = np.sqrt(uBR_snd**2 + vBR_snd**2)
angBR_snd = 180 + np.arctan2(vBR_snd, uBR_snd)*180/np.pi
VH06_snd = np.sqrt(u06_snd**2 + v06_snd**2) #0-6km mean wind
ang06_snd = 180 + np.arctan2(v06_snd, u06_snd)*180/np.pi

T_parcel_snd = mc.parcel_profile(sndsp*units.hPa, sndst[0]*units.K, sndtd[0]*units.K)
CC_snd = mc.cape_cin(sndsp*units.hPa, sndst*units.K, sndtd*units.K, T_parcel_snd)
cape_snd = CC_snd[0].magnitude
cin_snd = CC_snd[1].magnitude
MU_snd = mc.most_unstable_cape_cin(sndsp*units.hPa, sndst*units.K, sndtd*units.K, height=sndsh*units.m)
mucape_snd = MU_snd[0].magnitude
mucin_snd = MU_snd[1].magnitude
SB_snd = mc.surface_based_cape_cin(sndsp*units.hPa, sndst*units.K, sndtd*units.K)
sbcape_snd = SB_snd[0].magnitude
sbcin_snd = SB_snd[1].magnitude
ML_snd = mc.mixed_layer_cape_cin(sndsp*units.hPa, sndst*units.K, sndtd*units.K, height=sndsh*units.m)
mlcape_snd = ML_snd[0].magnitude
mlcin_snd = ML_snd[1].magnitude

lcl_snd = mc.lcl(sndsp[0]*units.hPa, sndst[0]*units.K, sndtd[0]*units.K)
Plcl_snd = lcl_snd[0].magnitude
ilcl = np.where(sndsp <= Plcl_snd)[0][0]
Zlcl_snd = sndsh[ilcl]

el_snd = mc.el(sndsp*units.hPa, sndst*units.K, sndtd*units.K)
Pel_snd = el_snd[0].magnitude
iel = np.where(sndsp <= Pel_snd)[0][0]
Zel_snd = sndsh[iel]

shr6km_snd = mc.bulk_shear(sndsp*units.hPa, sndsu*units('m/s'), sndsv*units('m/s'), height=sndsh*units.m, bottom=min(sndsh)*units.m, depth=6000*units.m)
ushr06_snd = shr6km_snd[0].magnitude
vshr06_snd = shr6km_snd[1].magnitude
shr06_snd = np.sqrt(ushr06_snd**2 + vshr06_snd**2)
angSHR_snd = 180 + np.arctan2(vshr06_snd, ushr06_snd)*180/np.pi
shrE_snd = mc.bulk_shear(sndsp*units.hPa, sndsu*units('m/s'), sndsv*units('m/s'), height=sndsh*units.m, bottom=ebot_snd*units.m, depth=(Zel_snd-ebot_snd)*units.m)
uEBS_snd = shrE_snd[0].magnitude
vEBS_snd = shrE_snd[1].magnitude
EBS_snd = np.sqrt(uEBS_snd**2 + vEBS_snd**2)
angEBS_snd = 180 + np.arctan2(vEBS_snd, uEBS_snd)*180/np.pi

srh500_snd = mc.storm_relative_helicity(sndsh*units.m, sndsu*units('m/s'), sndsv*units('m/s'), depth=500*units.m, storm_u=uBR_snd*units('m/s'), storm_v=vBR_snd*units('m/s'))
SRH500_snd = srh500_snd[2].magnitude
srh1km_snd = mc.storm_relative_helicity(sndsh*units.m, sndsu*units('m/s'), sndsv*units('m/s'), depth=1000*units.m, storm_u=uBR_snd*units('m/s'), storm_v=vBR_snd*units('m/s'))
SRH1km_snd = srh1km_snd[2].magnitude
srh3km_snd = mc.storm_relative_helicity(sndsh*units.m, sndsu*units('m/s'), sndsv*units('m/s'), depth=3000*units.m, storm_u=uBR_snd*units('m/s'), storm_v=vBR_snd*units('m/s'))
SRH3km_snd = srh3km_snd[2].magnitude
esrh_snd = mc.storm_relative_helicity(sndsh*units.m, sndsu*units('m/s'), sndsv*units('m/s'), depth=(etop_snd-ebot_snd)*units.m, bottom=ebot_snd*units.m, storm_u=uBR_snd*units('m/s'), storm_v=vBR_snd*units('m/s'))
ESRH_snd = esrh_snd[2].magnitude

crit_snd = mc.critical_angle(sndsp*units.hPa, sndsu*units('m/s'), sndsv*units('m/s'), sndsh*units.m, uBR_snd*units('m/s'), vBR_snd*units('m/s'))
crit_snd = crit_snd.magnitude
scp_snd = mc.supercell_composite(mucape_snd*units('J/kg'), ESRH_snd*units('m^2/s^2'), EBS_snd*units('m/s'))
SCP_snd = scp_snd[0].magnitude
stp_snd = mc.significant_tornado(sbcape_snd*units('J/kg'), Zlcl_snd*units.m, SRH500_snd*units('m^2/s^2'), shr06_snd*units('m/s'))
STP_snd = stp_snd[0].magnitude

#%
    

print("---RAW SOUNDING---")
print(f"Profile top:          {max(sndsh):.0f} m")
print(f"Sfc CAPE,CIN:         {cape_snd:.0f} J/kg, {cin_snd:.0f} J/kg")
print(f"MUCAPE,MUCIN:         {mucape_snd:.0f} J/kg, {mucin_snd:.0f} J/kg")
print(f"SBCAPE,SBCIN:         {sbcape_snd:.0f} J/kg, {sbcin_snd:.0f} J/kg")
print(f"MLCAPE,MLCIN:         {mlcape_snd:.0f} J/kg, {mlcin_snd:.0f} J/kg")
print(f"LCL height:           {Zlcl_snd:.0f} m ({Plcl_snd:.0f} hPa)")
print(f"Eff. inflow layer:    {ebot_snd:.0f} m-{etop_snd:.0f} m")
print(f"Bunkers RM:           {smBR_snd:.1f} m/s at {angBR_snd:.0f} deg (Vector: {uBR_snd:.1f} m/s, {vBR_snd:.1f} m/s)")
print(f"0-6 km mean wind:     {VH06_snd:.1f} m/s at {ang06_snd:.0f} deg (Vector: {u06_snd:.1f} m/s, {v06_snd:.1f} m/s)")
print(f"-----")
print(f"Bulk shear (0-6 km):  {shr06_snd:.1f} m/s at {angSHR_snd:.0f} deg (Vector: {ushr06_snd:.1f} m/s, {vshr06_snd:.1f} m/s)")
print(f"Bulk shear (Eff.):    {EBS_snd:.1f} m/s at {angEBS_snd:.0f} deg (Vector: {uEBS_snd:.1f} m/s, {vEBS_snd:.1f} m/s)")
print(f"SRH (0-500 m):        {SRH500_snd:.0f} m2/s2")
print(f"SRH (0-1 km):         {SRH1km_snd:.0f} m2/s2")
print(f"SRH (0-3 km):         {SRH3km_snd:.0f} m2/s2")
print(f"SRH (Eff.):           {ESRH_snd:.0f} m2/s2")
print(f"Critical angle:       {crit_snd:.1f} degrees")
print(f"Supercell composite:  {SCP_snd:.1f}")
print(f"Significant tornado:  {STP_snd:.1f}")


print(f"-----------------------------------")
#%% cm1 stable sounding

eff_cm1 = effective_layer(prs0/100*units.hPa, T0*units.K, Td0*units.K, z*1000*units.m, height_layer=True)
ebot = eff_cm1[0].magnitude
etop = eff_cm1[1].magnitude

bwnd = mc.bunkers_storm_motion(prs0/100*units.hPa, u0*units('m/s'), v0*units('m/s'), z*1000*units.m)
uBR = bwnd[0].magnitude[0]
vBR = bwnd[0].magnitude[1]
uBL = bwnd[1].magnitude[0]
vBL = bwnd[1].magnitude[1]
u06 = bwnd[2].magnitude[0]
v06 = bwnd[2].magnitude[1]
smBR = np.sqrt(uBR**2 + vBR**2)
angBR = 180 + np.arctan2(vBR, uBR)*180/np.pi
VH06 = np.sqrt(u06**2 + v06**2)
ang06 = 180 + np.arctan2(v06, u06)*180/np.pi

T_parcel = mc.parcel_profile(prs0/100*units.hPa, T0[0]*units.K, Td0[0]*units.K)
CC_cm1 = mc.cape_cin(prs0/100*units.hPa, T0*units.K, Td0*units.K, T_parcel)
cape = CC_cm1[0].magnitude
cin = CC_cm1[1].magnitude
MU_cm1 = mc.most_unstable_cape_cin(prs0/100*units.hPa, T0*units.K, Td0*units.K, height=z*1000*units.m)
mucape = MU_cm1[0].magnitude
mucin = MU_cm1[1].magnitude
SB_cm1 = mc.surface_based_cape_cin(prs0/100*units.hPa, T0*units.K, Td0*units.K)
sbcape = SB_cm1[0].magnitude
sbcin = SB_cm1[1].magnitude
ML_cm1 = mc.mixed_layer_cape_cin(prs0/100*units.hPa, T0*units.K, Td0*units.K, height=z*1000*units.m)
mlcape = ML_cm1[0].magnitude
mlcin = ML_cm1[1].magnitude

lcl_cm1 = mc.lcl(prs0/100*units.hPa, T0[0]*units.K, Td0[0]*units.K)
Plcl = lcl_cm1[0].magnitude[0]
ilcl = np.where(prs0/100 <= Plcl)[0][0]
Zlcl = z[ilcl]*1000

el_cm1 = mc.el(prs0/100*units.hPa, T0*units.K, Td0*units.K)
Pel = el_cm1[0].magnitude
iel = np.where(prs0/100 <= Pel)[0][0]
Zel = z[iel]*1000

shr6km = mc.bulk_shear(prs0/100*units.hPa, u0*units('m/s'), v0*units('m/s'), height=z*1000*units.m, bottom=10*units.m, depth=6000*units.m)
ushr06 = shr6km[0].magnitude
vshr06 = shr6km[1].magnitude
shr06 = np.sqrt(ushr06**2 + vshr06**2)
angSHR = 180 + np.arctan2(vshr06, ushr06)*180/np.pi
shrE = mc.bulk_shear(prs0/100*units.hPa, u0*units('m/s'), v0*units('m/s'), height=z*1000*units.m, bottom=ebot*units.m, depth=(Zel-ebot)*units.m)
uEBS = shrE[0].magnitude
vEBS = shrE[1].magnitude
EBS = np.sqrt(uEBS**2 + vEBS**2)
angEBS = 180 + np.arctan2(vEBS, uEBS)*180/np.pi

srh500 = mc.storm_relative_helicity(z*1000*units.m, u0*units('m/s'), v0*units('m/s'), depth=500*units.m, storm_u=uBR*units('m/s'), storm_v=vBR*units('m/s'))
SRH500 = srh500[2].magnitude
srh1km = mc.storm_relative_helicity(z*1000*units.m, u0*units('m/s'), v0*units('m/s'), depth=1000*units.m, storm_u=uBR*units('m/s'), storm_v=vBR*units('m/s'))
SRH1km = srh1km[2].magnitude
srh3km = mc.storm_relative_helicity(z*1000*units.m, u0*units('m/s'), v0*units('m/s'), depth=3000*units.m, storm_u=uBR*units('m/s'), storm_v=vBR*units('m/s'))
SRH3km = srh3km[2].magnitude
esrh = mc.storm_relative_helicity(z*1000*units.m, u0*units('m/s'), v0*units('m/s'), depth=(etop-ebot)*units.m, bottom=ebot*units.m, storm_u=uBR*units('m/s'), storm_v=vBR*units('m/s'))
ESRH = esrh[2].magnitude

crit_cm1 = mc.critical_angle(prs0/100*units.hPa, u0*units('m/s'), v0*units('m/s'), z*1000*units.m, uBR*units('m/s'), vBR*units('m/s'))
crit_angle = crit_cm1.magnitude
scp_cm1 = mc.supercell_composite(mucape*units('J/kg'), ESRH*units('m^2/s^2'), EBS*units('m/s'))
SCP = scp_cm1[0].magnitude
stp_cm1 = mc.significant_tornado(sbcape*units('J/kg'), Zlcl*units.m, SRH500*units('m^2/s^2'), shr06*units('m/s'))
STP = stp_cm1[0].magnitude


print("---STABLE PROFILE---")
print(f"Profile top:          {1000*max(z)+140:.0f} m")
print(f"Sfc CAPE,CIN:         {cape:.0f} J/kg, {cin:.0f} J/kg")
print(f"MUCAPE,MUCIN:         {mucape:.0f} J/kg, {mucin:.0f} J/kg")
print(f"SBCAPE,SBCIN:         {sbcape:.0f} J/kg, {sbcin:.0f} J/kg")
print(f"MLCAPE,MLCIN:         {mlcape:.0f} J/kg, {mlcin:.0f} J/kg")
print(f"LCL height:           {Zlcl:.0f} m ({Plcl:.0f} hPa)")
print(f"Eff. inflow layer:    {ebot:.0f} m-{etop:.0f} m")
print(f"Bunkers RM:           {smBR:.1f} m/s at {angBR:.0f} deg (Vector: {uBR:.1f} m/s, {vBR:.1f} m/s)")
print(f"0-6 km mean wind:     {VH06:.1f} m/s at {ang06:.0f} deg (Vector: {u06:.1f} m/s, {v06:.1f} m/s)")
print(f"-----")
print(f"Bulk shear (0-6 km):  {shr06:.1f} m/s at {angSHR:.0f} deg (Vector: {ushr06:.1f} m/s, {vshr06:.1f} m/s)")
print(f"Bulk shear (Eff.):    {EBS:.1f} m/s at {angEBS:.0f} deg (Vector: {uEBS:.1f} m/s, {vEBS:.1f} m/s)")
print(f"SRH (0-500 m):        {SRH500:.0f} m2/s2")
print(f"SRH (0-1 km):         {SRH1km:.0f} m2/s2")
print(f"SRH (0-3 km):         {SRH3km:.0f} m2/s2")
print(f"SRH (Eff.):           {ESRH:.0f} m2/s2")
print(f"Critical angle:       {crit_angle:.1f} degrees")
print(f"Supercell composite:  {SCP:.1f}")
print(f"Significant tornado:  {STP:.1f}")

#%%
figsave = False

fig = plt.figure(figsize=(8,8))

skew = SkewT(fig=fig)
skew.plot(sndsp, (sndst-273.15), '-r', linewidth=2)
skew.plot(sndsp, (sndtd-273.15), '-g', linewidth=2)
skew.plot(sndsp, np.array(T_parcel_snd.magnitude[:])-273.15, '-k', linewidth=2)
skew.plot_dry_adiabats()
skew.plot_moist_adiabats()
skew.plot_mixing_lines()
skew.ax.set_ylim(1000, 100)
skew.ax.set_xlim(-40, 30)
plt.title('Input sounding')
ax_hod = inset_axes(skew.ax, '42%', '42%', loc=1)
H = Hodograph(ax_hod, component_range=60.)
H.add_grid(increment=20)
H.plot(sndsu, sndsv, color='k', linewidth=1.5)

plt.show()
if figsave:
    plt.savefig(fp+'input_sounding.png', dpi=300)


fig = plt.figure(figsize=(8,8))

skew = SkewT(fig=fig)
skew.plot(prs0[0,:,0,0]/100., (T0-273.15), '-r', linewidth=2)
skew.plot(prs0[0,:,0,0]/100., (Td0-273.15), '-g', linewidth=2)
skew.plot(prs0[0,:,0,0]/100., np.array(T_parcel.magnitude[:])-273.15, '-k', linewidth=2)
skew.plot_dry_adiabats()
skew.plot_moist_adiabats()
skew.plot_mixing_lines()
skew.ax.set_ylim(1000, 100)
skew.ax.set_xlim(-40, 30)
plt.title('Base state')
ax_hod = inset_axes(skew.ax, '42%', '42%', loc=1)
H = Hodograph(ax_hod, component_range=60.)
H.add_grid(increment=20)
H.plot(u0[0,:,0,0]+umove, v0[0,:,0,0], color='k', linewidth=1.5)

plt.show()
if figsave:
    plt.savefig(fp+'base_state.png', dpi=300)

# fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(12,8))

# c1 = plot_cfill(ds.xf, ds.yf, th[0,:,:], 'th', ax1)
# ax1.set_title(f"mean th={thmn[0]:.1f} K, th0={th0[0]:.1f} K")
# c2 = plot_cfill(ds.xf, ds.yf, prs[0,:,:]/100, 'prs', ax2)
# ax2.set_title(f"mean prs={prsmn[0]/100:.1f} hPa, prs0={prs0[0]/100:.1f} hPa")
# c3 = plot_cfill(ds.xf, ds.yf, u[0,:,:], 'u', ax3, datalims=[-20.,20.])
# ax3.set_title(f"mean u={umn[0]:.1f} m/s, u0={u0[0]:.1f} m/s")
# c4 = plot_cfill(ds.xf, ds.yf, v[0,:,:], 'v', ax4, datalims=[-20.,20.])
# ax4.set_title(f"mean v={vmn[0]:.1f} m/s, v0={v0[0]:.1f} m/s")


# plt.show()


'''
fig, ((ax1,ax2,ax3),(ax4,ax5,ax6)) = plt.subplots(2,3,figsize=(15,8))

c1 = plot_cfill(ds.xf, ds.yf, th[0,:,:], 'th', ax1, datalims=[290.,310.])
ax1.set_title(f"mean th={thmn[0]:.1f} K, th0={th0[0]:.1f} K")
c2 = plot_cfill(ds.xf, ds.yf, prs[0,:,:]/100, 'prs', ax2, datalims=[991.,993.])
ax2.set_title(f"mean prs={prsmn[0]/100:.1f} hPa, prs0={prs0[0]/100:.1f} hPa")
c3 = plot_cfill(ds.xf, ds.yf, qv[0,:,:]*1000, 'qv', ax3, datalims=[10.,12.])
ax3.set_title(f"mean qv={qvmn[0]*1000:.1f} g/kg, qv0={qv0[0]*1000:.1f} g/kg")
c4 = plot_cfill(ds.xf, ds.yf, u[0,:,:], 'u', ax4, datalims=[-20.,20.])
ax4.set_title(f"mean u={umn[0]:.1f} m/s, u0={u0[0]:.1f} m/s")
c5 = plot_cfill(ds.xf, ds.yf, v[0,:,:], 'v', ax5, datalims=[-20.,20.])
ax5.set_title(f"mean v={vmn[0]:.1f} m/s, v0={v0[0]:.1f} m/s")
c6 = plot_cfill(ds.xf, ds.yf, ds.cref[0,:,:], 'dbz', ax6, datalims=[0.,80.])
ax6.set_title("cref")

plt.show()
'''

#%% Save to .txt file

if save_flag:
    fsave = snd_fn[62:-4]+'_stable.txt'
    hd = np.zeros((1,3))
    hd[0][0] = prsmn[0]/100
    hd[0][1] = thmn[0]
    hd[0][2] = qvmn[0]*1000
    np.savetxt(fp+fsave, hd, fmt='%f')
    
    zsave = list(z*1000)
    thsave = list(thmn)
    qvsave = list(qvmn*1000)
    usave = list(umn)
    vsave = list(vmn)
    
    dat = {'z': zsave, 'theta': thsave, 'qv': qvsave, 'u': usave, 'v': vsave}
    dfs = pd.DataFrame(data=dat, dtype=float)
    
    with open(fp+fsave, 'a') as ff:
        ff.write(dfs.to_string(header=False, index=False))


#%%

#print(f"---{file[-10:-4]}---")
print("---RAW SOUNDING---")
print(f"Profile top: {max(sndsh):.0f} m")
if max(sndsh) >= 6123:
    bwnd = mc.bunkers_storm_motion(sndsp*units.hPa, sndsu*units('m/s'), sndsv*units('m/s'), sndsh*units.m)
    uBR_snd = bwnd[0].magnitude[0]
    vBR_snd = bwnd[0].magnitude[1]
    uBL_snd = bwnd[1].magnitude[0]
    vBL_snd = bwnd[1].magnitude[1]
    u06_snd = bwnd[2].magnitude[0]
    v06_snd = bwnd[2].magnitude[1]
    print(f"0-6 shear vector: {u06_snd:.1f}, {v06_snd:.1f}")
    print(f"Bunkers RM vector: {uBR_snd:.1f}, {vBR_snd:.1f}")
else:
    print("0-6 shear vector: --")
    print("Bunkers RM vector: --")

T_parcel_snd = mc.parcel_profile(sndsp*units.hPa, sndst[0]*units.K, sndtd[0]*units.K)
SB_snd = mc.cape_cin(sndsp*units.hPa, sndst*units.K, sndtd*units.K, T_parcel_snd)
MU_snd = mc.most_unstable_cape_cin(sndsp*units.hPa, sndst*units.K, sndtd*units.K)
cape_snd = SB_snd[0].magnitude
cin_snd = SB_snd[1].magnitude
mucape_snd = MU_snd[0].magnitude
mucin_snd = MU_snd[1].magnitude
print(f"CAPE/CIN: {cape_snd:.0f}/{cin_snd:.0f}")
print(f"MUCAPE/MUCIN: {mucape_snd:.0f}/{mucin_snd:.0f}")


print("---STABLE PROFILE---")
print(f"Profile top: {max(z*1000):.0f} m")
if max(z*1000) >= 6123:
    bwnd = mc.bunkers_storm_motion(prs0/100*units.hPa, u0*units('m/s'), v0*units('m/s'), z*1000*units.m)
    uBR_cm1 = bwnd[0].magnitude[0]
    vBR_cm1 = bwnd[0].magnitude[1]
    uBL_cm1 = bwnd[1].magnitude[0]
    vBL_cm1 = bwnd[1].magnitude[1]
    u06_cm1 = bwnd[2].magnitude[0]
    v06_cm1 = bwnd[2].magnitude[1]
    print(f"0-6 shear vector: {u06_cm1:.1f}, {v06_cm1:.1f}")
    print(f"Bunkers RM vector: {uBR_cm1:.1f}, {vBR_cm1:.1f}")
else:
    print("0-6 shear vector: --")
    print("Bunkers RM vector: --")

T_parcel = mc.parcel_profile(prs0/100*units.hPa, t0[0]*units.K, (td0[0]+273.15)*units.K)
SB_cm1 = mc.cape_cin(prs0/100*units.hPa, t0*units.K, (td0+273.15)*units.K, T_parcel)
MU_cm1 = mc.most_unstable_cape_cin(prs0/100*units.hPa, t0*units.K, (td0+273.15)*units.K)
cape = SB_cm1[0].magnitude
cin = SB_cm1[1].magnitude
mucape = MU_cm1[0].magnitude
mucin = MU_cm1[1].magnitude
print(f"CAPE/CIN: {cape:.0f}/{cin:.0f}")
print(f"MUCAPE/MUCIN: {mucape:.0f}/{mucin:.0f}")





plt.rcParams['figure.figsize']=(9,9)

skew = SkewT()
skew.plot(sndsp, (sndtd-273.15), '-k', linewidth=2)
skew.plot(sndsp, (sndst-273.15), '-k', linewidth=2)
# skew.plot(prsmn/100., tdmn, '-g', linewidth=2.5)
skew.plot(prsmn/100., (tmn-273.15), '-r', linewidth=2.5)
# skew.plot(prsmn/100., T_parcel.magnitude[:]-273.15, '-b', linewidth=2.5)
# Add the relevant special lines
skew.plot_dry_adiabats()
skew.plot_moist_adiabats()
skew.plot_mixing_lines()
skew.ax.set_ylim(1000, 100)
skew.ax.set_xlim(-40, 30)
plt.title('CM1 dry stable state + input')
# Create a hodograph
ax_hod = inset_axes(skew.ax, '42%', '42%', loc=1)
h = Hodograph(ax_hod, component_range=40.)
h.add_grid(increment=20)
h.plot(sndsu, sndsv, color='k', linewidth=1.8)
h.plot(umn, vmn, color='r', linewidth=1.5)

plt.show()





plt.rcParams['figure.figsize']=(9,9)

skew2 = SkewT()
# skew2.plot(prs0/100., td0, '-k', linewidth=2)
skew2.plot(prs0/100., (t0-273.15), '-k', linewidth=2)
# skew2.plot(prsmn/100., tdmn, '-g', linewidth=2.5)
skew2.plot(prsmn/100., (tmn-273.15), '-r', linewidth=2.5)
# skew2.plot(prsmn/100., T_parcel.magnitude[:]-273.15, '-b', linewidth=2.5)
# Add special lines
skew2.plot_dry_adiabats()
skew2.plot_moist_adiabats()
skew2.plot_mixing_lines()
skew2.ax.set_ylim(1000, 100)
skew2.ax.set_xlim(-40, 30)
plt.title('CM1 dry stable state + base state')
# Create hodograph
ax_hod = inset_axes(skew2.ax, '42%', '42%', loc=1)
h2 = Hodograph(ax_hod, component_range=40.)
h2.add_grid(increment=20)
h2.plot(u0, v0, color='k', linewidth=1.8)
h2.plot(umn, vmn, color='r', linewidth=1.5)

plt.show()


#%%
save_flag = True
fp = f"/Users/morgan.schneider/Documents/cm1/PERiLS/IOP{iop}/{snd_time}/"
ff = glob(fp+'moist/cm1out_0000*.nc')
fn = ff[0]
uv = read_cm1out(fn,['zh','u0','v0','uinterp','vinterp'])

umove = 6.0

ff = glob(fp+'moist/cm1out_0000*.nc')
fn = ff[0]
pth = read_cm1out(fn,['th','qv','prs','th0','qv0','prs0'])

z = uv.z
u0 = uv.u0[0,:,0,0]
v0 = uv.v0[0,:,0,0]
th0 = pth.th0[0,:,0,0]
prs0 = pth.prs0[0,:,0,0]
t0 = th0*(prs0/100000.)**0.286
qv0 = pth.qv0[0,:,0,0]
e0 = (qv0*prs0/100) / (0.622+qv0)
td0 = 243.5 / ((17.67/(np.log(e0/6.112)))-1)

u = uv.uinterp[0,:,:,:]
v = uv.vinterp[0,:,:,:]
th = pth.th[0,:,:,:]
prs = pth.prs[0,:,:,:]
t = th*(prs/100000.)**0.286
qv = pth.qv[0,:,:,:]
e = (qv*prs/100) / (0.622+qv)
td = 243.5 / ((17.67/(np.log(e/6.112)))-1)

thmn = np.mean(th, axis=(1,2))
prsmn = np.mean(prs, axis=(1,2))
# umn = np.mean(u, axis=(1,2))
# vmn = np.mean(v, axis=(1,2))
umn = u[:,0,0]
vmn = v[:,0,0]
tmn = np.mean(t, axis=(1,2))
qvmn = np.mean(qv, axis=(1,2))
tdmn = np.mean(td, axis=(1,2))


#plt.close('all')

# plt.rcParams['figure.figsize']=(9,9)
# skew = SkewT()
# skew.plot(prsmn/100., tdmn, '-g', linewidth=2.5)
# skew.plot(prsmn/100., (tmn-273.15), '-r', linewidth=2.5)
# # Add the relevant special lines
# skew.plot_dry_adiabats()
# skew.plot_moist_adiabats()
# skew.plot_mixing_lines()
# skew.ax.set_ylim(1000, 100)
# skew.ax.set_xlim(-40, 30)
# plt.title('CM1 stable state')
# # Create a hodograph
# ax_hod = inset_axes(skew.ax, '42%', '42%', loc=1)
# h = Hodograph(ax_hod, component_range=40.)
# h.add_grid(increment=20)
# h.plot(umn, vmn, color='k', linewidth=1.8)

# plt.show()


plt.rcParams['figure.figsize']=(9,9)
skew2 = SkewT()
skew2.plot(prs0/100., td0, '-g', linewidth=2.5)
skew2.plot(prs0/100., (t0-273.15), '-r', linewidth=2.5)
# Add special lines
skew2.plot_dry_adiabats()
skew2.plot_moist_adiabats()
skew2.plot_mixing_lines()
skew2.ax.set_ylim(1000, 100)
skew2.ax.set_xlim(-40, 30)
plt.title('CM1 base state + stable winds')
# Create hodograph
ax_hod = inset_axes(skew2.ax, '42%', '42%', loc=1)
h2 = Hodograph(ax_hod, component_range=40.)
h2.add_grid(increment=20)
h2.plot(umn, vmn, color='k', linewidth=1.8)
plt.show()


plt.rcParams['figure.figsize']=(9,9)
skew3 = SkewT()
skew3.plot(sndsp, (sndtd-273.15), '-g', linewidth=2.5)
skew3.plot(sndsp, (sndst-273.15), '-r', linewidth=2.5)
# Add the relevant special lines
skew3.plot_dry_adiabats()
skew3.plot_moist_adiabats()
skew3.plot_mixing_lines()
skew3.ax.set_ylim(1000, 100)
skew3.ax.set_xlim(-40, 30)
plt.title('CM1 input state')
# Create a hodograph
ax_hod = inset_axes(skew3.ax, '42%', '42%', loc=1)
h3 = Hodograph(ax_hod, component_range=40.)
h3.add_grid(increment=20)
h3.plot(sndsu, sndsv, color='b', linewidth=1.8)
plt.show()


fig,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(12,10))

ax1.plot(umn, z*1000, '-k', linewidth=2)
ax1.plot(sndsu, sndsh, '-r')
ax1.set_xlabel('u (m/s)')
ax1.set_ylabel('Height (m)')
ax1.set_xlim(-20,50)
ax1.set_ylim(0,16000)
ax1.grid(b=True)

ax2.plot(vmn, z*1000, '-k', linewidth=2)
ax2.plot(sndsv, sndsh, '-r')
ax2.set_xlabel('v (m/s)')
ax2.set_xlim(-5,40)
ax2.set_ylim(0,16000)
ax2.grid(b=True)

ax3.plot((umn**2 + vmn**2)**0.5, zh*1000, '-k', linewidth=2)
ax3.plot(sndws, sndsh, '-r')
ax3.set_xlabel('wspd (m/s)')
ax3.set_xlim(0,60)
ax3.set_ylim(0,16000)
ax3.grid(b=True)
plt.show()



print("---STABLE PROFILE---")
print(f"Profile top: {max(zh*1000):.0f} m")
if max(z*1000) >= 6123:
    bwnd = mc.bunkers_storm_motion(prs0/100*units.hPa, umn*units('m/s'), vmn*units('m/s'), z*1000*units.m)
    uBR_cm1 = bwnd[0].magnitude[0]
    vBR_cm1 = bwnd[0].magnitude[1]
    uBL_cm1 = bwnd[1].magnitude[0]
    vBL_cm1 = bwnd[1].magnitude[1]
    u06_cm1 = bwnd[2].magnitude[0]
    v06_cm1 = bwnd[2].magnitude[1]
    print(f"0-6 shear vector: {u06_cm1:.1f}, {v06_cm1:.1f}")
    print(f"Bunkers RM vector: {uBR_cm1:.1f}, {vBR_cm1:.1f}")
else:
    print("0-6 shear vector: --")
    print("Bunkers RM vector: --")

T_parcel = mc.parcel_profile(prs0/100*units.hPa, t0[0]*units.K, (td0[0]+273.15)*units.K)
SB_cm1 = mc.cape_cin(prs0/100*units.hPa, t0*units.K, (td0+273.15)*units.K, T_parcel)
MU_cm1 = mc.most_unstable_cape_cin(prs0/100*units.hPa, t0*units.K, (td0+273.15)*units.K)
cape = SB_cm1[0].magnitude
cin = SB_cm1[1].magnitude
mucape = MU_cm1[0].magnitude
mucin = MU_cm1[1].magnitude
print(f"CAPE/CIN: {cape:.0f}/{cin:.0f}")
print(f"MUCAPE/MUCIN: {mucape:.0f}/{mucin:.0f}")



if save_flag:
    fsave = snd_fn[62:-4]+'_stable.txt'
    hd = np.zeros((1,3))
    hd[0][0] = prs0[0]/100
    hd[0][1] = th0[0]
    hd[0][2] = qv0[0]*1000
    np.savetxt(fp+fsave, hd, fmt='%f')
    
    zsave = list(z*1000)
    thsave = list(th0)
    qvsave = list(qv0*1000)
    usave = list(umn)
    vsave = list(vmn)
    
    dat = {'z': zsave, 'theta': thsave, 'qv': qvsave, 'u': usave, 'v': vsave}
    dfs = pd.DataFrame(data=dat, dtype=float)
    
    with open(fp+fsave, 'a') as ff:
        ff.write(dfs.to_string(header=False, index=False))
        







