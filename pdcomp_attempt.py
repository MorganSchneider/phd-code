#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 14:37:26 2023

@author: morgan.schneider

Unintegrated pressure perturbation decomposition (just the terms in the diagnostic equation)
"""

from CM1utils import *
import os


def poiss_solver(F,x,y,z,cfa,cfb,cfc,ad1):
    lgbth = np.zeros(shape=(len(x),len(y),len(z),), dtype=complex)
    lgbph = np.zeros(shape=(len(x),len(y),len(z),), dtype=complex)
    pdt = np.zeros(shape=(len(x),len(y),len(z),), dtype=complex)
    deft2 = np.zeros(shape=(len(z),))
    def1 = np.zeros(shape=(len(x),len(y),len(z),), dtype=float)
    rhs = np.zeros(shape=(len(x),len(y),), dtype=complex)
    
    print('FFT...')
    for k in range(len(z)):
        for jj in range(len(y)):
            for ii in range(len(x)):
                rhs[ii,jj] = F[0,k,jj,ii]*ad1[k] + 0j
        trans = np.fft.fft2(rhs)
        
        for jj in range(len(y)):
            for ii in range(len(x)):
                if k+1 == 1:
                    tem = 1/cfb[ii,jj,0]
                    lgbth[ii,jj,0] = -cfc[0] * tem
                    lgbph[ii,jj,0] = trans[ii,jj] * tem
                elif k+1 < len(z):
                    temc = 1/(cfa[k]*lgbth[ii,jj,k-1] + cfb[ii,jj,k])
                    lgbth[ii,jj,k] = -cfc[k] * temc
                    lgbph[ii,jj,k] = (trans[ii,jj] - cfa[k]*lgbph[ii,jj,k-1]) * temc
                elif k+1 == len(z):
                    temc = 1/(cfa[k]*lgbth[ii,jj,k-1] + cfb[ii,jj,k])
                    lgbth[ii,jj,k] = -cfc[k] * temc
                    lgbph[ii,jj,k] = (trans[ii,jj] - cfa[k]*lgbph[ii,jj,k-1]) * temc
                    pdt[ii,jj,k] = lgbph[ii,jj,k]
        deft2[k] = np.real(trans[0,0])

    r1 = np.zeros(shape=(len(z),), dtype=float)
    print('Backward tridiagonal solver...')
    for k in np.arange(len(z)-2,-1,-1):
        for jj in range(len(y)):
            for ii in range(len(x)):
                pdt[ii,jj,k] = lgbth[ii,jj,k]*pdt[ii,jj,k+1] + lgbph[ii,jj,k]
    r1[0] = 0.
    pdt[0,0,0] = 0.
    r1[1] = (deft2[0] - cfb[0,0,0]*r1[0])/cfc[0]
    pdt[0,0,1] = r1[1] + 0j

    for k in np.arange(1,len(z)-1):
        r1[k+1] = (deft2[k] - cfa[k]*r1[k-1] - cfb[0,0,k]*r1[k])/cfc[k]
        pdt[0,0,k+1] = r1[k+1] + 0j

    print('IFFT...')
    ppi = np.zeros(shape=(len(z),len(y),len(x),), dtype=float)
    for k in range(len(z)):
        rhs = pdt[:,:,k]
        trans = np.fft.ifft2(rhs)
        for jj in range(len(y)):
            for ii in range(len(x)):
                ppi[k,jj,ii] = np.real(trans[ii,jj])
    print('Done!')
    
    return ppi
    

def calc_cterms(x,y,z,zf,rho0,thr0,dx,dz):
    from scipy.interpolate import interp1d
    
    zzh,mh,mf = z_extend(z,zf,dz)
    
    fr = interp1d(z, rho0[0,:,1,1], fill_value='extrapolate')
    rf0 = fr(zf)
    ft = interp1d(z, thr0[0,:,1,1], fill_value='extrapolate')
    thrf0 = ft(zzh)

    cfa = np.zeros(shape=(len(z),), dtype=float)
    cfc = np.zeros(shape=(len(z),), dtype=float)
    cfb = np.zeros(shape=(len(x),len(y),len(z),), dtype=float)
    ad1 = np.zeros(shape=(len(z),), dtype=float)
    for k in np.arange(1,len(zf)):
        cfa[k-1] = mh[k]*mf[k]*rf0[k-1] * 0.5 * (thrf0[k-1]+thrf0[k]) / (rho0[0,k-1,1,1]*thrf0[k]*dz**2)
        cfc[k-1] = mh[k]*mf[k+1]*rf0[k] * 0.5 * (thrf0[k]+thrf0[k+1]) / (rho0[0,k-1,1,1]*thrf0[k]*dz**2)
        ad1[k-1] = 1/(1004.5*rho0[0,k-1,1,1]*thrf0[k])
    cfa[0] = 0
    cfc[-1] = 0

    for k in np.arange(0,len(z)):
        for j in np.arange(0,len(y)):
            for i in np.arange(0,len(x)):
                cfb[i,j,k] = 2 * (np.cos(2*np.pi*i/len(x)) + np.cos(2*np.pi*j/len(y)) - 2)/(dx**2) - cfa[k] - cfc[k]
    
    return cfa,cfb,cfc,ad1


def z_extend(z,zf,dz):
    zzh = np.zeros(shape=(len(z)+2,), dtype=float)
    zzf = np.zeros(shape=(len(zf)+2,), dtype=float)
    mf = np.zeros(shape=(len(zf)+2,), dtype=float)

    zzf[1:-1] = zf
    zzf[0] = -zf[1]
    zzf[-1] = zf[-1] + (zf[-1] - zf[-2])

    zzh[1:-1] = 0.5 * (zf[1:]+zf[0:-1])
    zzh[0] = -z[0]
    zzh[-1] = z[-1] + 2*(zf[-1]-z[-1])

    mh = dz/(np.diff(zzf))
    mf[1:-1] = dz/(np.diff(zzh))
    mf[0] = mf[1]
    mf[-1] = mf[-2]
    
    return zzh,mh,mf





#%%

fp = '/Volumes/Promise_Pegasus_70TB/2021-SupSims-final/'
# sim_dirs = os.listdir(fp)[1:-1]
fd = 'R27_rerun'
fn = fp+fd+'/cm1out_000064.nc'

figsave = False


ds = nc.Dataset(fp+fd+'/cm1out_000001.nc')
xh = ds.variables['xh'][:].data
xf = ds.variables['xf'][:].data
yh = ds.variables['yh'][:].data
yf = ds.variables['yf'][:].data
z = ds.variables['z'][:].data
zf = ds.variables['zf'][:].data

prs0 = ds.variables['prs'][:].data
pi0 = ds.variables['pi'][:].data
th0 = ds.variables['th'][:].data
qv0 = ds.variables['qv'][:].data
u0 = ds.variables['uinterp'][:].data
v0 = ds.variables['vinterp'][:].data
ds.close()


ds = nc.Dataset(fn)
time = ds.variables['time'][:].data

prs = ds.variables['prs'][:].data
pi = ds.variables['pi'][:].data
th = ds.variables['th'][:].data
qv = ds.variables['qv'][:].data
qx = (ds.variables['qc'][:].data + ds.variables['qr'][:].data + 
      ds.variables['qi'][:].data + ds.variables['qs'][:].data + 
      ds.variables['qg'][:].data + ds.variables['qhl'][:].data)
u = ds.variables['uinterp'][:].data
v = ds.variables['vinterp'][:].data
w = ds.variables['winterp'][:].data
zvort = ds.variables['zvort'][:].data
ds.close()

prspert = prs - prs0
pipert = pi - pi0
thpert = th - th0
qvpert = qv - qv0
thr0 = th0 * (1 + 0.61*qv0)
thr = th * (1 + 0.61*qv - qx)
thrpert = thr - thr0
B = 9.8 * (thpert/th0 + 0.61*qvpert - qx)
upert = u - u0
vpert = v - v0

du_dx = np.gradient(upert, xh*1000, axis=3)
du_dy = np.gradient(upert, yh*1000, axis=2)
du_dz = np.gradient(upert, z*1000, axis=1)
dum_dz = np.gradient(u0, z*1000, axis=1)

dv_dx = np.gradient(vpert, xh*1000, axis=3)
dv_dy = np.gradient(vpert, yh*1000, axis=2)
dv_dz = np.gradient(vpert, z*1000, axis=1)
dvm_dz = np.gradient(v0, z*1000, axis=1)

dw_dx = np.gradient(w, xh*1000, axis=3)
dw_dy = np.gradient(w, yh*1000, axis=2)
dw_dz = np.gradient(w, z*1000, axis=1)

dB_dz = np.gradient(B, z*1000, axis=1)

divh = du_dx + dv_dy

F_fe = du_dx**2 + dv_dy**2 + dw_dz**2
F_ns = 2*((dv_dx*du_dy) + (dw_dx*du_dz) + (dw_dy*dv_dz))
F_nd = F_fe + F_ns
F_ld = 2*((dw_dx*dum_dz) + (dw_dy*dvm_dz))
F_b = -dB_dz

F = struct()
F.SetAttr('nd', F_nd[0,:,:,:])
F.SetAttr('ld', F_ld[0,:,:,:])
F.SetAttr('b', F_b[0,:,:,:])



fig,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(12,8))

iz = np.where(z>=0)[0][0]
lims = [-0.02,0.02]
cm = 'bwr'
qix = 6
xlims = [-10,10]
ylims = [-10,10]

c1 = plot_cfill(xf, yf, F_nd[0,iz,:,:]+F_ld[0,iz,:,:]+F_b[0,iz,:,:], 'pipert', ax1, datalims=lims, cmap=cm)
# ax1.quiver(xh[::qix], yh[::qix], u[0,iz,::qix,::qix], v[0,iz,::qix,::qix], scale=450, width=0.003, pivot='middle')
ax1.set_title('Forcing - total')
ax1.set_xlim(xlims)
ax1.set_ylim(ylims)

c2 = plot_cfill(xf, yf, F_nd[0,iz,:,:], 'pipert', ax2, datalims=lims, cmap=cm)
# ax2.quiver(xh[::qix], yh[::qix], u[0,iz,::qix,::qix], v[0,iz,::qix,::qix], scale=450, width=0.003, pivot='middle')
ax2.set_title('Forcing - nonlinear dynamic')
ax2.set_xlim(xlims)
ax2.set_ylim(ylims)

c3 = plot_cfill(xf, yf, F_ld[0,iz,:,:], 'pipert', ax3, datalims=lims, cmap=cm)
# ax3.quiver(xh[::qix], yh[::qix], u[0,iz,::qix,::qix], v[0,iz,::qix,::qix], scale=450, width=0.003, pivot='middle')
ax3.set_title('Forcing - linear dynamic')
ax3.set_xlim(xlims)
ax3.set_ylim(ylims)

c4 = plot_cfill(xf, yf, F_b[0,iz,:,:], 'pipert', ax4, datalims=lims, cmap=cm)
# ax4.quiver(xh[::qix], yh[::qix], u[0,iz,::qix,::qix], v[0,iz,::qix,::qix], scale=450, width=0.003, pivot='middle')
ax4.set_title('Forcing - buoyancy')
ax4.set_xlim(xlims)
ax4.set_ylim(ylims)

plt.suptitle(f"z = {z[iz]*1000:.0f} m, t = {time[0]/60:.0f} min")
if figsave:
    plt.savefig('/Volumes/Promise_Pegasus_70TB/general/'+fd+f"/TLV_plots/ppforcing_{z[iz]*1000:.0f}m_{time[0]/60:.0f}min.png", dpi=400)


fig,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(12,8))

c1 = plot_cfill(xf, yf, thrpert[0,iz,:,:], 'thrpert', ax1, datalims=[-10,10])
ax1.quiver(xh[::qix], yh[::qix], u[0,iz,::qix,::qix], v[0,iz,::qix,::qix], scale=450, width=0.003, pivot='middle')
ax1.set_title('thrpert')
ax1.set_xlim(xlims)
ax1.set_ylim(ylims)

c2 = plot_cfill(xf, yf, zvort[0,iz,:,:], 'zvort', ax2, datalims=[-0.3,0.3])
ax2.quiver(xh[::qix], yh[::qix], u[0,iz,::qix,::qix], v[0,iz,::qix,::qix], scale=450, width=0.003, pivot='middle')
ax2.set_title('zvort')
ax2.set_xlim(xlims)
ax2.set_ylim(ylims)

c3 = plot_cfill(xf, yf, prspert[0,iz,:,:]/100, 'prspert', ax3, datalims=[-10,10])
ax3.quiver(xh[::qix], yh[::qix], u[0,iz,::qix,::qix], v[0,iz,::qix,::qix], scale=450, width=0.003, pivot='middle')
ax3.set_title('prspert')
ax3.set_xlim(xlims)
ax3.set_ylim(ylims)

c4 = plot_cfill(xf, yf, divh[0,iz,:,:], 'divh', ax4, datalims=[-0.3,0.3])
ax4.quiver(xh[::qix], yh[::qix], u[0,iz,::qix,::qix], v[0,iz,::qix,::qix], scale=450, width=0.003, pivot='middle')
ax4.set_title('divergence')
ax4.set_xlim(xlims)
ax4.set_ylim(ylims)

plt.suptitle(f"z = {z[iz]*1000:.0f} m, t = {time[0]/60:.0f} min")
if figsave:
    plt.savefig('/Volumes/Promise_Pegasus_70TB/general/'+fd+f"/TLV_plots/thr+vort+div_{z[iz]*1000:.0f}m_{time[0]/60:.0f}min.png", dpi=400)




#%%

# Only need to do this if the thermodynamic base state changes - takes a long time to run
recalc_c = True

fp = '/Volumes/Promise_Pegasus_70TB/2021-SupSims-final/'
# sim_dirs = os.listdir(fp)[1:-1]
fd = 'R10_rerun'
fn = fp+fd+'/cm1out_000060.nc'


ds = nc.Dataset(fp+fd+'/cm1out_000001.nc')
xh = ds.variables['xh'][:].data
xf = ds.variables['xf'][:].data
yh = ds.variables['yh'][:].data
yf = ds.variables['yf'][:].data
z = ds.variables['z'][:].data
zf = ds.variables['zf'][:].data

prs0 = ds.variables['prs'][:].data
pi0 = ds.variables['pi'][:].data
th0 = ds.variables['th'][:].data
qv0 = ds.variables['qv'][:].data
u0 = ds.variables['u'][:].data
v0 = ds.variables['v'][:].data
ui0 = ds.variables['uinterp'][:].data
vi0 = ds.variables['vinterp'][:].data
ds.close()


ds = nc.Dataset(fn)
time = ds.variables['time'][:].data

prs = ds.variables['prs'][:].data
pi = ds.variables['pi'][:].data
th = ds.variables['th'][:].data
qv = ds.variables['qv'][:].data
qx = (ds.variables['qc'][:].data + ds.variables['qr'][:].data + 
      ds.variables['qi'][:].data + ds.variables['qs'][:].data + 
      ds.variables['qg'][:].data + ds.variables['qhl'][:].data)
u = ds.variables['u'][:].data
v = ds.variables['v'][:].data
w = ds.variables['w'][:].data
ui = ds.variables['uinterp'][:].data
vi = ds.variables['vinterp'][:].data
wi = ds.variables['winterp'][:].data
ds.close()

# Matt's pressure perturbations
ds = nc.Dataset(fp+fd+'/pp/dyn_000060.nc')
p_dl = ds.variables['p_dl'][:].data
pi_dl = ds.variables['pi_dl'][:].data
p_b = ds.variables['p_b'][:].data
pi_b = ds.variables['pi_b'][:].data
p_dn = ds.variables['p_dn'][:].data
pi_dn = ds.variables['pi_dn'][:].data
ds.close()


# Random calculations
prspert = prs - prs0
pipert = pi - pi0
thpert = th - th0
qvpert = qv - qv0
thr0 = th0 * (1 + 0.61*qv0)
thr = th * (1 + 0.61*qv - qx)
thrpert = thr - thr0
#%
dx = 0.5 # set in namelist.input
dz = 0.2 # set in namelist.input
t0 = th0*pi0
rho0 = prs0/(287.5*t0)

# Numerical junk
if recalc_c:
    cfa,cfb,cfc,ad1 = calc_cterms(1000*xh,1000*yh,1000*z,1000*zf,rho0,thr0,1000*dx,1000*dz)

#%
# Calculate forcing terms
B = 9.8 * (thpert/th0 + 0.61*qvpert - qx)
upert = ui - ui0
vpert = vi - vi0

du_dx = np.gradient(upert, 1000*xh, axis=3)
du_dy = np.gradient(upert, 1000*yh, axis=2)
du_dz = np.gradient(upert, 1000*z, axis=1)
dum_dz = np.gradient(ui0, 1000*z, axis=1)

dv_dx = np.gradient(vpert, 1000*xh, axis=3)
dv_dy = np.gradient(vpert, 1000*yh, axis=2)
dv_dz = np.gradient(vpert, 1000*z, axis=1)
dvm_dz = np.gradient(vi0, 1000*z, axis=1)

dw_dx = np.gradient(wi, 1000*xh, axis=3)
dw_dy = np.gradient(wi, 1000*yh, axis=2)
dw_dz = np.gradient(wi, 1000*z, axis=1)

dB_dz = np.gradient(B, 1000*z, axis=1)

divh = du_dx + dv_dy

# F_fe = -rho0*(du_dx**2 + dv_dy**2 + dw_dz**2)
# F_ns = -2*rho0*((dv_dx*du_dy) + (dw_dx*du_dz) + (dw_dy*dv_dz))
F_nd = -rho0*(du_dx**2 + dv_dy**2 + dw_dz**2) - 2*rho0*((dv_dx*du_dy) + (dw_dx*du_dz) + (dw_dy*dv_dz))
F_ld = -2*rho0*((dw_dx*dum_dz) + (dw_dy*dvm_dz))
F_b = rho0*dB_dz





print('Solving linear dynamic!')
ppi_ld = poiss_solver(F_ld,xh,yh,z,cfa,cfb,cfc,ad1)

print('Solving buoyancy!')
ppi_b = poiss_solver(F_b,xh,yh,z,cfa,cfb,cfc,ad1)

print('Solving nonlinear dynamic!')
ppi_nd = poiss_solver(F_nd,xh,yh,z,cfa,cfb,cfc,ad1)



#%

pp_ld = 1004.5*thr0[0,:,:,:]*rho0[0,:,:,:]*ppi_ld
pp_nd = 1004.5*thr0[0,:,:,:]*rho0[0,:,:,:]*ppi_nd
pp_b = 1004.5*thr0[0,:,:,:]*rho0[0,:,:,:]*ppi_b


lims = [-10,10]


fig,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(12,9))

plot_cfill(xf, yf, pp_nd[0,:,:]/100, 'prspert', ax1, datalims=lims, xlims=[-10,10], ylims=[-10,10])
ax1.set_title("nonlinear dyn p'")

plot_cfill(xf, yf, pp_ld[0,:,:]/100, 'prspert', ax2, datalims=lims, xlims=[-10,10], ylims=[-10,10])
ax2.set_title("linear dyn p'")

plot_cfill(xf, yf, pp_b[0,:,:]/100, 'prspert', ax3, datalims=lims, xlims=[-10,10], ylims=[-10,10])
ax3.set_title("buoyancy p'")

# plot_cfill(xf, yf, (pp_ld[0,:,:]+pp_nd[0,:,:]+pp_b[0,:,:])/100, 'prspert', ax4, datalims=lims, xlims=[-10,10], ylims=[-10,10])
plot_cfill(xf, yf, (prspert[0,0,:,:]-pp_ld[0,:,:]-pp_b[0,:,:])/100, 'prspert', ax4, datalims=lims, xlims=[-10,10], ylims=[-10,10])
ax4.set_title("nonlinear residual p'")



# fig,ax = plt.subplots(1,1,figsize=(8,6))
# plot_cfill(xf, yf, prspert[0,0,:,:]/100, 'prspert', ax, datalims=lims, xlims=[-10,10], ylims=[-10,10])
# ax.set_title("cm1 p'")




fig,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(12,9))

plot_cfill(xf, yf, p_dn[0,:,:]/100, 'prspert', ax1, datalims=lims, xlims=[-10,10], ylims=[-10,10])
ax1.set_title("nonlinear dyn p'")

plot_cfill(xf, yf, p_dl[0,:,:]/100, 'prspert', ax2, datalims=lims, xlims=[-10,10], ylims=[-10,10])
ax2.set_title("linear dyn p'")

plot_cfill(xf, yf, p_b[0,:,:]/100, 'prspert', ax3, datalims=lims, xlims=[-10,10], ylims=[-10,10])
ax3.set_title("buoyancy p'")

plot_cfill(xf, yf, (p_dn[0,:,:]+p_dl[0,:,:]+p_b[0,:,:])/100, 'prspert', ax4, datalims=lims, xlims=[-10,10], ylims=[-10,10])
ax4.set_title("total p'")







#%%

lims = [-0.01,0.01]

fig,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(12,9))

plot_cfill(xf, yf, ppi_nd[0,:,:], 'prspert', ax1, datalims=lims, xlims=[-10,10], ylims=[-10,10])
ax1.set_title("nonlinear dyn pi'")

plot_cfill(xf, yf, ppi_ld[0,:,:], 'prspert', ax2, datalims=lims, xlims=[-10,10], ylims=[-10,10])
ax2.set_title("linear dyn pi'")

plot_cfill(xf, yf, ppi_b[0,:,:], 'prspert', ax3, datalims=lims, xlims=[-10,10], ylims=[-10,10])
ax3.set_title("buoyancy pi'")

# plot_cfill(xf, yf, (pp_ld[0,:,:]+pp_nd[0,:,:]+pp_b[0,:,:])/100, 'prspert', ax4, datalims=lims, xlims=[-10,10], ylims=[-10,10])
plot_cfill(xf, yf, (pipert[0,0,:,:]-ppi_ld[0,:,:]-ppi_b[0,:,:]), 'prspert', ax4, datalims=lims, xlims=[-10,10], ylims=[-10,10])
ax4.set_title("nonlinear residual pi'")




fig,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(12,9))

plot_cfill(xf, yf, pi_dn[0,:,:], 'prspert', ax1, datalims=lims, xlims=[-10,10], ylims=[-10,10])
ax1.set_title("nonlinear dyn pi'")

plot_cfill(xf, yf, pi_dl[0,:,:], 'prspert', ax2, datalims=lims, xlims=[-10,10], ylims=[-10,10])
ax2.set_title("linear dyn pi'")

plot_cfill(xf, yf, pi_b[0,:,:], 'prspert', ax3, datalims=lims, xlims=[-10,10], ylims=[-10,10])
ax3.set_title("buoyancy pi'")

plot_cfill(xf, yf, (pi_dn[0,:,:]+pi_dl[0,:,:]+pi_b[0,:,:]), 'prspert', ax4, datalims=lims, xlims=[-10,10], ylims=[-10,10])
ax4.set_title("total pi'")


















