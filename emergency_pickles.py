#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 14:15:52 2024

@author: morgan.schneider
"""

from CM1utils import *


def make_pickle(pkl_name, data, key, new_pkl=False):
    if new_pkl:
        tmp = {f"{key}":data}
        dbfile = open(pkl_name, 'wb')
        pickle.dump(tmp, dbfile)
        dbfile.close()
    else:
        dbfile = open(pkl_name, 'rb')
        tmp = pickle.load(dbfile)
        dbfile.close()
        
        tmp.update({f"{key}":data})
        dbfile = open(pkl_name, 'wb')
        pickle.dump(tmp, dbfile)
        dbfile.close()

#%%


fp = '/Volumes/Promise_Pegasus_70TB/merger/merger-125m/'
ip = '/Users/morgan.schneider/Documents/merger/emergency-pickles/'

fnums1 = [28, 33, 38, 39, 40, 41, 42, 43, 48, 53, 54, 55, 56, 57, 58]
fnums2 = np.linspace(14, 73, 60)


# fn0 = 'coords.pkl'
# fn1 = 'base.pkl'
# fn2 = 'u.pkl'
# fn3 = 'v.pkl'
# fn4 = 'w.pkl'
# fn5 = 'xvort.pkl'
# fn6 = 'yvort.pkl'
# fn7 = 'zvort.pkl'
# fn8 = 'prs.pkl'
# fn9 = 'rho.pkl'
# fn10 = 'dbz.pkl'
# fn11 = 'thrpert.pkl'
# fn12 = 'buoy.pkl'


ds = nc.Dataset(fp+'base/cm1out_000013.nc')
time = ds.variables['time'][:].data[0]
xh = ds.variables['xh'][:].data
xf = ds.variables['xf'][:].data
yh = ds.variables['yh'][:].data
yf = ds.variables['yf'][:].data
z = ds.variables['z'][:].data
zf = ds.variables['zf'][:].data

xlims = [-50, 10]
ylims = [-90, -30]
zlims = [0, 5]
ix = slice(np.where(xh >= xlims[0])[0][0], np.where(xh >= xlims[1])[0][1])
iy = slice(np.where(yh >= ylims[0])[0][0], np.where(yh >= ylims[1])[0][1])
iz = slice(0, np.where(zh >= zlims[1])[0][1])

u0 = ds.variables['u0'][:].data[0,iz,0,0]
v0 = ds.variables['v0'][:].data[0,iz,0,0]
prs0 = ds.variables['prs0'][:].data[0,iz,0,0]
ds.close()


ds = nc.Dataset(fp+'pp/dyn_000014.nc')
thr0 = ds.variables['thr0'][:].data[iz]
ds.close()


ds = nc.Dataset(fp+'cm1out_pdata.nc')
ptime = ds.variables['time'][:].data
ds.close()


print('...coords.pkl...')
dbf0 = open(ip+"coords.pkl", 'wb')
tmp = {'time':time, 'xh':xh, 'xf':xf, 'yh':yh, 'yf':yf, 'z':z, 'zf':zf,
       'xlims':xlims, 'ylims':ylims, 'zlims':zlims}
pickle.dump(tmp, dbf0)
dbf0.close()

print('...base.pkl...')
dbf1 = open(ip+"base.pkl", 'wb')
tmp = {'u0':u0, 'v0':v0, 'prs0':prs0, 'thr0':thr0}
pickle.dump(tmp, dbf1)
dbf1.close()

print('...ptime.pkl...')
dbft = open(ip+"ptime.pkl", 'wb')
tmp = {'ptime':ptime}
pickle.dump(tmp, dbft)
dbft.close()

#%%
dbf2 = open(ip+"u.pkl", 'wb'); pickle.dump({},dbf2); dbf2.close()
dbf3 = open(ip+"v.pkl", 'wb'); pickle.dump({},dbf3); dbf3.close()
dbf4 = open(ip+"w.pkl", 'wb'); pickle.dump({},dbf4); dbf4.close()
dbf5 = open(ip+"xvort.pkl", 'wb'); pickle.dump({},dbf5); dbf5.close()
dbf6 = open(ip+"yvort.pkl", 'wb'); pickle.dump({},dbf6); dbf6.close()
dbf7 = open(ip+"zvort.pkl", 'wb'); pickle.dump({},dbf7); dbf7.close()
dbf8 = open(ip+"prs.pkl", 'wb'); pickle.dump({},dbf8); dbf8.close()
dbf9 = open(ip+"rho.pkl", 'wb'); pickle.dump({},dbf9); dbf9.close()
dbf10 = open(ip+"dbz.pkl", 'wb'); pickle.dump({},dbf10); dbf10.close()
dbf11 = open(ip+"thrpert.pkl", 'wb'); pickle.dump({},dbf11); dbf11.close()
dbf12 = open(ip+"buoy.pkl", 'wb'); pickle.dump({},dbf12); dbf12.close()



for fnum in fnums1:
    print(f"...pickles for cm1out_{fnum:06d}...")
    
    ds = nc.Dataset(fp+f"cm1out_{fnum:06d}.nc")
    time = ds.variables['time'][:].data[0]
    
    u = ds.variables['uinterp'][:].data[0,iz,iy,ix]
    make_pickle(ip+"u.pkl", u, time/60)
    del u
    
    v = ds.variables['vinterp'][:].data[0,iz,iy,ix]
    make_pickle(ip+"v.pkl", v, time/60)
    del v
    
    w = ds.variables['winterp'][:].data[0,iz,iy,ix]
    make_pickle(ip+"w.pkl", w, time/60)
    del w
    
    xvort = ds.variables['xvort'][:].data[0,iz,iy,ix]
    make_pickle(ip+"xvort.pkl", xvort, time/60)
    del xvort
    
    yvort = ds.variables['yvort'][:].data[0,iz,iy,ix]
    make_pickle(ip+"yvort.pkl", yvort, time/60)
    del yvort
    
    zvort = ds.variables['zvort'][:].data[0,iz,iy,ix]
    make_pickle(ip+"zvort.pkl", zvort, time/60)
    del zvort
    
    prs = ds.variables['prs'][:].data[0,iz,iy,ix]
    make_pickle(ip+"prs.pkl", prs, time/60)
    del prs
    
    rho = ds.variables['rho'][:].data[0,iz,iy,ix]
    make_pickle(ip+"rho.pkl", rho, time/60)
    del rho
    
    if fnum == 43:
        ixdbz = slice(np.where(xh >= -25)[0][0], np.where(xh >= 0)[0][1])
        iydbz = slice(np.where(yh >= -100)[0][0], np.where(yh >= -50)[0][1])
        dbz = ds.variables['dbz'][:].data[0,0,iydbz,ixdbz]
        make_pickle(ip+"dbz.pkl", dbz, time/60)
        del dbz
    elif fnum == 58:
        ixdbz = slice(np.where(xh >= -17)[0][0], np.where(xh >= 8)[0][1])
        iydbz = slice(np.where(yh >= -85)[0][0], np.where(yh >= -35)[0][1])
        dbz = ds.variables['dbz'][:].data[0,0,iydbz,ixdbz]
        make_pickle(ip+"dbz.pkl", dbz, time/60)
        del dbz
    ds.close()
    
    ds = nc.Dataset(fp+f"pp/dyn_{fnum:06d}.nc")
    thrpert = ds.variables['thrpert'][:].data[iz,iy,ix]
    make_pickle(ip+"thrpert.pkl", thrpert, time/60)
    
    buoy = 9.8 * (thrpert / np.moveaxis(np.tile(thr0, (len(xh[ix]), len(xh[ix]), 1)), -1, 0))
    make_pickle(ip+"buoy.pkl", buoy, time/60)
    del thrpert,buoy
    ds.close()
    
    
    
    
    
    










