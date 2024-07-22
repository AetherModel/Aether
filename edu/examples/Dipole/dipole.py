#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

def get_lshell(x, y, z):
    # x, y, z need to be normalized to radius of planet!

    xy = np.sqrt(x * x + y * y)
    xyz = np.sqrt(x * x + y * y + z * z)
    cosLat = xy / xyz
    lshell = 1.0 / (cosLat * cosLat)
    return lshell

def get_lshell_sphere(lat, r):
    # r needs to be normalized to radius of planet!
    cosLat = np.cos(lat)
    lshell = r / (cosLat * cosLat)
    return lshell

def get_lat_from_r_and_lshell(r, lshell):
    cosLat = np.sqrt(r / lshell)
    cosLat[cosLat > 1] = 1.0
    cosLat[cosLat < -1] = -1.0
    lat = np.arccos(cosLat)
    return lat

def get_bfield(lat, r):
    # r should be normalized to radius of planet
    # b is normalized to Bequator at surface
    oor3 = (1.0 / r)**3
    br = -2 * oor3 * np.cos(lat)
    bt = -1 * oor3 * np.sin(lat)
    return bt, br

def get_uniform_spacing(lat, rmin, rmax, nPts):
    lshell = get_lshell_sphere(lat, 1.0)
    print('building field line for lat : ', lat)
    print(' --> lshell : ', lshell)
    if (lshell < rmax):
        rmax = lshell
        print(' --> lshell is limiter!')
    else:
        print(' --> rmax is limiter : ', rmax)
    r = np.linspace(rmin, rmax, num = nPts)
    return r

def get_r3_spacing(lat, rmin, rmax, nPts):
    lshell = get_lshell_sphere(lat, 1.0)
    print('building field line for lat : ', lat)
    print(' --> lshell : ', lshell)
    if (lshell < rmax):
        rmax = lshell
        print(' --> lshell is limiter!')
    else:
        print(' --> rmax is limiter : ', rmax)

    r1 = rmin**(1.0/3.0)
    r2 = rmax**(1.0/3.0)
    r3 = np.linspace(r1, r2, num = nPts)
    r = r3**3
    return r

#------------------------------------------------------------------------
# Main code is here:
#------------------------------------------------------------------------

nPts = 20
rEarth = 6372.0
altMin = 100.0
altMax = 3*rEarth

rMin = (rEarth + altMin) / rEarth
rMax = (rEarth + altMax) / rEarth

fig = plt.figure(figsize = (10,10))
ax = fig.add_axes([0.1,0.1,0.8,0.8])


baseLats = np.arange(-80, 90, 10)

for baseLat in baseLats:

    lat = baseLat * np.pi / 180.0
    lshell = get_lshell_sphere(lat, 1.0)
    print('lshell : ', lshell)
    r = get_r3_spacing(lat, rMin, rMax, nPts)

    lats = get_lat_from_r_and_lshell(r, lshell)
    if (baseLat < 0):
        lats = -lats

    x = r * np.cos(lats)
    z = r * np.sin(lats)

    ax.plot(x, z)

ax.set_ylim([-rMax*1.1, rMax*1.1])
ax.set_xlim([-rMax*1.1, rMax*1.1])
ax.set_aspect(1.0)

plotfile = 'test.png'
print('Writing plot : ', plotfile)
fig.savefig(plotfile)
plt.close()
