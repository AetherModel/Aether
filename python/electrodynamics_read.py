#!/usr/bin/env python

import netCDF4 as nc
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
from pylab import cm
from datetime import datetime
from datetime import timedelta
import sys

args = sys.argv

file = args[1]
ncfile = nc.Dataset(file, 'r')

lats = np.array(ncfile.variables['MagneticLatitude'])
mlts = np.array(ncfile.variables['MagneticLocalTime'])
potential = np.array(ncfile.variables['Potential'])
eflux = np.array(ncfile.variables['EFlux'])
times_in_s = np.array(ncfile.variables['Time'])

nTimes = len(times_in_s)
base = datetime(1965,1,1,0,0,0)
times = []
for t in times_in_s:
    times.append(base + timedelta(seconds = t))

theta, r = np.meshgrid(mlts * np.pi/12.0 - np.pi/2.0, 90.0 - lats)

maxi = np.max(np.abs(eflux))
mini = 0.0

fig = plt.figure(figsize = (10,10))
ax = fig.add_subplot(projection = 'polar')

iT = 20
pot2d = potential[iT]
eflux2d = eflux[iT]

norm = cm.colors.Normalize(vmax=mini, vmin=maxi)
if (mini >= 0):
    cmap = cm.plasma
else:
    cmap = cm.bwr

cax = ax.pcolor(theta, r, eflux2d, \
                vmin = mini, vmax = maxi, cmap = cmap)
ax.contour(theta, r, pot2d, 6, colors = 'k')
cbar = fig.colorbar(cax)

fig.savefig('netcdf_test.png')
plt.close()
