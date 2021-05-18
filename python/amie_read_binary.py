#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
from pylab import cm
from amie_routines import *
import sys

args = sys.argv

file = args[1]
data = amie_read_binary(file)

lats = data["lats"]
mlts = data["mlts"]
vars = data["Vars"]

theta, r = np.meshgrid(mlts * np.pi/12.0 - np.pi/2.0, 90.0 - lats)

iT = 0
pot2d = data[vars[0]][iT]
eflux2d = data[vars[1]][iT]

maxi = np.max(np.abs(eflux2d))
mini = 0.0

fig = plt.figure(figsize = (10,10))
ax = fig.add_subplot(projection = 'polar')

norm = cm.colors.Normalize(vmax=mini, vmin=maxi)
if (mini >= 0):
    cmap = cm.plasma
else:
    cmap = cm.bwr

cax = ax.pcolor(theta, r, eflux2d, \
                vmin = mini, vmax = maxi, cmap = cmap)
ax.contour(theta, r, pot2d, 6, colors = 'k')
cbar = fig.colorbar(cax)

i = file.find('.bin')
outfile = file[0:i]+".png"

fig.savefig(outfile)
plt.close()


