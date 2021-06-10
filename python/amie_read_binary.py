#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
from pylab import cm
from amie_routines import *
import sys
import argparse

# ----------------------------------------------------------------------
# Function to parse input arguments
# ----------------------------------------------------------------------

def parse_args():

    parser = argparse.ArgumentParser(description = 'Plot AMIE files')
    parser.add_argument('files', metavar = 'file', nargs = '+', \
                        help = 'Files to process')
    parser.add_argument('-start', nargs=1, \
                        help='start iteration to plot', \
                        default=0, type = int)
    parser.add_argument('-end', nargs=1, \
                        help='end iteration to plot (default to max in file)', \
                        default=10000, type = int)
    parser.add_argument('-step', nargs=1, \
                        help='step in iterations to plot', \
                        default=5, type = int)
    parser.add_argument('-vars', \
                        help='list vars in file', \
                        action="store_true")
    
    args = parser.parse_args()

    return args

# ----------------------------------------------------------------------
# Main Code
# ----------------------------------------------------------------------

args = parse_args()

start = args.start
if (not np.isscalar(start)):
    start = start[0]
end = args.end
if (not np.isscalar(end)):
    end = end[0]
di = args.step
if (not np.isscalar(di)):
    di = di[0]

file = args.files[0]
data = amie_read_binary(file)

lats = data["lats"]
mlts = data["mlts"]
vars = data["Vars"]

if (args.vars):
    for i, v in enumerate(vars):
        print(v)
    exit()

theta, r = np.meshgrid(mlts * np.pi/12.0 - np.pi/2.0, 90.0 - lats)

nTimes = len(data["times"])
if (end < nTimes):
    nTimes = end

if (di < 1):
    di = 1

print(start, nTimes, di)
    
ind = np.arange(start, nTimes, di)

eflux3d = np.array(data[vars[1]])
potential3d = np.array(data[vars[0]])

maxi = np.max(np.abs(eflux3d[ind]))
mini = 0.0

potmax = np.max(np.abs(potential3d[ind]))
potmin = -potmax

dl = (potmax-potmin)/15.0
levels = np.arange(potmin, potmax, dl)

iT = 0

for iT in ind:

    pot2d = potential3d[iT]
    eflux2d = eflux3d[iT]

    time = data["times"][iT]
    print(time)
    title = time.strftime('%b %d, %Y %H:%M:%S')
    
    fig = plt.figure(figsize = (10,10))
    ax = fig.add_subplot(projection = 'polar')

    norm = cm.colors.Normalize(vmax=mini, vmin=maxi)
    if (mini >= 0):
        cmap = cm.plasma
    else:
        cmap = cm.bwr

    cax = ax.pcolor(theta, r, eflux2d, \
                    vmin = mini, vmax = maxi, cmap = cmap)
    ax.contour(theta, r, pot2d, levels, colors = 'w')

    xlabels = ['', '12', '18', '00']
    ylabels = ['80', '70', '60', '50']
    ax.set_xticklabels(xlabels)
    ax.set_yticklabels(ylabels)
    ax.grid(linestyle=':', color='black')
    ax.set_xticks(np.arange(0,2*np.pi,np.pi/2))
    ax.set_yticks(np.arange(10,50,10))
    ax.set_title(title)

    cbar = fig.colorbar(cax, shrink = 0.5, pad=0.01)

    i = file.find('.bin')
    outfile = file[0:i] + "_%4.4d.png" % iT

    print("Writing file : ", outfile)
    fig.savefig(outfile)
    plt.close()


