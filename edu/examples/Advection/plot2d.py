#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import argparse

# ----------------------------------------------------------------------
# Function to parse input arguments
# ----------------------------------------------------------------------

def parse_args():

    parser = argparse.ArgumentParser(description = 'Plot advection model results')
    parser.add_argument('filein', metavar = 'filein', nargs = 1, \
                        help = 'input file')

    args = parser.parse_args()

    return args

def read_file(file):

    fpin = open(file, 'r')
    lines = fpin.readlines()
    aline = lines[0].split()
    nX = int(aline[0])
    nY = int(aline[1])

    nLines = len(lines)-1
    nTimes = int(nLines / nX)
    print(nX, nY, nLines, nTimes)
    
    values = np.zeros((nTimes, nX, nY))

    for iTime in range(nTimes):
        for i in range(nX):
            iL = iTime * nX + i + 1
            line = lines[iL]
            aline = line.split()
            v = np.array(aline).astype(float)
            values[iTime][i][:] = v

    return values

def plot_data(values, fileout, mini, maxi, cmap):

    fig = plt.figure(figsize = (10,10))
    ax = fig.add_subplot(111)
    pc = ax.pcolor(v, cmap=cmap, vmin = mini, vmax = maxi)
    plt.colorbar(pc)

    fig.savefig(fileout)
    fig.clf()
    plt.close()
        

args = parse_args()
filein = args.filein
if (not np.isscalar(filein)):
    filein = filein[0]

values = read_file(filein)

n = len(values)
if (n > 20):
    nSkip = int(n / 10)
else:
    nSkip = 1

if (nSkip > 10):
    nSkip = 10

mini = np.nanmin(values)
maxi = np.nanmax(abs(values))

if (mini < 0):
    mini = -maxi
    cmap = 'bwr'
else:
    cmap = 'plasma'

print(mini,maxi)

for i in range(0, n, nSkip):
    print("Plotting time : ", i, " of ", n)
    v = values[i].transpose()
    fileout = filein+'.%04d.png' % i
    plot_data(v, fileout, mini, maxi, cmap=cmap)
    


