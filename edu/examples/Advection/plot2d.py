#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import argparse
import matplotlib.cm as cm

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

def plot_data(values, fileout, mini, maxi):

    fig = plt.figure(figsize = (10,10))
    ax = fig.add_subplot(111)
    cmap = cm.plasma if mini >= 0 else cm.bwr
    ax.pcolor(v, vmin = mini, vmax = maxi, cmap = cmap)

    fig.savefig(fileout)
    plt.close()
        

args = parse_args()
filein = args.filein
if (not np.isscalar(filein)):
    filein = filein[0]

values = read_file(filein)

n = len(values)
if (n > 20):
    nSkip = int(n / 5)
else:
    nSkip = 1

if (nSkip > 10):
    nSkip = 5

mini = np.min(values)
maxi = np.max(abs(values))

if (mini < 0):
    mini = -maxi

print(mini,maxi)

for i in range(0, n, nSkip):
    print("Plotting time : ", i, " of ", n)
    v = values[i].transpose()
    fileout = filein+'.%04d.png' % i
    plot_data(v, fileout, mini, maxi)
    


