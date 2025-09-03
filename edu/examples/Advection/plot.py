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
    parser.add_argument('fileout', metavar = 'fileout', nargs = 1, \
                        help = 'output file name')

    args = parser.parse_args()

    return args

def read_file(file):

    fpin = open(file, 'r')

    values = []
    for line in fpin:
        aline = line.split()
        v = np.array(aline).astype(float)

        values.append(v)

    return values

def plot_data(values, fileout, x = []):

    if (len(x) == 0):
        nPts = len(values[0])
        x = np.arange(0, nPts)
    fig = plt.figure(figsize = (10,10))
    ax = fig.add_subplot(111)
    n = len(values)-1
    if (n > 20):
        nSkip = int(n / 10)
    else:
        nSkip = 1
        if (n == 0):
            n = 1
    for i, v in enumerate(values):
        per = 0.1 + 0.9 * float(i) / float(n)
        if (i % nSkip == 0):
            print('Plotting... ', i, values[i][10])
            ax.plot(x, values[i])
    fig.savefig(fileout)
    plt.close()
        

args = parse_args()
filein = args.filein
if (not np.isscalar(filein)):
    filein = filein[0]
fileout = args.fileout
if (not np.isscalar(fileout)):
    fileout = fileout[0]

values = read_file(filein)
x = read_file('x.txt')
x = x[0]
plot_data(values, fileout, x = x)

