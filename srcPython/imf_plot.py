#!/usr/bin/env python

from swmf_imf import *
import argparse
import numpy as np


# ----------------------------------------------------------------------
# Function to parse input arguments
# ----------------------------------------------------------------------

def parse_args():

    # If you use nargs = 1, it produces a list.  If not included, it
    # produces a scalar
    
    parser = argparse.ArgumentParser(description = 'Plot IMF data')
    parser.add_argument('infile', metavar = 'infile', \
                        help = 'File to read and rotate')
    parser.add_argument('-start', metavar = 'start', \
                        help = 'start time relative to start of file (hours)', \
                        default = 0, type = float)
    parser.add_argument('-end', metavar = 'end', \
                        help = 'end time relative to start of file (hours)', \
                        default = -1, type = float)
    parser.add_argument('-plotfile', metavar = 'plotfile', \
                        help = 'plot file name', default = '')

    args = parser.parse_args()

    return args

# ----------------------------------------------------------------------
# main code
# ----------------------------------------------------------------------

args = parse_args()
data = read_swmf_file(args.infile)

# calculate time since the start of the time (in hours):

time_hrs = []

for t in data['times']:
    time_hrs.append((t - data['times'][0]).total_seconds()/3600.0)

time_hrs = np.array(time_hrs)

start = args.start
if (args.end < 0):
    end = np.max(time_hrs)
else:
    end = args.end

print(data['Vars'])

data_new = {}
data_new['Vars'] = data['Vars']
data_new['nVars'] = data['nVars']

for v in data['Vars']:
    temp = np.array(data[v])
    temp_sub = temp[((time_hrs >= start) & (time_hrs <= end))]
    data_new[v] = temp_sub

t = np.array(data['times'])
ts = t[((time_hrs >= start) & (time_hrs <= end))]
data_new['times'] = ts

    
if (len(args.plotfile) < 1):
    plotfile = args.infile+'.png'
else:
    plotfile = args.plotfile
plot_imf(data_new, plotfile = plotfile)
