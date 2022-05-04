#!/usr/bin/env python
""" Standard model visualization routines
"""

from glob import glob
import re
from aetherpy.io import read_routines
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
from netCDF4 import Dataset
from aetherpy.utils.time_conversion import datetime_to_epoch

#----------------------------------------------------------------------------
# This returns the core of the filename without the _g????.nc
#----------------------------------------------------------------------------

def get_core_file(filename):
    coreFile = ''
    m = re.match('.*([0123]D.*)(_g\d*)(\..*)',filename)
    if m:
        coreFile = m.group(1)            
    return coreFile

#----------------------------------------------------------------------------
# Add to the list of strings if there isn't already an identical string
#----------------------------------------------------------------------------

def append_if_unique(list, newItem):
    newList = list
    IsFound = False
    for item in list:
        if (item == newItem):
            IsFound = True
    if (not IsFound):
        newList.append(newItem)
    return newList

#----------------------------------------------------------------------------
# This looks at all of the netcdf files, figures out the core name, and
# adds them to the list of files to process
#----------------------------------------------------------------------------

def get_base_files():
    filelist = sorted(glob('?????*.nc'))
    files = []
    for file in filelist:
        coreFile = get_core_file(file)
        if (len(coreFile) > 0):
            files = append_if_unique(files, coreFile)
    return files

#----------------------------------------------------------------------------
# Simply return the index of the matching string from a list
#----------------------------------------------------------------------------

def find_var_index(allVars, varToFind):
    index = -1
    for i,var in enumerate(allVars):
        if (var == varToFind):
            index = i
            break
    return index

#----------------------------------------------------------------------------
# Find the minimum and maximum of a variable in a bunch of blocks
#----------------------------------------------------------------------------

def determine_min_max(allBlockData, varToPlot, altToPlot):
    mini = 1e32
    maxi = -1e32
    for data in allBlockData:
        iVar = find_var_index(data['vars'], varToPlot)
        if (iVar > -1):
            v = data[iVar][2:-2, 2:-2, altToPlot].transpose()
            if (np.min(v) < mini):
                mini = np.min(v)
            if (np.min(v) > maxi):
                maxi = np.max(v)
    if (mini < 0):
        if (np.abs(mini) > maxi):
            maxi = np.abs(mini)
        mini = -maxi
    return mini, maxi
    
#----------------------------------------------------------------------------
# Plot a single block of a variable at a given altitude
#----------------------------------------------------------------------------

def plot_block(data, varToPlot, altToPlot, ax, mini, maxi):
    iLon = find_var_index(data['vars'], 'lon')
    iLat = find_var_index(data['vars'], 'lat')
    iAlt = find_var_index(data['vars'], 'z')
    iVar = find_var_index(data['vars'], varToPlot)
        
    alts = data[iAlt][0][0] / 1000.0  # Convert from m to km
    lons = data[iLon][2:-2, 0, 0]
    lats = data[iLat][0, 2:-2, 0]

    x_pos = lons
    y_pos = lats
    v = data[iVar][2:-2, 2:-2, altToPlot].transpose()
    dx = (x_pos[1] - x_pos[0]) / 2.0
    xp = np.append(x_pos - dx, x_pos[-1:] + dx)
    dy = (y_pos[1] - y_pos[0]) / 2.0
    yp = np.append(y_pos - dy, y_pos[-1] + dy)

    if (mini < 0):
        cmap = cm.bwr
    else:
        cmap = cm.plasma
    
    ax.pcolormesh(xp, yp, v, vmin = mini, vmax = maxi, cmap = cmap)

#----------------------------------------------------------------------------
# make a figure with all of the block plotted for the specified variable
# and altitude
#----------------------------------------------------------------------------

def plot_all_blocks(allBlockData, varToPlot, altToPlot, plotFile):
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)
    mini, maxi = determine_min_max(allBlockData, varToPlot, altToPlot)
    for data in allBlockData:
        plot_block(data, varToPlot, altToPlot, ax, mini, maxi)
    print('  Outputting file : ', plotFile)
    fig.savefig(plotFile)
    plt.close(fig)
    
#----------------------------------------------------------------------------
# Read all of the block files in, given the core filename
#----------------------------------------------------------------------------

def read_block_files(coreFile):
    print('Figuring out coreFile : ', coreFile)
    fileList = sorted(glob(coreFile + '_g*.nc'))
    header = read_routines.read_aether_headers(fileList)
    allBlockData = []
    for file in fileList:
        print('  Reading file : ', file)
        data = read_routines.read_aether_file(file, file_vars=header['vars'])
        allBlockData.append(data)
    return allBlockData

#----------------------------------------------------------------------------
# return the nLons (nX), nLats (nY), and nAlts (nZ) of a block's data
#----------------------------------------------------------------------------

def get_sizes(allBlockData):
    lon3d = np.asarray(allBlockData[0][0])
    s = lon3d.shape
    nX = s[0]
    nY = s[1]
    nZ = s[2]
    return nX, nY, nZ
    
#----------------------------------------------------------------------------
# Write a NetCDF file from the data
#----------------------------------------------------------------------------

def write_netcdf(allBlockData, fileName):

    print('  Outputting file : ', fileName)
    ncfile = Dataset(fileName, 'w')

    nBlocks = len(allBlockData)
    nLons, nLats, nZ = get_sizes(allBlockData)

    lon_dim = ncfile.createDimension('lon', nLons)
    lat_dim = ncfile.createDimension('lat', nLats)
    z_dim = ncfile.createDimension('z', nZ)
    block_dim = ncfile.createDimension('block', None)
    time_dim = ncfile.createDimension('time', 1)

    t = ncfile.createVariable('time', np.float64, ('time',))
    
    oneBlock = allBlockData[0]

    t = datetime_to_epoch(oneBlock["time"])
        
    allNetCDFVars = []
    # create all of the variables
    varList = []
    for iV, v in enumerate(oneBlock['vars']):
        varList.append(v)
        longName = oneBlock['long_name'][iV]
        unitName = oneBlock['units'][iV]
        allNetCDFVars.append(ncfile.createVariable(v, np.float32, \
                                                   ('block', 'lon', 'lat', 'z')))
        allNetCDFVars[-1].units = unitName
        allNetCDFVars[-1].long_name = longName

    # loop through variables now
    for iV, v in enumerate(varList):
        for iB, oneBlock in enumerate(allBlockData):
            tmp = np.asarray(oneBlock[iV])
            allNetCDFVars[iB][:,:,:] = tmp
        
    ncfile.close()

#----------------------------------------------------------------------------
# main code
#----------------------------------------------------------------------------
        
files = get_base_files()

for coreFile in files:
    allBlockData = read_block_files(coreFile)
    netcdfName = coreFile + '.nc'
    write_netcdf(allBlockData, netcdfName)

    plotFile = coreFile + '.png'
    var = allBlockData[0]['vars'][3]
    plot_all_blocks(allBlockData, var, 10, plotFile)
