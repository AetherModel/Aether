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
            v = data[iVar][2:-2, 2:-2, altToPlot]
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

def plot_block(data, varToPlot, altToPlot, ax, mini, maxi, i):
    iLon = find_var_index(data['vars'], 'lon')
    iLat = find_var_index(data['vars'], 'lat')
    iAlt = find_var_index(data['vars'], 'z')
    iVar = find_var_index(data['vars'], varToPlot)

    if (mini < 0):
        cmap = cm.bwr
    else:
        cmap = cm.plasma
    
    alts = data[iAlt][0][0] / 1000.0  # Convert from m to km

    # Change to 2d representation:

    #lons = data[iLon][2:-2, 2:-2, 0]
    #lats = data[iLat][2:-2, 2:-2, 0]
    #v = data[iVar][2:-2, 2:-2, altToPlot]
    lons = data[iLon][:, :, 0]
    lats = data[iLat][:, :, 0]
    v = data[iVar][:, :, altToPlot]

    #print("lons : ", lons[:, 5])
    #print("lats : ", lats[:, 5])
    #print("value : ", i, v[:, 5]*180.0/np.pi)
    
    nLons, nLats = np.shape(lons)
    xp = np.zeros((nLons+1, nLats+1))
    yp = np.zeros((nLons+1, nLats+1))

    for iX in np.arange(1, nLons):
        for iY in np.arange(1, nLats):
            xp[iX, iY] = (lons[iX-1, iY] + lons[iX, iY])/2.0
            yp[iX, iY] = (lats[iX, iY-1] + lats[iX, iY])/2.0

        xp[iX, 0] = (lons[iX-1, 0] + lons[iX, 0])/2.0
        xp[iX, nLats] = xp[iX, nLats-1]

        yp[iX, 0] = 2 * lats[iX, 0] - yp[iX, 0]
        yp[iX, nLats] = 2 * yp[iX, nLats-1] - lats[iX, nLats-1]
    
    # more xp
    for iY in np.arange(1, nLats):
        xp[0, iY] = 2 * lons[0, iY] - xp[1, iY]
        xp[nLons, iY] = 2 * xp[nLons-1, iY] - lons[nLons-1, iY]
    xp[0, 0] = xp[0, 1]
    xp[nLons, nLats] = xp[nLons, nLats-1]

    xp[0, nLats] = xp[0, nLats-1]
    xp[nLons, 0] = 2 * xp[nLons-1, 0] - lons[nLons-1, 0]
            
    # more yp
    for iY in np.arange(1, nLats):
        yp[0, iY] = (lats[0, iY-1] + lats[0, iY])/2.0
        yp[nLons, iY] = (lats[nLons-1, iY-1] + lats[nLons-1, iY])/2.0
    yp[0, 0] = yp[1, 0]
    yp[nLons, nLats] = yp[nLons-1, nLats]

    yp[0, nLats] = 2 * yp[0, nLats-1] - lats[0, nLats-1]
    yp[nLons, 0] = yp[nLons-1, 0]

    
    #print(yp)

    #ax.pcolor(xp, yp, v, vmin = mini, vmax = maxi, cmap = cmap)

    #n = np.asarray((v - mini) / (maxi - mini) * 255.0, dtype = np.int)
    ax.scatter(lons, lats, c = v, cmap=cmap, vmin = mini, vmax = maxi)

    #lons = data[iLon][2:-2, 0, 0]
    #lats = data[iLat][0, 2:-2, 0]
    #
    #x_pos = lons
    #y_pos = lats
    #v = data[iVar][2:-2, 2:-2, altToPlot].transpose()
    #dx = (x_pos[1] - x_pos[0]) / 2.0
    #xp = np.append(x_pos - dx, x_pos[-1:] + dx)
    #dy = (y_pos[1] - y_pos[0]) / 2.0
    #yp = np.append(y_pos - dy, y_pos[-1] + dy)
    #
    #ax.pcolormesh(xp, yp, v, vmin = mini, vmax = maxi, cmap = cmap)

#----------------------------------------------------------------------------
# make a figure with all of the block plotted for the specified variable
# and altitude
#----------------------------------------------------------------------------

def plot_all_blocks(allBlockData, varToPlot, altToPlot, plotFile):
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)
    mini, maxi = determine_min_max(allBlockData, varToPlot, altToPlot)
    i = 0
    for data in allBlockData:
        if (i < 100):
            plot_block(data, varToPlot, altToPlot, ax, mini, maxi, i)
        i = i+1
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
    #print(allBlockData[0]['vars'])
    var = allBlockData[0]['vars'][14]
    plot_all_blocks(allBlockData, var, 10, plotFile)
    
