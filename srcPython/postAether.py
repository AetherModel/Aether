#!/usr/bin/env python
""" Standard model visualization routines
"""

import datetime as dt
from glob import glob
import re
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
from netCDF4 import Dataset
import argparse
import os
import json
from struct import unpack

# ----------------------------------------------------------------------
# Want to eliminate need for aetherpy to be installed for ease of use.
# Therefore a bunch of this stuff is from aetherpy.
# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
# convert time from datetime to double

def datetime_to_epoch(dtime):
    """Convert datetime to epoch seconds.

    Parameters
    ----------
    dtime : dt.datetime or dt.date
        Datetime object

    Returns
    -------
    epoch_time : float
        Seconds since 1 Jan 1965

    """

    epoch_time = (dtime - dt.datetime(1965, 1, 1)).total_seconds()

    return epoch_time

# ----------------------------------------------------------------------
# convert time format from double to datetime

def epoch_to_datetime(epoch_time):
    """Convert from epoch seconds to datetime.

    Parameters
    ----------
    epoch_time : int
        Seconds since 1 Jan 1965

    Returns
    -------
    dtime : dt.datetime
        Datetime object corresponding to `epoch_time`

    Notes
    -----
    Epoch starts at 1 Jan 1965.

    """

    dtime = dt.datetime(1965, 1, 1) + dt.timedelta(seconds=epoch_time)

    return dtime


class DataArray(np.ndarray):
    def __new__(cls, input_array, attrs={}):
        obj = np.asarray(input_array).view(cls)
        obj.attrs = attrs
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.attrs = getattr(obj, 'attrs', {
            'units': None,
            'long_name': None
        })

# ----------------------------------------------------------------------
# read aether header - figure out whether it is json or netcdf
        
def read_aether_headers(filelist, finds=None, filetype="netcdf"):
    """Obtain ancillary information from Aether files.

    Parameters
    ----------
    filelist : array-like
        Array-like object of names for Aether files
    finds : int, list, slice, or NoneType
        Index(es) for file(s) from which the header information will be read.
        If None, reads headers from all files, ensuring data across all files
        are consistent. (default=None)
    filetype : str
        Input file type (must be the same for all files). Expects one of
        'ascii' or 'netcdf' (default='netcdf')

    Returns
    -------
    header : dict
        A dictionary containing header information from the files,
        including:
        nlons - number of longitude grids
        nlats - number of latitude grids
        nalts - number of altitude grids
        vars - list of data variable names
        time - list of datetimes for the processed file start times
        filename - list of the input filenames

    Raises
    ------
    ValueError
        If there are an unexpected number or name of variables or dimensions
        for any file in the file list.

    See Also
    --------
    read_aether_ascii_header and read_aether_netcdf_header

    Notes
    -----
    Reading may be sped up by loading data from a single header, but loading
    from all will ensure all files in the list have the same format.

    """

    # Initialize the output
    header = {"vars": [], "time": [], "filename": np.asarray(filelist)}

    # Ensure the filelist is array-like, allowing slicing of input
    if header['filename'].shape == ():
        header['filename'] = [filelist]
    else:
        header['filename'] = list(header['filename'])

    # Ensure selected files are list-like
    if finds is None:
        hfiles = np.asarray(header['filename'])
    else:
        hfiles = np.asarray(header['filename'])[finds]
        if hfiles.shape == ():
            hfiles = np.asarray([hfiles])

    # Read the header info from the desired files
    for iFile, filename in enumerate(hfiles):
        if filetype.lower() == "netcdf":
            fheader = read_aether_netcdf_header(filename)
        else:
            fheader = read_aether_json_header(filename)
            header['filename'][iFile] = filename.replace(".json", ".bin")

        for hkey in fheader.keys():
            if hkey in header.keys():
                if hkey == 'time':
                    if fheader[hkey] not in header[hkey]:
                        header[hkey].append(fheader[hkey])
                elif hkey == 'vars' and np.any(header[hkey] != fheader[hkey]):
                    #header[hkey] = list(np.unique(header[hkey] + fheader[hkey]))
                    header[hkey] = fheader[hkey]
                elif hkey != 'filename' and header[hkey] != fheader[hkey]:
                    raise ValueError(''.join(['header values changed between',
                                              ' files: {:} != {:}'.format(
                                                  header[hkey],
                                                  fheader[hkey])]))
            else:
                header[hkey] = fheader[hkey]

    return header

# ----------------------------------------------------------------------
# read netcdf header

def read_aether_netcdf_header(filename, epoch_name='time'):
    """Read header information from an Aether netCDF file.

    Parameters
    ----------
    filename : str
        An Aether netCDF filename
    epoch_name : str
        Epoch variable name used in the data file (default='time')

    Returns
    -------
    header : dict
        A dictionary containing header information from the netCDF files,
        including:
        nlons - number of longitude grids
        nlats - number of latitude grids
        nalts - number of altitude grids
        vars - list of data variable names
        time - list of datetimes with data
        filename - filename of file containing header data

    See Also
    --------
    read_aether_header

    """

    # Test the filename and add it to the header
    if not os.path.isfile(filename):
        raise IOError('unknown aether netCDF file: {:}'.format(filename))

    header = {'filename': filename,
              'nlons': 0,
              'nlats': 0,
              'nalts': 0,
              'nblocks': 1}

    # Open the file and read the header data
    with Dataset(filename, 'r') as ncfile:
        ncvars = list()
        IsFirst = True
        for var in ncfile.variables.values():
            if (len(var.shape) >= 3):
                iOff_ = 0
                if (len(var.shape) == 4):
                    iOff_ = 1
                    header["nblocks"] = var.shape[0]
                nlons = var.shape[0 + iOff_]
                nlats = var.shape[1 + iOff_]
                nalts = var.shape[2 + iOff_]
                ncvars.append(var.name)

                if (IsFirst):
                    header["nlons"] = nlons
                    header["nlats"] = nlats
                    header["nalts"] = nalts
                    IsFirst = False

        # Save the unique variable names
        ncvars = np.unique(ncvars)
        if "vars" not in header.keys() or len(header['vars']) == 0:
            header["vars"] = list(ncvars)
        elif np.any(header["vars"] != ncvars):
            raise IOError(''.join(['unexpected number or name of variables in',
                                   ' file: ', filename]))

        # Add the time for this file
        epoch = np.double(ncfile.variables[epoch_name][0])
        header["time"] = epoch_to_datetime(epoch)

    return header

# ----------------------------------------------------------------------
# read json header

def read_aether_json_header(filename):
    """Read information from an Aether ascii header file.

    Parameters
    ----------
    filename : str
        An ascii header file name or binary filename with associated header
        file present in the same directory.

    Returns
    -------
    header : dict
        A dictionary containing header information from the netCDF files,
        including:
        nlons - number of longitude grids
        nlats - number of latitude grids
        nalts - number of altitude grids
        nvars - number of data variable names
        time - datetime for the file
        vars - variables in the file
        long_name - variables in the file
        units - units of variables in the file
        version - version number of the file
        filename - filename of file containing header data

    See Also
    --------
    read_aether_header

    """

    header = {"filename": filename}

    # returns JSON object as 
    # a dictionary
    fpin = open(filename)
    jsonHeader = json.load(fpin)
    fpin.close()

    # future proof:
    if ('long_name' in jsonHeader.keys()):
        needLongName = False
    else:
        needLongName = True
        
    for key in jsonHeader.keys():
        newkey = key.lower()
        if (newkey == 'variables'):
            newkey = 'vars'
        if (newkey == 'time'):
            t = jsonHeader['time']
            time = dt.datetime(t[0], t[1], t[2],
                               t[3], t[4], t[5], t[6])
            jsonHeader['time'] = time
        header[newkey] = jsonHeader[key]
        if ((newkey == 'vars') and needLongName):
            header['long_name'] = jsonHeader[key]
            
    return header
    
# ----------------------------------------------------------------------
# read aether file that is in binary format

def read_aether_one_binary_file(header, ifile, vars_to_read):
    """Read in list of variables from a single netCDF file.

    Parameters
    ----------
    header : dict
        Dict of header data returned from read_aether_ascii_header
    ifile : int
        Integer corresponding to filename in the header dict
    vars_to_read : list
        List of desired variable names to read

    Returns
    -------
    data : dict
        Dict with keys 'time', which contains a datetime object specifying the
        time of the file and indices ranging from 0 to `len(vars_to_read) - 1`,
        corresponding to the variable names in `vars_to_read` that holds arrays
        of the specified data

    See Also
    --------
    read_aether_ascii_header

    """

    file_to_read = header["filename"][ifile]
    print("Reading file : ", file_to_read)

    data = {hkey: header[hkey] for hkey in ["version",
                                            "nlons", "nlats",
                                            "nalts", "nvars", "vars",
                                            "units", "long_name"]}
    data["time"] = header["time"][0]
    with open(file_to_read, 'rb') as fin:
        # Read in the header data
        header_len = 0
        num_tot = data["nlons"] * data["nlats"] * data["nalts"]

        # Data were saved as floats, not doubles
        data_len = num_tot * 4
        endChar = '<'
        for ivar in vars_to_read:
            fin.seek(header_len + ivar * data_len)
            data[ivar] = np.array(unpack(endChar + '%if' % num_tot,
                                         fin.read(data_len)))
            data[ivar] = data[ivar].reshape(
                (data["nlons"], data["nlats"], data["nalts"]), order="F")

    return data

# ----------------------------------------------------------------------
# read aether netcdf file

def read_aether_file(filename, file_vars=None, epoch_name='time'):
    """Read in list of variables from a netCDF file.

    Parameters
    ----------
    filename : str
        Name of netCDF file to read
    file_vars : list or NoneType
        List of desired variable names to read, or None to read all
        (default=None)
    epoch_name : str
        Epoch variable name (default='time')

    Returns
    -------
    data : dict or xr.Dataset
        If `out_type` is 'dict', `data` is a dict with keys 'time', which
        contains a datetime object specifying the time of the file, 'vars',
        which contains a list of variable names, and zero-offset indices,
        corresponding to the variable names in 'vars' key, and 'units',
        a list of unit strings for each of the variables.

    """
    data = dict()
    print("Reading file : ", filename)

    if not os.path.isfile(filename):
        raise IOError('input file does not exist')

    # Set the default attributes
    def_attr = {'units': '', 'long_name': None}

    with Dataset(filename, 'r') as ncfile:
        # Save a list of all desired variable names
        data['vars'] = [var for var in ncfile.variables.keys()
                        if file_vars is None or var in file_vars]
        for attr in def_attr.keys():
            data[attr] = list()

        # Save the data as numpy arrays, using variable index as a key
        for i_var, var in enumerate(data['vars']):
            data[i_var] = np.array(ncfile.variables[var])

            # Save the attributes
            for attr in def_attr.keys():
                if hasattr(ncfile.variables[var], attr):
                    data[attr].append(getattr(ncfile.variables[var], attr))
                else:
                    if def_attr[attr] is None:
                        data[attr].append(var)
                    else:
                        data[attr].append(def_attr[attr])

        # Calculate the date and time for this data
        time = np.array(ncfile.variables[epoch_name])
        data['time'] = epoch_to_datetime(time[0])

    return data

# ----------------------------------------------------------------------
# Function to parse input arguments
# ----------------------------------------------------------------------

def parse_args():

    parser = argparse.ArgumentParser(description = 'Post process Aether files')
    parser.add_argument('-rm', \
                        help='removes processed files', \
                        action="store_true")

    args = parser.parse_args()

    return args

#----------------------------------------------------------------------------
# This returns the core of the filename without the _g????.nc
#----------------------------------------------------------------------------

def get_core_file(filename):
    coreFile = ''
    isEnsemble = False
    ensembleFile = ''
    ensembleNumber = -1
    m = re.match('.*([0123]D.*)(_g\d*)(\..*)',filename)
    if m:
        coreFile = m.group(1)
        # check if file is a member of an ensemble:
        check = re.match('.*([0123]D.*)(_m)(\d*)',coreFile)
        if (check):
            ensembleFile = check.group(1)
            isEnsemble = True
            ensembleNumber = int(check.group(3)) + 1

    fileInfo = {'coreFile': coreFile,
                'isEnsemble': isEnsemble,
                'ensembleFile': ensembleFile,
                'ensembleNumber': ensembleNumber,
                'ensembleMembers': -1}
    return fileInfo

#----------------------------------------------------------------------------
# Add to the list of strings if there isn't already an identical string
#----------------------------------------------------------------------------

def if_unique(list, newItem):
    IsFound = False
    index = -1
    for i, item in enumerate(list):
        if (item == newItem):
            IsFound = True
            index = i
    return IsFound, index

#----------------------------------------------------------------------------
# This looks at all of the netcdf files, figures out the core name, and
# adds them to the list of files to process
#----------------------------------------------------------------------------

def get_base_files():
    isNetCDF = False
    filelist = sorted(glob('?????*.bin'))
    if (len(filelist) == 0):
        filelist = sorted(glob('?????*.nc'))
        isNetCDF = True

    print(filelist)

    files = []
    filesInfo = []
    for file in filelist:
        fileInfo = get_core_file(file)
        print(file, '->> ', fileInfo)
        fileInfo["isNetCDF"] = isNetCDF
        coreFile = fileInfo['coreFile']
        if (len(coreFile) > 0):
            IsFound, i = if_unique(files, coreFile)
            if (not IsFound):
                files.append(coreFile)
                filesInfo.append(fileInfo)

    # Figure out ensembles:
    # (1) get list of unique ensemble filess
    # (2) count how many files have this unique ensemble name
    
    ensembleFiles = []
    ensembleCounter = []
    for fileInfo in filesInfo:
        if (len(fileInfo['ensembleFile']) > 0):
            IsFound, i = if_unique(ensembleFiles, fileInfo['ensembleFile'])
            if (IsFound):
                ensembleCounter[-1] += 1
            else:
                ensembleFiles.append(fileInfo['ensembleFile'])
                ensembleCounter.append(1)
    
    # (3) store the number of ensemble members:            
    for i, fileInfo in enumerate(filesInfo):
        if (len(fileInfo['ensembleFile']) > 0):
            IsFound, item = if_unique(ensembleFiles, fileInfo['ensembleFile'])
            if (IsFound):
                filesInfo[i]['ensembleMembers'] = ensembleCounter[item]
    
    return filesInfo


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
    if (iAlt < 0):
        iAlt = find_var_index(data['vars'], 'alt')
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
    
    plot = ax.scatter(lons, lats, c = v, cmap=cmap, vmin = mini, vmax = maxi)
    
    return plot
    
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
            plot = plot_block(data, varToPlot, altToPlot, ax, mini, maxi, i)
        i = i+1
    cbar = fig.colorbar(plot, ax = ax)
    cbar.set_label(varToPlot)
    ax.set_title(varToPlot)
    ax.set_xlabel('Longitude (deg)')
    ax.set_ylabel('Latitude (deg)')
    print('  Outputting file : ', plotFile)
    fig.savefig(plotFile)
    plt.close(fig)
    
#----------------------------------------------------------------------------
# Read all of the block files in, given the core filename
#----------------------------------------------------------------------------

def read_block_files(coreFile, isNetCDF):
    print('Figuring out coreFile : ', coreFile)
    if (isNetCDF):
        fileList = sorted(glob(coreFile + '_g*.nc'))
        filetype = "netcdf"
    else:
        fileList = sorted(glob(coreFile + '_g*.json'))
        filetype = "json"
    header = read_aether_headers(fileList, filetype=filetype)
    allBlockData = []
    for iFile, file in enumerate(fileList):
        print('  Reading file : ', file)
        if (isNetCDF):
            varsToRead = header['vars']
            data = read_aether_file(file, file_vars=varsToRead)
        else:
            varsToRead = range(len(header['vars']))
            data = read_aether_one_binary_file(header,
                                               iFile,
                                               varsToRead)
            
        allBlockData.append(data)
    return allBlockData, fileList

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
    time_dim = ncfile.createDimension('time', None)

    oneBlock = allBlockData[0]

    time_out = ncfile.createVariable('time', np.float64, ('time',))
    time_out[0] = datetime_to_epoch(oneBlock["time"])
        
    allNetCDFVars = []
    # create all of the variables
    varList = []
    for iV, v in enumerate(oneBlock['vars']):
        varList.append(v)
        if ('long_name' in oneBlock):
            longName = oneBlock['long_name'][iV]
        else:
            longName = v
        unitName = oneBlock['units'][iV]
        allNetCDFVars.append(ncfile.createVariable(v, np.float32, \
                                                   ('block', 'lon', 'lat', 'z')))
        allNetCDFVars[-1].units = unitName
        allNetCDFVars[-1].long_name = longName

        for iB, oneBlock in enumerate(allBlockData):
            tmp = np.asarray(oneBlock[iV])
            allNetCDFVars[-1][iB,:,:,:] = tmp
        
    ncfile.close()

#----------------------------------------------------------------------------
# copy block data in one file
#----------------------------------------------------------------------------

def copy_or_add_block_data(allBlockData,
                           oldBlockData = [],
                           factor = 1.0):

    if (len(oldBlockData) > 0):
        doAdd = True
    else:
        doAdd = False

    newBlockData = []

    for ib, oneBlock in enumerate(allBlockData):
        obCopy = {}
        for key in oneBlock.keys():
            if (type(key) == int):
                if (doAdd):
                    obCopy[key] = oldBlockData[ib][key] + \
                        oneBlock[key] * factor
                else:
                    obCopy[key] = oneBlock[key] * factor
            else:
                obCopy[key] = oneBlock[key]
        newBlockData.append(obCopy)
        
    return newBlockData

#----------------------------------------------------------------------------
# copy block data in one file
#----------------------------------------------------------------------------

iCopy_ = 0
iAdd_ = 1
iSub_ = 2
iMult_ = 3
iPower_ = 4

def do_math_on_block_data(blockData1,
                          blockData2 = [],
                          math = iCopy_,
                          factor = 1.0,
                          iLowestVar = 3):

    newBlockData = []

    for ib, oneBlock in enumerate(blockData1):
        obCopy = {}
        for key in oneBlock.keys():
            if (type(key) == int):
                if (key >= iLowestVar):
                    if (math == iCopy_):
                        obCopy[key] = oneBlock[key] * 1.0
                    if (math == iAdd_):
                        obCopy[key] = oneBlock[key] + blockData2[ib][key]
                    if (math == iSub_):
                        obCopy[key] = oneBlock[key] - blockData2[ib][key]
                    if (math == iMult_):
                        obCopy[key] = oneBlock[key] * factor
                    if (math == iPower_):
                        obCopy[key] = oneBlock[key] ** factor
                else:
                    obCopy[key] = oneBlock[key]
            else:
                obCopy[key] = oneBlock[key]                    
        newBlockData.append(obCopy)
    return newBlockData

#----------------------------------------------------------------------------
# Calc Standard Deviation of Ensemble Run
#----------------------------------------------------------------------------

def calc_std_of_ensembles(filesInfo,
                          ensembleIndexList,
                          ensembleMean):
    
    for i, iF in enumerate(ensembleIndexList):
        cF = filesInfo[iF]['coreFile']
        print('---> Going back through corefiles: ', cF)
        allBlockData, filelist = read_block_files(cF)

        # subtract
        diff = do_math_on_block_data(allBlockData,
                                     blockData2 = ensembleMean,
                                     math = iSub_)
        # square
        diffs = do_math_on_block_data(diff,
                                      factor = 2,
                                      math = iPower_)
                
        if (i == 0):
            sums = do_math_on_block_data(diffs,
                                         math = iCopy_)
        else:
            sums = do_math_on_block_data(sums,
                                         blockData2 = diffs,
                                         math = iAdd_)
    factor = 1.0 / fileInfo['ensembleMembers']
    sumsD = do_math_on_block_data(sums,
                                  factor = factor,
                                  math = iMult_)
    stdData = do_math_on_block_data(sumsD,
                                    factor = 0.5,
                                    math = iPower_)

    return stdData


#----------------------------------------------------------------------------
# write and plot data
#----------------------------------------------------------------------------

def write_and_plot_data(dataToWrite,
                        fileStart,
                        fileAddon,
                        iVar,
                        iAlt):

    netcdfFile = fileStart + fileAddon + '.nc'
    write_netcdf(dataToWrite, netcdfFile)

    plotFile = fileStart + fileAddon + '.png'
    var = dataToWrite[0]['vars'][iVar]
    plot_all_blocks(dataToWrite, var, iAlt, plotFile)

    return


#----------------------------------------------------------------------------
# main code
#----------------------------------------------------------------------------

if __name__ == '__main__':  # main code block

    args = parse_args()

    filesInfo = get_base_files()

    print(filesInfo)

    iVar = 3
    iAlt = 10

    for iFile, fileInfo in enumerate(filesInfo):
        coreFile = fileInfo['coreFile']
        isNetCDF = fileInfo['isNetCDF']
        allBlockData, filelist = read_block_files(coreFile, isNetCDF)
        write_and_plot_data(allBlockData, coreFile, '', iVar, iAlt)
    
        if (fileInfo['isEnsemble']):
            factor = 1.0 / float(fileInfo['ensembleMembers'])
            if (fileInfo['ensembleNumber'] == 1):
                ensembleData = copy_or_add_block_data(allBlockData,
                                                      factor = factor)
                ensembleIndexList = [iFile]
            else:
                ensembleData = copy_or_add_block_data(allBlockData,
                                                      oldBlockData = ensembleData,
                                                      factor = factor)
                ensembleIndexList.append(iFile)
            if (fileInfo['ensembleNumber'] == fileInfo['ensembleMembers']):
    
                write_and_plot_data(ensembleData, fileInfo['ensembleFile'],
                                    '_mean', iVar, iAlt)
                
                stdData = calc_std_of_ensembles(filesInfo,
                                                ensembleIndexList,
                                                ensembleData)
                write_and_plot_data(stdData, fileInfo['ensembleFile'],
                                    '_std', iVar, iAlt)
                            
        if (args.rm):
            print('  ---> Removing files...')
            for file in filelist:
                command = 'rm -f '+file
                os.system(command)
    
