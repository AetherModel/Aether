#!/usr/bin/env python

from netCDF4 import Dataset
import datetime as dt
import numpy as np

def read_aether_header(filelist):
    """ Obtain ancillary information from the netCDF file

    Get some ancillary information from the netcdf file.  This simply
    mirrors the gitm header info.

    Parameters
    ----------
    filelist : list
        A list of netcdf file names.  The number of files is recorded, but only the last
        file is used.

    Returns
    -------
    header : dict
        A dictionary containing header information from the netCDF files, including:
        nFilles - number of files in list
        nLons - number of longitude grids
        nLats - number of latitude grids
        nAlts - number of altitude grids
        nVars - number of data variable names
        time - list of datetimes with data
    """

    header = {"nFiles": len(filelist),
              "version": 0.1, \
              "nLons": 0, \
              "nLats": 0, \
              "nAlts": 0, \
              "nVars": 0, \
              "vars": [], \
              "time": [], \
              "filename": [filelist[-1]] }
    
    with Dataset(filelist[-1], 'r') as ncfile:
        header["nVars"] = 0
        for var in ncfile.variables.values():
            if (len(var.shape) == 3):
                header["nLons"] = var.shape[0]
                header["nLats"] = var.shape[1]
                header["nAlts"] = var.shape[2]
                header["vars"].append(var.name)
                header["nVars"] += 1

        base = dt.datetime(1965, 1, 1, 0, 0, 0)
        time = np.array(ncfile.variables['Time'])
        current = base + dt.timedelta(seconds=time[0])
        header["time"].append(current)

    return header
    
#-----------------------------------------------------------------------------
# Read in some variables from one file. Inputs:
#   - file: netcdf file to read
#   - vars: list of variable NAMES to read
# Output:
#   - data["time"]: datetime of the file
#   - data[NUMBER]: data that is read in.
#                   NUMBER does from 0 - number of vars read in (0-3 typical)
#-----------------------------------------------------------------------------

def read_aether_one_file(file, vars):
    r""" Read in list of variables from a single netCDF file

    Parameters
    ----------
    file: netcdf file to read
    vars: list of variable NAMES to read

    Returns
    -------
    data["time"]: datetime of the file
    data[NUMBER]: data that is read in.
                  NUMBER goes from 0 - number of vars read in (0-3 typical)

    """

    data = {}
    
    ncfile = Dataset(file, 'r')

    base = dt.datetime(1965, 1, 1, 0, 0, 0)
    time = np.array(ncfile.variables['Time'])
    current = base + dt.timedelta(seconds = time[0])
    data["time"] = current

    for iVar, var in enumerate(vars):
        data[iVar] = np.array(ncfile.variables[var])
    ncfile.close()

    return data
    
