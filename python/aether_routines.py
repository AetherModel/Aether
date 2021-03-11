#!/usr/bin/env python

from netCDF4 import Dataset
import datetime as dt
import numpy as np

def read_aether_header(filelist):
    r""" Grab ancillary information from the netcdf file

    Get some ancillary information from the netcdf file.  This simply
    mirrors the gitm header info.

    Parameters
    ----------
    filelist : list
        A list of netcdf file names.  The number of files is recorded, but only the last
        file is used.

    Returns
    -------
    header: A dictionary containing information about the netcdf file, such
            as nLons, nLons, nAlts, nVars, variable names, time(s)
    """

    header = {"nFiles": len(filelist), \
              "version": 0.1, \
              "nLons": 0, \
              "nLats": 0, \
              "nAlts": 0, \
              "nVars": 0, \
              "vars": [], \
              "time": [], \
              "filename": [filelist[-1]] }
    
    ncfile = Dataset(filelist[-1], 'r')
    header["nVars"] = 0
    for var in ncfile.variables.values():
        s = var.shape
        if (len(s) == 3):
            header["nLons"] = s[0]
            header["nLats"] = s[1]
            header["nAlts"] = s[2]
            header["vars"].append(var.name)
            header["nVars"] += 1

    base = dt.datetime(1965, 1, 1, 0, 0, 0)
    time = np.array(ncfile.variables['Time'])
    current = base + dt.timedelta(seconds = time[0])
    header["time"].append(current)
            
    ncfile.close()

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
    

