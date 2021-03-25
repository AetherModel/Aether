#!/usr/bin/env python
"""Routines to read Aether files
"""

from netCDF4 import Dataset
import numpy as np

from time_routines import epoch_to_datetime


def read_aether_header(filelist, find=-1):
    """ Obtain ancillary information from the netCDF file

    Parameters
    ----------
    filelist : list
        A list of netcdf file names.  The number of files is recorded, but only
        one file is read.
    find : int
        Index indicating which file in the filelist should be read (default=-1)

    Returns
    -------
    header : dict
        A dictionary containing header information from the netCDF files,
        including:
        nfiles - number of files in list
        nlons - number of longitude grids
        nlats - number of latitude grids
        nalts - number of altitude grids
        nvars - number of data variable names
        time - list of datetimes with data

    Notes
    -----
    This routine obtains the same info as the GITM header routine

    """

    header = {"nfiles": len(filelist), "version": 0.1, "nlons": 0, "nlats": 0,
              "nalts": 0, "nvars": 0, "vars": [], "time": [],
              "filename": [filelist[find]]}

    with Dataset(filelist[find], 'r') as ncfile:
        header["nvars"] = 0
        for var in ncfile.variables.values():
            if len(var.shape) == 3:
                header["nlons"] = var.shape[0]
                header["nlats"] = var.shape[1]
                header["nalts"] = var.shape[2]
                header["vars"].append(var.name)
                header["nvars"] += 1

        time = np.array(ncfile.variables['Time'])
        header["time"].append(epoch_to_datetime(time[0]))

    return header


def read_aether_file(filename, file_vars):
    """ Read in list of variables from a netCDF file

    Parameters
    ----------
    filename : str
        Name of netCDF file to read
    file_vars : list
        List of desired variable names to read

    Returns
    -------
    data : dict
        Dict with keys 'time', which contains a datetime object specifying the
        time of the file and zero-offset indices, corresponding to the
        variable names in `file_vars` that holds arrays of the specified data

    """

    with Dataset(filename, 'r') as ncfile:
        # Save the data as numpy arrays, using variable index as a key
        data = {i_var: np.array(ncfile.variables[var])
                for i_var, var in enumerate(file_vars)}

        # Calculate the date and time for this data
        time = np.array(ncfile.variables['Time'])
        data['time'] = epoch_to_datetime(time[0])

    return data
