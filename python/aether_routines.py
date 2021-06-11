#!/usr/bin/env python

from netCDF4 import Dataset
import datetime as dt
import numpy as np
import re
from struct import unpack

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
    

def parse_line_into_int_and_string(line):
    x = line.split(" ")
    number = int(x[0])
    string = x[1]
    if (len(x) > 2):
        for s in x[2:]:
            string = string + " " + s
    return [number, string]

def read_aether_ascii_header(filelist):
    """ Obtain ancillary information from the netCDF file

    Get information from the ascii header file.  This simply
    mirrors the gitm header info.

    Parameters
    ----------
    filelist : list
        A list of ascii header file names.  The number of files is recorded, 
        but only the last file is used.

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
              "nBlocksLons": 0, \
              "nBlocksLats": 0, \
              "nBlocksAlts": 0, \
              "nVars": 0, \
              "vars": [], \
              "time": [], \
              "filename": [] }

    for file in filelist:
        
        header["filename"].append(file)
        l = file.find(".bin")
        file = file[0:l]+".header"
        fpin = open(file, 'r')

        for line in fpin:

            m = re.match(r'BLOCKS',line)
            if m:
                header["nBlocksLons"] = int(fpin.readline())
                header["nBlocksLats"] = int(fpin.readline())
                header["nBlocksAlts"] = int(fpin.readline())
            
            m = re.match(r'NUMERICAL',line)
            if m:
                header["nVars"], s = parse_line_into_int_and_string(fpin.readline())
                header["nLons"], s = parse_line_into_int_and_string(fpin.readline())
                header["nLats"], s = parse_line_into_int_and_string(fpin.readline())
                header["nAlts"], s = parse_line_into_int_and_string(fpin.readline())

            m = re.match(r'VERSION',line)
            if m:
                header["version"] = float(fpin.readline())

            m = re.match(r'TIME',line)
            if m:
                year = int(fpin.readline())
                month = int(fpin.readline())
                day = int(fpin.readline())
                hour = int(fpin.readline())
                minute = int(fpin.readline())
                second = int(fpin.readline())
                msec = int(fpin.readline())
                header["time"].append(dt.datetime(year, month, day, \
                                                  hour, minute, second, msec))

            m = re.match(r'VARIABLE',line)
            if ((m) and (len(header["vars"]) < 1)):
                for i in range(header["nVars"]):
                    n, s = parse_line_into_int_and_string(fpin.readline())
                    header["vars"].append(s.strip())

            
        fpin.close()

    return header
    
#--------------------------------------------------------------------------------
# Read binary file:
#--------------------------------------------------------------------------------

def read_aether_one_binary_file(header, iFile, vars_to_read):
    """ Read in list of variables from a single netCDF file

    Parameters
    ----------
    header : dict
        header returned from read_aether_ascii_header
    vars_to_read : list
        List of desired variable names to read

    Returns
    -------
    data : dict
        Dict with keys 'time', which contains a datetime object specifying the time of
        the file and indices ranging from 0 to `len(vars_to_read) - 1`, corresponding to the
        variable names in `vars_to_read` that holds arrays of the specified data
    
    """

    file_to_read = header["filename"][iFile]
    print("Reading file : "+file_to_read)

    data = {"version": header["version"], \
            "nLons": header["nLons"], \
            "nLats": header["nLats"], \
            "nAlts": header["nAlts"], \
            "nVars": header["nVars"], \
            "time": header["time"][iFile], \
            "vars": header["vars"]}

    f=open(file_to_read, 'rb')

    iHeaderLength = 0
    nTotal = data["nLons"]*data["nLats"]*data["nAlts"]
    # floats and not doubles!
    iDataLength = nTotal*4
    endChar='<'
    for iVar in vars_to_read:
        f.seek(iHeaderLength+iVar*iDataLength)
        data[iVar] = np.array(unpack(endChar+'%if'%(nTotal),f.read(iDataLength)))
        data[iVar] = data[iVar].reshape( 
            (data["nLons"],data["nLats"],data["nAlts"]),order="F")
    f.close()
    
    return data


        
#--------------------------------------------------------------------------------
# Read netcdf file:
#--------------------------------------------------------------------------------
    
def read_aether_one_file(file, vars):
    """ Read in list of variables from a single netCDF file

    Parameters
    ----------
    file : str
        Name of netCDF file to read
    vars : list
        List of desired variable names to read

    Returns
    -------
    data : dict
        Dict with keys 'time', which contains a datetime object specifying the time of
        the file and indices ranging from 0 to `len(vars) - 1`, corresponding to the
        variable names in `vars` that holds arrays of the specified dat
    
    """

    with Dataset(file, 'r') as ncfile:
        # Save the data as numpy arrays, using variable index as a key
        data = {i_var: np.array(ncfile.variables[var]) for i_var, var in enumerate(vars)}
        
        # Calculate the date and time for this data
        base = dt.datetime(1965, 1, 1)
        time = np.array(ncfile.variables['Time'])
        data['time'] = base + dt.timedelta(seconds=time[0])

    return data
    
