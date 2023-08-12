#!/usr/bin/env python
# Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
# Full license can be found in License.md
"""Routines to read Aether files."""

# TODO(#19): Most parts of this file will be changed when switching to xarray

import datetime as dt
from netCDF4 import Dataset
import numpy as np
import os
from struct import unpack
import re

from aetherpy.utils.time_conversion import epoch_to_datetime
from aetherpy import logger


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


def parse_line_into_int_and_string(line, parse_string=True):
    """Parse a data string into integer and string components.

    Parameters
    ----------
    line : str
        Data line with a known format, described in Notes section
    parse_string : bool
        Construct a single white space separated string consisting of the
        remaining portion of `line`

    Returns
    -------
    line_num : int
        Integer corresponding to the first number in `line`
    line_str : str
        Whitespace separated string containing the rest of the `line` data if
        `parse_string` is True, otherwise it returns the original line

    Notes
    -----
    Splits data line using a single space.  The first element is
    expected to be an integer. If desired, the remaining data is returned as a
    single space-separated string.

    """

    split_line = line.split(" ")
    line_num = int(split_line[0])

    if parse_string:
        line_str = " ".join(split_line[1:])
    else:
        line_str = str(line)

    return line_num, line_str


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
    for filename in hfiles:
        if filetype.lower() == "netcdf":
            fheader = read_aether_netcdf_header(filename)
        else:
            fheader = read_aether_ascii_header(filename)

        for hkey in fheader.keys():
            if hkey in header.keys():
                if hkey == 'time':
                    if fheader[hkey] not in header[hkey]:
                        header[hkey].append(fheader[hkey])
                elif hkey == 'vars' and np.any(header[hkey] != fheader[hkey]):
                    header[hkey] = list(np.unique(header[hkey] + fheader[hkey]))
                elif hkey != 'filename' and header[hkey] != fheader[hkey]:
                    raise ValueError(''.join(['header values changed between',
                                              ' files: {:} != {:}'.format(
                                                  header[hkey],
                                                  fheader[hkey])]))
            else:
                header[hkey] = fheader[hkey]

    return header


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

    header = {'filename': filename}

    # Open the file and read the header data
    with Dataset(filename, 'r') as ncfile:
        ncvars = list()
        for var in ncfile.variables.values():
            if len(var.shape) == 3:
                nlons = var.shape[0]
                nlats = var.shape[1]
                nalts = var.shape[2]
                ncvars.append(var.name)

                # Test the dimensions
                if np.any([dim_var not in header.keys()
                           for dim_var in ['nlons', 'nlats', 'nalts']]):
                    header["nlons"] = nlons
                    header["nlats"] = nlats
                    header["nalts"] = nalts
                elif(header['nlons'] != nlons or header['nlats'] != nlats
                     or header['nalts'] != nalts):
                    raise IOError(''.join(['unexpected dimensions for ',
                                           'variable ', var.name, ' in file ',
                                           filename]))

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


def read_aether_ascii_header(filename):
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
        version - version number of the file
        filename - filename of file containing header data

    See Also
    --------
    read_aether_header

    """

    headerfile = filename.replace(".bin", ".header")

    # Test the filename and add it to the header
    if not os.path.isfile(filename) and os.path.isfile(headerfile):
        raise IOError('unknown aether header file: {:}'.format(headerfile))

    header = {"filename": filename}

    with open(headerfile, 'r') as fpin:
        for line in fpin:
            if re.match(r'BLOCKS', line):
                header["nblockslons"] = int(fpin.readline())
                header["nblockslats"] = int(fpin.readline())
                header["nblocksalts"] = int(fpin.readline())

            if re.match(r'NUMERICAL', line):
                header["nvars"], _ = parse_line_into_int_and_string(
                    fpin.readline(), False)
                header["nlons"], _ = parse_line_into_int_and_string(
                    fpin.readline(), False)
                header["nlats"], _ = parse_line_into_int_and_string(
                    fpin.readline(), False)
                header["nalts"], _ = parse_line_into_int_and_string(
                    fpin.readline(), False)

            if re.match(r'VERSION', line):
                header["version"] = float(fpin.readline())

            if re.match(r'TIME', line):
                year = int(fpin.readline())
                month = int(fpin.readline())
                day = int(fpin.readline())
                hour = int(fpin.readline())
                minute = int(fpin.readline())
                second = int(fpin.readline())
                msec = int(fpin.readline())
                header["time"] = dt.datetime(year, month, day, hour,
                                             minute, second, msec)

            mvar = re.match(r'VARIABLE', line)
            if mvar and len(header["vars"]) < 1:
                for i in range(header["nvars"]):
                    nlin, sline = parse_line_into_int_and_string(
                        fpin.readline(), True)
                    header["vars"] = sline.strip()

    return header


def read_blocked_netcdf_header(filename):
    """Read header information from a blocked Aether netcdf file.

    Parameters
    ----------
    filename : str
        An Aether netCDF filename

    Returns
    -------
    header : dict
        A dictionary containing header information from the netCDF file,
        including:
        filename - filename of file containing header data
        nlons - number of longitude grids per block
        nlats - number of latitude grids per block
        nalts - number of altitude grids per block
        nblocks - number of blocks in file
        vars - list of data variable names
        time - datetime for time of file

    Raises
    --------
    IOError
        If the input file does not exist
    KeyError
        If any expected dimensions of the input netCDF file are not present

    Notes
    -----
    This routine only works with blocked Aether netCDF files.

    """

    # Checks for file existence
    if not os.path.isfile(filename):
        raise IOError(f"unknown aether netCDF blocked file: {filename}")

    header = {'filename': filename}  # Included for compatibility

    with Dataset(filename, 'r') as ncfile:
        # Process header information: nlons, nlats, nalts, nblocks
        header['nlons'] = len(ncfile.dimensions['lon'])
        header['nlats'] = len(ncfile.dimensions['lat'])
        header['nalts'] = len(ncfile.dimensions['z'])
        header['nblocks'] = len(ncfile.dimensions['block'])

        # Included for compatibility ('vars' slices time out for some reason)
        header['vars'] = list(ncfile.variables.keys())[1:]
        header['time'] = epoch_to_datetime(
            np.array(ncfile.variables['time'])[0])

    return header


def read_blocked_netcdf_file(filename, file_vars=None):
    """Read all data from a blocked Aether netcdf file.

    Parameters
    ----------
    filename : str
        An Aether netCDF filename
    file_vars : list or NoneType
        List of desired variable neames to read, or None to read all
        (default=None)

    Returns
    -------
    data : dict
        A dictionary containing all data from the netCDF file, including:
        filename - filename of file containing header data
        nlons - number of longitude grids per block
        nlats - number of latitude grids per block
        nalts - number of altitude grids per block
        nblocks - number of blocks in file
        vars - list of data variable names
        time - datetime for time of file
        The dictionary also contains a read_routines.DataArray keyed to the
        corresponding variable name. Each DataArray carries both the variable's
        data from the netCDF file and the variable's corresponding attributes.

    Raises
    --------
    IOError
        If the input file does not exist
    KeyError
        If any expected dimensions of the input netCDF file are not present

    Notes
    -----
    This routine only works with blocked Aether netCDF files.

    """

    # Checks for file existence
    if not os.path.isfile(filename):
        raise IOError(f"unknown aether netCDF blocked file: {filename}")

    # NOTE: Includes header information for easy access until
    #       updated package structure is confirmed
    # Initialize data dict with defaults (will remove these defaults later)
    data = {'filename': filename,
            'units': '',
            'long_name': None}

    with Dataset(filename, 'r') as ncfile:
        # Process header information: nlons, nlats, nalts, nblocks
        data['nlons'] = len(ncfile.dimensions['lon'])
        data['nlats'] = len(ncfile.dimensions['lat'])
        data['nalts'] = len(ncfile.dimensions['z'])
        data['nblocks'] = len(ncfile.dimensions['block'])

        # Included for compatibility
        data['vars'] = [var for var in ncfile.variables.keys()
                        if file_vars is None or var in file_vars]

        # Fetch requested variable data
        for key in data['vars']:
            var = ncfile.variables[key]  # key is var name
            data[key] = DataArray(np.array(var), var.__dict__)

        data['time'] = epoch_to_datetime(np.array(ncfile.variables['time'])[0])

    return data


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
    logger.info("Reading file : ", file_to_read)

    data = {hkey: header[hkey] for hkey in ["version", "nlons", "nlats",
                                            "nalts" "nvars", "vars"]}
    data["time"] = header["time"][ifile]

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


def read_gitm_headers(filelist, finds=-1):
    """Read ancillary information from a GITM file.

    Parameters
    ----------
    filelist : array-like
        Array-like object of names for netCDF files
    finds : int, list, or slice
        Index(es) for file(s) from which the header information will be read.
        (default=-1)

    Returns
    -------
    header : dict
        A dictionary containing header information from the netCDF files,
        including:
        nlons - number of longitude grids
        nlats - number of latitude grids
        nalts - number of altitude grids
        vars - list of data variable names
        time - list of datetimes for the processed file start times
        filename - list of the processed filenames
        version - file version number

    Raises
    ------
    IOError
        If any one of the input files encounters an unexpected value or
        dimension.

    Notes
    -----
    This routine obtains the same info as `read_aether_headers`

    """

    # Initialize the output
    header = {"vars": [], "time": [], "filename": np.asarray(filelist)}

    # Ensure the filelist is array-like, allowing slicing of input
    if header['filename'].shape == ():
        header['filename'] = np.asarray([filelist])
    else:
        header['filename'] = list(header['filename'])

    # Ensure selected files are list-like
    hfiles = np.asarray(header['filename'][finds])
    if hfiles.shape == ():
        hfiles = np.asarray([hfiles])

    for filename in hfiles:
        # Read in the header from the binary file
        file_vars = list()

        with open(filename, 'rb') as fin:
            # Test to see if the correct endian is being used
            end_char = '>'
            raw_rec_len = fin.read(4)
            rec_len = (unpack(end_char + 'l', raw_rec_len))[0]
            if rec_len > 10000 or rec_len < 0:
                # Ridiculous record length implies wrong endian, fix here
                end_char = '<'
                rec_len = (unpack(end_char + 'l', raw_rec_len))[0]

            # Read version; read fortran footer+header.
            file_version = unpack(end_char + 'd', fin.read(rec_len))[0]
            _, rec_len = unpack(end_char + '2l', fin.read(8))

            # Test the version number
            if 'version' not in header.keys():
                header['version'] = file_version
            elif header['version'] != file_version:
                raise IOError(''.join(['unexpected version number in file ',
                                       filename]))

            # Read grid size information.
            nlons, nlats, nalts = unpack(end_char + 'lll', fin.read(rec_len))
            _, rec_len = unpack(end_char + '2l', fin.read(8))

            # Test the dimensions
            if np.any([dim_var not in header.keys()
                       for dim_var in ['nlons', 'nlats', 'nalts']]):
                header["nlons"] = nlons
                header["nlats"] = nlats
                header["nalts"] = nalts
            elif(header['nlons'] != nlons or header['nlats'] != nlats
                 or header['nalts'] != nalts):
                raise IOError(''.join(['unexpected dimensions in file ',
                                       filename]))

            # Read number of variables.
            num_vars = unpack(end_char + 'l', fin.read(rec_len))[0]
            _, rec_len = unpack(end_char + '2l', fin.read(8))

            # Collect variable names.
            for ivar in range(num_vars):
                vcode = unpack(end_char + '%is' % (rec_len),
                               fin.read(rec_len))[0]
                var = vcode.decode('utf-8').replace(" ", "")
                file_vars.append(var)
                _, rec_len = unpack(end_char + '2l', fin.read(8))

            # Test the variable names
            if len(header["vars"]) == 0:
                header["vars"] = list(file_vars)
            elif header["vars"] != file_vars:
                raise IOError(''.join(['unexpected number or name of ',
                                       'variables in file ', filename]))

            # Extract time
            out_time = np.array(unpack(end_char + 'lllllll', fin.read(rec_len)))
            out_time[-1] *= 1000  # Convert from millisec to microsec
            header["time"].append(dt.datetime(*out_time))

    return header


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
    logger.info("Reading file : ", filename)

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


def read_gitm_file(filename, file_vars=None):
    """Read list of variables from one GITM file.

    Parameters
    ----------
    filename : str
        GITM file to read
    file_vars : list or NoneType
        List of desired variable names to read or None to read all
        (default=None)

    Returns
    -------
    data : dict
        Dict with keys 'time', which contains a datetime object specifying the
        time of the file and zero-offset indices, corresponding to the
        variable names in `file_vars` that holds arrays of the specified data.
        Also contains version number, dimensions, and a list of the variable
        names.

    """

    data = {"vars": []}
    logger.info("Reading file : ", filename)

    if not os.path.isfile(filename):
        raise IOError('input file does not exist')

    with open(filename, 'rb') as fin:
        # Determine the correct endian
        end_char = '>'
        raw_rec_len = fin.read(4)
        rec_len = (unpack(end_char + 'l', raw_rec_len))[0]
        if rec_len > 10000 or rec_len < 0:
            # Ridiculous record length implies wrong endian.
            end_char = '<'
            rec_len = (unpack(end_char + 'l', raw_rec_len))[0]

        # Read version; read fortran footer+data.
        data["version"] = unpack(end_char + 'd', fin.read(rec_len))[0]

        _, rec_len = unpack(end_char + '2l', fin.read(8))

        # Read grid size information.
        data["nlons"], data["nlats"], data["nalts"] = unpack(
            end_char + 'lll', fin.read(rec_len))
        _, rec_len = unpack(end_char + '2l', fin.read(8))

        # Read number of variables.
        num_vars = unpack(end_char + 'l', fin.read(rec_len))[0]
        _, rec_len = unpack(end_char + '2l', fin.read(8))

        if file_vars is None:
            file_vars = np.arange(0, num_vars, 1)

        # Collect variable names in a list
        for ivar in range(num_vars):
            vcode = unpack(end_char + '%is' % (rec_len),
                           fin.read(rec_len))[0]
            var = vcode.decode('utf-8').replace(" ", "")
            data['vars'].append(var)
            dummy, rec_lec = unpack(end_char + '2l', fin.read(8))

        # Extract time
        rec_time = np.array(unpack(end_char + 'lllllll', fin.read(28)))
        rec_time[-1] *= 1000  # convert from millisec to microsec
        data["time"] = dt.datetime(*rec_time)

        # Header is this length:
        # Version + start/stop byte
        # nlons, nlats, nalts + start/stop byte
        # num_vars + start/stop byte
        # variable names + start/stop byte
        # time + start/stop byte

        iheader_length = 84 + num_vars * 48

        ntotal = data["nlons"] * data["nlats"] * data["nalts"]
        idata_length = ntotal * 8 + 8

        # Save the data for the desired variables
        for ivar in file_vars:
            fin.seek(iheader_length + ivar * idata_length)
            sdata = unpack(end_char + 'l', fin.read(4))[0]
            data[ivar] = np.array(
                unpack(end_char + '%id' % (ntotal), fin.read(sdata))).reshape(
                    (data["nlons"], data["nlats"], data["nalts"]), order="F")

    return data
