#!/usr/bin/env python3
""" 
Routines to read Armadillo objects
==================================


This is just a couple of functions which will read Armadillo exports into Python.
Directly exporting Armadillo cubes can be very useful during development and/or
debugging, but is not suitable for production runs.

To export an Armadillo cube from Aether:
(substitute the cube name & file name. This does support multiple processors)

grid.geoLon_scgc.save("geoLon_" + tostr(iProc, 3) + ".txt", arma_ascii));

Notes:
- tostr() is defined in src/tools.cpp of Aether & zero-pads an int to return a str.
- This will rewrite the existing file each time it is called.
- Output is to the same directory the executable is called from.
- This uses the arma_ascii format, which is way less efficient than HDF5 or binary. 
  I have found this format is the easiest to work with, but your mileage may vary.
- The several python armadillo implementations look abandoned and/or did not work for me.
- See the armadillo documentation for more information on saving cubes or other data types:
  https://arma.sourceforge.net/docs.html#save_load_mat

"""

import numpy as np
from glob import glob
import os, errno

def check_file_inputs(files):
    """ Make sorted list of files (that exist) from a str or list

    Inputs
    ------
        files (str or list) Can be list of files, single file, directory, or a pattern to glob

    Returns
    -------
        list: sorted list of files that so indeed exist
    
    """
    
    if isinstance(files, str): # Probably need to glob
        if "*" in files: # Definitely need to glob
            files2read = np.sort(glob(files))
        elif os.path.isfile(files):
            files2read = [files] # Single file needs to be made into list
        elif os.path.isdir(files):
            # We were given a directory. Read all .txt files without log in name
            files_ = np.sort(glob(os.path.join(files, "*.txt")))
            files2read = [f for f in files_ if "log" not in f]
        else: # pretty error message from stack overflow
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), *files)
        
        if len(files2read) == 0:
            raise ValueError(
                f"Could not find any armadillo cubes from '{files}'."
                " Check path or provide files.\n")

        return files2read

    # Sort list & check if all files exist. error if not.
    try: # attempt to handle anything listlike (np arrays, dict keys, etc.)
        files2read = [f for f in np.sort(files) if os.path.isfile(f)]
        if len(files2read) != len(files):
            bad_files = [f for f in files if f not in files2read]
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), *bad_files)

    except: # Not expected types
        raise TypeError("Need list or str. Could not handle type: " + type(files))
    
    return files2read


def cube2np(files2read):
    """ Read armadillo cubes from .txt files, automatically globs input.
        return np array of shape (nFiles, n_x, n_y, n_z)

    Inputs
    ------
        files (str or list-like): either path to files or list of files. If it's a str,
            the pattern is globbed & sorted, or the directory's .txt files are sorted.
            If it's list-like, the list is sorted.

    Outputs
    -------
        np.array of shape (nFiles, n_x, n_y, n_z) & dtype float. If we are only reading
            one file, return shape is just (n_x, n_y, n_z)

    Usage
    -----

    lons = cube2np("../run/geolon_*.txt")
    lons = cube2np(np.sort(glob.glob("../run/geolon_*.txt")))

    """

    # Sanitize input
    files2read = check_file_inputs(files2read)

    out = [] # output holder
    for thisf in files2read:
        with open(thisf, 'r') as f:
            _ = f.readline() # first line is a header, not needed
            shape = f.readline().strip() # next line holds the shape of the cube
            shape = shape.split(' ')
            if len(shape) != 3:
                raise ValueError(
                    f"File ({thisf}) does not appear to be an armadillo cube.\n"
                    f"Found shape: {shape}")
            shape = np.array(shape, dtype=int) # convert shape to np array of int's
        # Read in cube, transform shape, use same indexing as arma does
        one_cube = np.loadtxt(thisf, skiprows=2, ).reshape(shape)
        
        out.append(one_cube) # speed not a huge issue, work with lists

    # remove 0th dimension if we only are reading one file
    if len(files2read) == 1:
        out = out[0] 

    return np.array(out)




def mat2np(files2read):
    """ Read armadillo matrices from .txt files, automatically globs input.
        return np array of shape (nFiles, n_x, n_y)

    Inputs
    ------
        files (str or list-like): either path to files or list of files. If it's a str,
            the pattern is globbed & sorted, or the directory's .txt files are sorted.
            If it's list-like, the list is sorted.

    Outputs
    -------
        np.array of shape (nFiles, n_x, n_y) & dtype float. If we are only reading
            one file, return shape is just (n_x, n_y)

    Usage
    -----

    lons = mat2np("../run/geolon_*.txt")
    lons = mat2np(np.sort(glob.glob("../run/geolon_*.txt")))

    """

    # Sanitize input
    files2read = check_file_inputs(files2read)

    out = [] # output holder
    for thisf in files2read:
        with open(thisf, 'r') as f:
            _ = f.readline() # first line is a header, not needed
            shape = f.readline().strip() # next line holds the shape of the cube
            shape = shape.split(' ')
            if len(shape) != 2:
                raise ValueError(
                    f"File ({thisf}) does not appear to be an armadillo matrix.\n"
                    f"Found shape: {shape}")
            shape = np.array(shape, dtype=int) # convert shape to np array of int's
        
        one_mat = np.loadtxt(thisf, skiprows=2).reshape(shape)

        out.append(ls) # speed not a huge issue, work with lists

    # remove 0th dimension if we only are reading one file
    if len(files2read) == 1:
        out = out[0] 

    return np.array(out)
