#!/usr/bin/env python
# Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
# Full license can be found in License.md
""" Utilities for outputing movies
"""

from glob import glob
import os


def setup_movie_dir(movie_dir, overwrite=True):
    """Set up a directory for movie files

    Parameters
    ----------
    movie_dir : str
        Output filename with directory, but no extention
    overwrite : bool
        Overwrite an existing movie of the same name (default=True)

    Returns
    -------
    img_names : str
        Image name formatting string

    """
    
    # Test the output directory for existence and existing image files
    if os.path.isdir(movie_dir):
        oldfiles = glob(os.path.join(movie_dir, "image_????.png"))
        if len(oldfiles) > 0:
            if overwrite:
                for ofile in oldfiles:
                    os.remove(ofile)
            else:
                raise IOError('files present in movie directory: {:}'.format(
                    movie_dir))
    else:
        os.makedirs(movie_dir)

    # Create the movie image naming string
    img_names = os.path.join(movie_dir, "image_{:04d}.png")

    return img_names


def save_movie(fileroot, ext='mp4', rate=30, overwrite=True):
    """Save the output as a movie

    Parameters
    ----------
    fileroot : str
        Output filename with directory, but no extention
    ext : str
        Movie extention (default='mp4')
    rate : int
        Movie frame rate (default=30)
    overwrite : bool
        Overwrite an existing movie of the same name (default=True)

    Notes
    -----
    Uses ffmpeg to create the movie, this must be installed for success

    """
    # Construct the output filenames
    outfile = ".".join(fileroot, ext)
    image_files = os.path.join(fileroot, 'image_%04d.png')
    
    # Test the output file
    if os.path.isfile(outfile):
        if overwrite:
            os.remove(outfile)
        else:
            raise IOError('movie file {:} already exists'.format(outfile))

    # Construct the movie commannd
    command = "ffmpeg -r {:d} -i {:s} {:s}".format(rate, image_files, outfile)
    os.system(command)
    return
