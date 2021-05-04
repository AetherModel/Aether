#!/usr/bin/env python
# Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
# Full license can be found in License.md
"""Routines to find and retrieve files
"""

from glob import glob
import os

from aetherpy import logger


def get_gitm_filelist(file_dir):
    """Get a list of GITM files from a specified directory

    Parameters
    ----------
    file_dir : str
        File directory to search or a glob search string

    Returns
    -------
    filelist : list
        List of GITM filenames

    Notes
    -----
    Only retrieves 1D files if 3D files are not present

    """
    # Check to see if the directory exists
    if os.path.isdir(file_dir):
        # Get a list of 3DALL files or 1DALL files
        filelist = glob(os.path.join(file_dir, '3DALL*.bin'))

        if len(filelist) == 0:
            logger.info("".join(["No 3DALL files found in ", file_dir,
                                 ", checking for 1DALL."]))
            filelist = glob(os.path.join(file_dir, '1DALL*.bin'))

            if len(filelist) == 0:
                logger.warning('No 1DALL files found in {:s}'.format(file_dir))
    else:
        # This may be a glob search string
        filelist = glob(file_dir)

        if len(filelist) == 0:
            logger.warning('No files found using search string: {:s}'.format(
                file_dir))

    return filelist
