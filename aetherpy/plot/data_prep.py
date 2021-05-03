#!/usr/bin/env python
# Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
# Full license can be found in License.md
""" Utilities for slicing and preparing data for plotting
"""

import numpy as np

from aetherpy import logger


def get_cut_index(lons, lats, alts, cut_val, cut_coord='alt'):
    """ Get the cut index needed to obtain a slice in the remaining two coords

    Parameters
    ----------
    lons : array-like
        1D array of longitudes in degrees
    lats : array-like
        1D array of latitudes in degrees
    alts : array-like
        1D array of altitudes in km
    cut_val : int or float
        Data value or grid number (alt only) along which icut will be set
    cut_coord : str
        Expects one of 'lat', 'lon', or 'alt' and will return an index for
        that data, allowing a 2D slice to be created along other two
        coordinates (default='alt')

    Returns
    -------
    icut : int
        Cut index
    cut_data : tuple
        Tuple to select data of the expected dimensions along `icut`
    x_coord : array-like
        Array of data to include along the x-axis
    y_coord : array-like
        Array of data to include along the y-axis
    z_val : float
        Data value for cut index

    Notes
    -----
    `lons`, `lats`, and `alts` do not need to be the same shape

    """
    # Ensure inputs are array-like
    lons = np.asarray(lons)
    lats = np.asarray(lats)
    alts = np.asarray(alts)

    # Set the x and y coordinates
    x_coord = lons if cut_coord in ['alt', 'lat'] else lats
    y_coord = alts if cut_coord in ['lat', 'lon'] else lats

    if cut_coord == 'alt':
        if cut_val < 50:
            # The model doesn't go this low, it must be a grid value
            if len(alts) <= 1 and cut_val != 0:
                logger.warning(''.join(['Requested altitude slice ',
                                        '{:d}'.format(cut_val),
                                        ' is absent, using 0']))
                icut = 0
            else:
                icut = cut_val
        else:
            # This is a data value, find the index closest to the data value
            icut = abs(alts - cut_val).argmin()
            
            # Set a false limit near the top, 3 cells down
            if icut > alts.shape[0] - 3:
                logger.warning(''.join(['Requested altitude slice is above ',
                                        'the recommended upper limit, setting',
                                        ' to recommended upper limit']))
                icut = alts.shape[0] - 3

        z_val = alts[icut]
        cut_data = (slice(0, lons.shape[0], 1), slice(0, lats.shape[0], 1),
                    icut)
    if cut_coord == 'lat':
        # Find closest value for latitude in between sensible limits
        if cut_val < lats[1] or cut_val > lat[-2]:
            logger.warning(''.join(['Requested latitude slice outside the',
                                    ' coordinate range, using midpoint']))
            icut = int(lats.shape[0] / 2)
        else:
            icut = abs(lats - cut_val).argmin()

        z_val = lats[icut]
        cut_data = (slice(0, lons.shape[0], 1), icut,
                    slice(0, alts.shape[0], 1))
    if cut_coord == 'lon':
        # Find closest value for longitude in between sensible limits
        if cut_val < lons[1] or cut_val > lons[-2]:
            logger.warning(''.join(['Requested longitude slice outside the',
                                    ' coordinate range, using midpoint']))
            icut = int(lons.shape[0] / 2)
        else:
            icut = abs(lons - cut_val).argmin()

        z_val = lons[icut]
        cut_data = (icut, slice(0, lats.shape[0], 1),
                    slice(0, alts.shape[0], 1))

    return icut, cut_data, x_coord, y_coord, z_val


def calc_tec(alt, ne, ialt_min=2, ialt_max=-4):
    """Calculate the TEC for the specified altitude range from electron density

    Parameters
    ----------
    alt : array-like
        1D Altitude data in km
    ne : array-like
        Electron density in cubic meters, with altitude as the last axis
    ialt_min : int
        Lowest altitude index to use (default=2)
    ialt_max : int
        Highest altitude index to use, may be negative (default=-4)
    
    Returns
    -------
    tec : array-like
        Array of TEC values in TECU

    Notes
    -----
    TEC = Total Electron Content
    TECU = 10^16 m^-2

    """
    alt = np.asarray(alt)
    ne = np.asarray(ne)

    # Initialize the TEC
    tec = np.zeros(shape=ne[..., 0].shape())


    # Get the range of altitudes to cycle over
    if ialt_max < 0:
        ialt_max += alt.shape[0]

    if ialt_max <= ialt_min:
        raise ValueError('`ialt_max` must be greater than `ialt_min`')

    # Cycle through each altitude bin, summing the contribution from the
    # electron density
    for ialt in np.arange(ialt_min, ialt_max, 1):
        tec += 1000.0 * ne[..., ialt] * (alts[ialt + 1] - alts[ialt - 1]) / 2.0

    # Convert TEC from per squared meters to TECU
    tec /= 1.0e16

    return tec
