#!/usr/bin/env python
# Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
# Full license can be found in License.md
"""Routines to perform temporal calculations
"""

import datetime as dt


def epoch_to_datetime(epoch_time):
    """Convert from epoch seconds to datetime

    Parameters
    ----------
    epoch_time : int
        Seconds since 1 Jan 1965

    Returns
    -------
    dtime : dt.datetime
        Datetime object corresponding to `epoch_time`

    """

    dtime = dt.datetime(1965, 1, 1) + dt.timedelta(seconds=epoch_time)

    return dtime


def datetime_to_epoch(dtime):
    """Convert datetime to epoch seconds

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
