#!/usr/bin/env python
# Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
# Full license can be found in License.md
""" Utilities for handling command line inputs
"""

import sys


def bool_string(line):
    """ Determine whether a string should be True or False

    Parameters
    ----------
    line : string
        Line to be tested

    Returns
    -------
    bout : bool
        Boolean output (True/False)

    Raises
    ------
    ValueError
        If the value cannot be interpreted as True or False

    Notes
    -----
    Accepts empty, true, false, t, f, 1, and 0 in any capitalization combo

    """

    line = line.lower()

    bout = None
    if line in ['true', 't', '1', '']:
        bout = True
    elif line in ['false', 'f', '0']:
        bout = False

    if bout is None:
        raise ValueError('input not interpretable as a boolean')

    return bout


def none_string(line):
    """ Determine whether a string should be None

    Parameters
    ----------
    line : string
        Line to be tested

    Returns
    -------
    out_line : string or NoneType
        None if all-lowercase version of line is "none" or line is zero length.
        Otherwise returns original value of line

    """
    out_line = None if line.lower() == "none" or len(line) == 0 else line
    return out_line


def process_command_line_input():
    """ Process command line input, needed to possible ipython use

    Returns
    -------
    input_args : list
        List of input arguements

    """

    input_args = sys.argv
    if input_args[0].find('ipython') >= 0:
        input_args = list()
    else:
        input_args.pop(0)

    return input_args
