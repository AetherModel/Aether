#!/usr/bin/env python
# Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
# Full license can be found in License.md

# Define a logger object to allow easier log handling
import logging

logging.raiseExceptions = False
logger = logging.getLogger('aetherpy_logger')

# Import the sub-modules
from aetherpy import io, utils, plot

# Define global variables
__version__ = '0.0.1a'
