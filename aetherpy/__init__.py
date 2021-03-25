#!/usr/bin/env python
# Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
# Full license can be found in License.md

import logging

from aetherpy import io, utils, plot

# Define global variables
__version__ = '0.0.1a'

# Define a logger object to allow easier log handling
logging.raiseExceptions = False
logger = logging.getLogger('aetherpy_logger')
