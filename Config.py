#!/usr/bin/env python

import argparse
from glob import glob
import os

# ----------------------------------------------------------------------
# Function to parse input arguments
# ----------------------------------------------------------------------

def parse_args(oss):

    parser = argparse.ArgumentParser(description = 'Aether Configure Script')

    parser.add_argument('-netcdf', \
                        help='Set Aether to use netCDF output', \
                        action='store_true')

    parser.add_argument('--verbose', '-v', \
                        help='Turn on Verbose', \
                        action='store_true')
    
    parser.add_argument('-os', \
                        help='Set OS to use (copies correct Makefile.OS)', \
                        choices = oss, \
                        default='none')

    args = parser.parse_args()
    
    return args



# ----------------------------------------------------------------------
# Execute a command
# ----------------------------------------------------------------------

def execute_command(command, verbose):

    if (verbose):
        print(command)
    os.system(command)

# ----------------------------------------------------------------------
# Main Code
# ----------------------------------------------------------------------

makefiles = glob('src/Makefile.NetCDF.*')
oss = []
for file in makefiles:
    oper = file[20:]
    if (oper != 'OS'):
        oss.append(oper)

args = parse_args(oss)

output_dir = "src"
output_target = "output.cpp"
if (args.netcdf):
    output_source = "output_netcdf.cpp"
else:
    output_source = "output_binary.cpp"

command = "cd " + output_dir + " ; /bin/rm -f " + output_target + " ; cd -"
execute_command(command, args.verbose)

command = "cd " + output_dir + \
    " ; ln -s " + output_source + " " + output_target + " ; cd -"
execute_command(command, args.verbose)

# ------------------------------------------------
# Makefile
# ------------------------------------------------

output_target = "Makefile.NetCDF"
if (args.netcdf):
    output_source = "Makefile.NetCDF." + args.os
else:
    output_source = "Makefile.NoNetCDF"

command = "cd " + output_dir + " ; /bin/rm -f " + output_target + " ; cd -"
execute_command(command, args.verbose)

command = "cd " + output_dir + \
    " ; ln -s " + output_source + " " + output_target + " ; cd -"
execute_command(command, args.verbose)

