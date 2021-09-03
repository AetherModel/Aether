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

    parser.add_argument('-rundir', \
                        help='Make a run directory', \
                        nargs = '?', \
                        const='run')

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

makefiles = glob('src/Makefile.*')
oss = []
for file in makefiles:
    oper = file[13:]
    if (oper != 'OS'):
        oss.append(oper)

args = parse_args(oss)

print(args.netcdf)
print(args.os)
print(args.verbose)
print(args.rundir)

if (args.rundir):
    dir = args.rundir
    isDir = os.path.isdir(dir)
    if (isDir):
        print("Directory " + dir + " exists!")
        print("Please do something about this!")
        exit()
    else:
        # Let's make the run directory and put files in the right place!

        aetherExe = os.getcwd() + "/build/aether"

        print(aetherExe)
        
        command = "mkdir "+dir
        execute_command(command, args.verbose)

        command = "mkdir " + dir + "/UA"
        execute_command(command, args.verbose)

        command = "mkdir " + dir + "/UA/inputs"
        execute_command(command, args.verbose)

        command = "mkdir " + dir + "/UA/output"
        execute_command(command, args.verbose)

        command = "mkdir " + dir + "/UA/restartOut"
        execute_command(command, args.verbose)

        command = "cd " + dir + "/UA; ln -s restartOut restartIn ; cd -"
        execute_command(command, args.verbose)

        command = "cp inputs/aether.in " + dir
        execute_command(command, args.verbose)

        command = "cp inputs/UA/*.csv " + dir + "/UA/inputs"
        execute_command(command, args.verbose)

        command = "cp inputs/UA/*.txt " + dir + "/UA/inputs"
        execute_command(command, args.verbose)

        command = "cp inputs/*.csv " + dir + "/UA/inputs"
        execute_command(command, args.verbose)

        command = "cp inputs/*.txt " + dir + "/UA/inputs"
        execute_command(command, args.verbose)

        command = "cd " + dir + "/UA; ln -s restartOut restartIn ; cd -"
        execute_command(command, args.verbose)

#command = "cd " + output_dir + " ; /bin/rm -f " + output_target + " ; cd -"
#execute_command(command, args.verbose)


