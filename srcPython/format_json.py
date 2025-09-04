#!/usr/bin/env python3

import os
import argparse
import json

# ----------------------------------------------------------------------------
# Get arguments as inputs into the code
#-----------------------------------------------------------------------------

def get_args():

    parser = argparse.ArgumentParser(
        description = 'reformat json files')
    
    # Get the files to plot:
    parser.add_argument('filelist', nargs='+', \
                        help = 'list of files for formatting')

    args = parser.parse_args()

    return args

# ----------------------------------------------------------------------
# do system command
# ----------------------------------------------------------------------

def run_command(command, verbose = False):
    if (verbose):
        print(" -> Running Command : ")
        print("    ", command)
    os.system(command)
    return True

# Needed to run main script as the default executable from the command line
if __name__ == '__main__':

    # Get the input arguments
    args = get_args()
    filelist = args.filelist

    for file in filelist:

        fileSave = file + '.orig'
        command = 'mv ' + file + ' ' + fileSave
        run_command(command, verbose = True)
        
        with open(fileSave, 'r') as handle:
            print('-> Reading : ', fileSave)
            parsed = json.load(handle)        

            fpOut = open(file, 'w')
            print('-> Writing : ', file)
            json.dump(parsed, fpOut, indent=4)
            fpOut.close()

