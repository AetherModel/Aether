#!/usr/bin/env python3

import sys
import os
import contextlib
import heapq

def cmp_time(line:str):
    """Return a list of int of size 7 representing the time of the input stirng."""
    time = line.split()
    return [int(time[i]) for i in range(7)]

def main():
    """Combine the temporary log files from different grids into one file."""
    try:
        with open("SAT_COMBINE_README.txt", 'r', encoding="utf-8") as config:
            # The first line is the header
            header = config.readline()
            # The second line is the number of grids
            nGrid = int(config.readline())
            # The third line is the final log file names
            log_files = config.readline().strip().split(' ')

            # Iterate on each satellite
            for log_file in log_files:
                # All files related to this satellite
                filenames = [f"{log_file}.g{iGrid:04d}" for iGrid in range(nGrid)]
                # Open file stream
                with contextlib.ExitStack() as stack:
                    try:
                        files = [stack.enter_context(open(filename, 'r', encoding="utf-8"))
                                    for filename in filenames]
                    except FileNotFoundError:
                        print("Fail to find one or more of the required files:")
                        print(filenames)
                        sys.exit(1)
                    # The input from each file is sorted, so use heapq to merge
                    with open(log_file, 'w', encoding="utf-8") as out:
                        # Write the header
                        out.write(header)
                        # Write each line in the order of time
                        for line in heapq.merge(*files, key=cmp_time):
                            out.write(line)
                # Delete the temporary files related to this satellite
                for filename in filenames:
                    os.remove(filename)
        # Delete SAT_COMBINE_README.txt
        os.remove("SAT_COMBINE_README.txt")
    except FileNotFoundError:
        print("Can not find file: SAT_COMBINE_README.txt in directory", os.getcwd())
        print("Use `pwd` to check which directory you are running this script in.")

if __name__ == "__main__":
    main()
