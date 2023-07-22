#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

def sat_line(line):
    """Convert a line of satellite file into time and position."""
    if line is None:
        return None
    line = line.strip().split(", ")
    return [int(line[i]) for i in range(6)] + [float(line[i]) for i in range(6, 9, 1)]

def log_line(line):
    """Convert a line of log file into time and position."""
    if line is None:
        return None
    line = line.strip().split(' ')
    # Need to neglect milliseconds here
    return [int(line[i]) for i in range(6)] + [float(line[i]) for i in range(7, 10, 1)]

def calc_diff(x:list, y:list):
    """Calculate the error between two positions."""
    if abs(x[0] - y[0]) > 180:
        # One point close to 0 and the other close to 360
        return 360 - abs(x[0] - y[0]) + abs(x[1] - y[1]) + abs(x[2] - y[2])
    else:
        return sum([abs(a - b) for a, b in zip(x, y)])


def main():
    """Check whether the output is correct by visualizing the difference of location."""
    # Open the satellite file and log file
    with open("UA/inputs/sat_20110320.csv", 'r', encoding="utf-8") as sat:
        with open("UA/output/sat_20110320_log.txt", 'r', encoding="utf-8") as log:
            # Skip the first two lines of satellite file and the first line of log file
            next(sat)
            next(sat)
            next(log)
            # Initialize the times and percentage of error to plot
            times_plot = []
            diff_plot = []
            sv = []
            lv = []
            # Use the iterator counter as time
            iter_count = 0
            # Start the two-pointer approach
            sat_val = sat_line(next(sat, None))
            log_val = log_line(next(log, None))
            while sat_val is not None and log_val is not None:
                sat_time = sat_val[:6]
                log_time = log_val[:6]
                if sat_time == log_time:
                    # Add plot list, counter, and go to next for both
                    times_plot.append(iter_count)
                    sv.append(sat_val[6:])
                    lv.append(log_val[6:])
                    diff_plot.append(calc_diff(sat_val[6:], log_val[6:]))
                    iter_count += 1
                    sat_val = sat_line(next(sat, None))
                    log_val = log_line(next(log, None))
                elif sat_time < log_time:
                    # Sat go to next
                    sat_val = sat_line(next(sat, None))
                else:
                    # Log go to next, add counter
                    iter_count += 1
                    log_val = log_line(next(log, None))
            # Plot error vs time

            fig = plt.figure(figsize = (10,10))
            ax1 = fig.add_subplot(211)
            ax2 = fig.add_subplot(212)

            sv = np.array(sv)
            lv = np.array(lv)
            ax1.scatter(times_plot, sv[:,0])
            ax1.scatter(times_plot, lv[:,0])
            ax1.set_xlabel("Time (15s)")
            ax1.set_ylabel("Longitudes (deg)")
            ax2.plot(times_plot, diff_plot)
            ax1.set_xlabel("Time (15s)")
            ax2.set_ylabel("Total Difference (deg + km)")
            
            fig.savefig("Satellite_log.png")
            plt.close()


if __name__ == "__main__":
    main()
