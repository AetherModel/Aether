#!/usr/bin/env python

from glob import glob
import csv
import matplotlib.pyplot as plt
import datetime as dt

filelist = sorted(glob('UA/output/*log*.txt'))

file = filelist[0]

# Initialize lists to hold time and temperature values
time_values = []
temp_values = []

# Read the file and process the data
with open(file, 'r') as file:
    reader = csv.DictReader(file, delimiter=' ')
    start = None
    first_line = True
    for row in reader:
        if first_line:
            start = dt.datetime(int(row['year']),int(row['month']),int(row['day']),int(row['hour']),int(row['minute']),int(row['second']))
            first_line = False
        current = dt.datetime(int(row['year']),int(row['month']),int(row['day']),int(row['hour']),int(row['minute']),int(row['second']))
        duration = current - start

        time_values.append(int(duration.seconds))
        temp_values.append(float(row['max_temp']))

# Create the plot
plt.plot(time_values, temp_values)
plt.xlabel('Time (s)')
plt.ylabel('Temperature')
plt.title('Time vs Temperature')

plotFile = 'log.png'
print('  Outputting file : ', plotFile)
plt.savefig(plotFile)
