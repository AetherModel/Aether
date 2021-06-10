#!/usr/bin/env python

import sys
import datetime as dt
import numpy as np
from omniweb import *

def write_swmf_imf_file(data, fileout):

    fp = open(fileout, 'wb')

    fp.write("\n".encode())
    fp.write("Data downloaded from OMNIWeb and processed by omniweb_read.py\n".encode())
    fp.write("\n".encode())
    fp.write("No AlphaRatio was considered\n".encode())
    fp.write("\n".encode())
    fp.write("#TIMEDELAY\n".encode())
    fp.write("0.0\n".encode())
    fp.write("\n".encode())
    fp.write("#START\n".encode())
    
    i = 0
    while (np.abs(data["vx"][i]) > 10000.0):
        i += 1
    
    vx0 = data["vx"][i]
    vy0 = data["vy"][i]
    vz0 = data["vz"][i]
    n0 = data["n"][i]
    temp0 = data["t"][i]
    
    for i, t in enumerate(data["times"]):
        bx = data["bx"][i]
        by = data["by"][i]
        bz = data["bz"][i]
        vx = data["vx"][i]
        vy = data["vy"][i]
        vz = data["vz"][i]
        n = data["n"][i]
        temp = data["t"][i]

        if (np.abs(vx) > 10000):
            vx = vx0
            vy = vy0
            vz = vz0
            n = n0
            temp = temp0
        else:
            vx0 = vx
            vy0 = vy
            vz0 = vz
            n0 = n
            temp0 = temp
            
        if (np.abs(bx) < 1000):
            sTime = t.strftime(' %Y %m %d %H %M %S 000')
            sImf = "%8.2f %8.2f %8.2f" % (bx, by, bz)
            sSWV = "%9.2f %9.2f %9.2f" % (vx, vy, vz)
            sSWnt = "%8.2f %11.1f" % (n, temp)

        line = sTime + sImf + sSWV + sSWnt + "\n"

        fp.write(line.encode())

    fp.close()
        
#--------------------------------------------------------------------------------
#SCRIPT USE
#example command line input: python omniweb_read.py 20110620 20110623 -all

args = sys.argv
#assuming first two args are the start/end dates, then info desired

results = download_omni_data(args[1], args[2], "-all")
data = parse_omni_data(results)


print(data["Vars"])

fileout = data["times"][0].strftime('imf%Y%m%d.dat')
write_swmf_imf_file(data, fileout)
