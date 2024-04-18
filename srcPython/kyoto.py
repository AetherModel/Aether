#!/usr/bin/env python

import requests
import re
import datetime as dt
import numpy as np

# This data is downloaded from:
# http://wdc.kugi.kyoto-u.ac.jp/dst_final/YYYYMM/index.html

def download_kyoto_dst(date):

    # assume date is YYYYMMDD, so need the substring YYYYMM:
    yyyymm = date[0:6]
    url = "http://wdc.kugi.kyoto-u.ac.jp/dst_final/"+yyyymm+"/index.html"

    res = requests.get(url = url)
    return res.text     #this returns the info file from the get request

def parse_kyoto_dst(filestring, date):

    year = int(date[0:4])
    month = int(date[4:6])

    lines = filestring.splitlines()

    IsFound = 0
    iLine = 0
    while (not IsFound):
        line = lines[iLine]
        m = re.match(r'DAY.*',line)
        if m:
            IsFound = 1
        iLine += 1

    data = {}
    data["Vars"] = ["dst"]
    data["nVars"] = 1
    data["times"] = []
    data["dst"] = []
        
    IsFound = 0
    while (not IsFound):
        line = lines[iLine]
        m = re.match(r'<!-- vvvvv S yyyymm_part3.html vvvvv -->.*',line)
        if m:
            IsFound = 1
        else:
            if (len(line) > 40):
                day = int(line[0:2])
                iStart = 3
                for hour in np.arange(0,24):
                    iEnd = iStart + 4
                    dst = line[iStart:iEnd]
                    data["dst"].append(int(dst))
                    time = dt.datetime(year, month, day, hour, 0, 0)
                    data["times"].append(time)
                    iStart = iEnd
                    if ((hour == 7) or (hour == 15)):
                        iStart += 1

        iLine += 1
        
    return data

    
    

