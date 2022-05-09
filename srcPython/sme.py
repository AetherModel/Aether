#!/usr/bin/env python                                                           


from supermag_api import *
import argparse
from datetime import datetime
from datetime import timedelta


# ----------------------------------------------------------------------        
# Function to parse input arguments                                             
# ----------------------------------------------------------------------        

def parse_args():

    parser = argparse.ArgumentParser(description = 'Download SME data')

    parser.add_argument('start', metavar = 'start', \
                        help = 'start date as YYYYMMDD[.HHMM]')
    parser.add_argument('end', metavar = 'end', \
			help = 'end date as YYYYMMDD[.HHMM]')

    args = parser.parse_args()

    return args
    
# ----------------------------------------------------------------------        
#                                                                               
# ----------------------------------------------------------------------        

def download_sme_data(start, end):

    startStr = start.strftime('%Y-%m-%dT%M:%D')
    length = (end - start).total_seconds()
    print('--> Downloading AE data')
    userid = 'ridley'
    (status,idxdata) = SuperMAGGetIndices(userid, start, length, 'all')

    nPts = len(idxdata['tval'])
    print('  --> Found %d points' % nPts)

    if (nPts > 0):

        base = datetime(1970,1,1,0,0,0)
        current = base + timedelta(seconds = idxdata['tval'][0])

        ymd = current.strftime('%Y%m%d')
        fileout = 'ae_' + ymd + '.txt'
        print('  --> Writing file ' + fileout)

        l0 = 'File created by python code using SuperMAGGetIndices\n'
        l1 = '============================================================\n'
        l2 = '<year>  <month>  <day>  <hour>  <min>  <sec>  '
        l2 = l2 + '<SME (nT)>  <SML (nT)>  <SMU (nT)>\n'

        fp = open(fileout, 'wb')
        fp.write(l0.encode())
        fp.write("\n".encode())
        fp.write(l1.encode())
        fp.write(l2.encode())

        for i, t in enumerate(idxdata['tval']):
            current = base + timedelta(seconds = t)
            ae = idxdata['SME'][i]
            al = idxdata['SML'][i]
            au = idxdata['SMU'][i]
            out = " %8.2f %8.2f %8.2f" % (ae, al, au)
            ymdhms = current.strftime('%Y  %m  %d  %H  %M  %S')
            line = ymdhms + out + "\n"
            fp.write(line.encode())

        fp.close()

    else:

        fileout = 'none.txt'
        print('  --> ERROR!! No AE data downloaded, no file written!!')

    return fileout


def time_string_to_datetime(inTime):

    yr = inTime[0:4]
    mo = inTime[4:6]
    da = inTime[6:8]

    if (len(inTime) >= 11):
        hr = inTime[9:11]
        if (len(inTime) >= 13):
            mi = inTime[11:13]
        else:
            mi = '00'
    else:
        hr = '00'
        mi = '00'

    outTime = datetime(int(yr), int(mo), int(da), int(hr), int(mi), 0)
    
    return outTime

#------------------------------------------------------------------------------ 
# main code:                                                                    
#------------------------------------------------------------------------------ 

args = parse_args()

# ----------------------------------------                                      
# Set Times:                                                                    

start = time_string_to_datetime(args.start)
end = time_string_to_datetime(args.end)

download_sme_data(start, end)
