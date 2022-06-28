#!/usr/bin/env python

import netCDF4 as nc
from scipy.io import FortranFile
import numpy as np
from datetime import datetime
from struct import unpack

# This should work on most systems:
endChar='<'

def read_record_float(fp, nPoints):
    recLen=unpack(endChar+'l',fp.read(4))[0]
    values = np.array(unpack(endChar+'%if'%(nPoints),fp.read(recLen)))
    endLen = np.array(unpack(endChar+'l',fp.read(4))[0])
    return values
    
def read_record_long(fp, nPoints):
    recLen = unpack(endChar+'l',fp.read(4))[0]
    values = np.array(unpack(endChar+'%il'%(nPoints),fp.read(recLen)))
    endLen = unpack(endChar+'l',fp.read(4))[0]
    return values
    
def read_record_string(fp):
    recLen = unpack(endChar+'l',fp.read(4))[0]
    value = unpack(endChar+'%is'%(recLen),fp.read(recLen))[0]
    endLen = unpack(endChar+'l',fp.read(4))[0]
    return value.decode("utf-8")

def amie_read_binary(file):

    f = open(file, 'rb')

    data = {}
    
    iVals = read_record_long(f, 3)
    data["nLats"] = iVals[0]
    data["nMlts"] = iVals[1]
    data["nTimes"] = iVals[2]

    coLats = read_record_float(f, data["nLats"])
    mlts = read_record_float(f, data["nMlts"])

    data["lats"] = 90.0 - coLats
    data["mlts"] = mlts

    data["nVars"] = read_record_long(f, 1)[0]

    data["Vars"] = []

    for i in np.arange(data["nVars"]):
        v = read_record_string(f).strip()
        data["Vars"].append(v)
        data[v] = []
    
    data["times"] = []
    data["imf"] = []
    data["ae"] = []
    data["dst"] = []
    data["hp"] = []
    data["cpcp"] = []
    
    nPts = data["nLats"] * data["nMlts"]

    for iT in np.arange(data["nTimes"]):

        (iStep, yy, mm, dd, hh, mn)= read_record_long(f, 6)
        data["times"].append(datetime(yy, mm, dd, hh, mn, 0, 0))

        indices = read_record_float(f, 13)
        data["imf"].append(indices[0:4])
        data["ae"].append(indices[4:8])
        data["dst"].append(indices[8:10])
        data["hp"].append(indices[10:12])
        data["cpcp"].append(indices[12])

        for iVar in np.arange(data["nVars"]):
            vals = read_record_float(f, nPts).reshape((data["nLats"], data["nMlts"]))
            v = data["Vars"][iVar]
            data[v].append(vals)
            
    data["version"] = read_record_float(f, 1)[0]
        
    f.close()

    return data

def amie_write_binary(file, data):

    fp = FortranFile(file, 'w')

    iVals = [data["nLats"], data["nMlts"], data["nTimes"]]
    
    fp.write_record(np.array(iVals, dtype = np.int32))

    coLats = np.array(90.0 - data["lats"], dtype = np.float32)
    mlts = np.array(data["mlts"], dtype = np.float32)

    fp.write_record(coLats)
    fp.write_record(mlts)

    # Write out variables!

    val = np.array([data["nVars"]], dtype = np.int32)
    fp.write_record(val)
    
    for var in data["Vars"]:
        varpad = var.ljust(30).encode('utf-8')
        fp.write_record(varpad)

    for iT in np.arange(data["nTimes"]):
        iymdhm = np.array([iT, \
                           data["times"][iT].year, \
                           data["times"][iT].month, \
                           data["times"][iT].day, \
                           data["times"][iT].hour, \
                           data["times"][iT].minute], \
                          dtype=np.int32)
        fp.write_record(iymdhm)

        indices = np.array(np.concatenate((data["imf"][iT], \
                                           data["ae"][iT], \
                                           data["dst"][iT], \
                                           data["hp"][iT], \
                                           [data["cpcp"][iT]])), dtype = np.float32)
        fp.write_record(indices)
        
        for iVar in np.arange(data["nVars"]):
            v = data["Vars"][iVar]
            vals = np.array(data[v][iT], dtype = np.float32)
            fp.write_record(vals)

    val = np.array([data["version"]], dtype = np.float32)
    fp.write_record(val)
            
    fp.close()
    
def amie_write_netcdf(file, data):
    
    ds = nc.Dataset(file, 'w', format = 'NETCDF4')

    cTime = 'Time'

    cLat = 'MagneticLatitude'
    cMlt = 'MagneticLocalTime'
    cPot = 'Potential'
    cEFlux = 'EFlux'
    cAveE = 'AveE'
    cBx = 'IMFBx'
    cBy = 'IMFBy'
    cBz = 'IMFBz'
    cFloat4 = 'f4'
    cInt2 = 'i2'

    time = ds.createDimension(cTime, None)
    lat = ds.createDimension(cLat, data["nLats"])
    mlt = ds.createDimension(cMlt, data["nMlts"])

    timesNC = ds.createVariable(cTime, 'f8', (cTime, ))
    timesNC.units = 's'

    latsNC = ds.createVariable(cLat, cFloat4, (cLat, ))
    latsNC.units = 'deg'
    mltsNC = ds.createVariable(cMlt, cFloat4, (cMlt, ))
    mltsNC.units = 'hours'

    potentialNC = ds.createVariable(cPot, cFloat4, (cTime, cLat, cMlt, ))
    potentialNC.units = 'Volts'

    efluxNC = ds.createVariable(cEFlux, cFloat4, (cTime, cLat, cMlt, ))
    efluxNC.units = 'ergs/cm2/s'
    aveeNC = ds.createVariable(cAveE, cFloat4, (cTime, cLat, cMlt, ))
    aveeNC.units = 'keV'

    latsNC[:] = data["lats"]
    mltsNC[:] = data["mlts"]

    timearray = [] 
    time_1965 = datetime(1965, 1, 1, 0, 0, 0)
    for i in np.arange(0,data["nTimes"]):
        time_delta = data["times"][i] - time_1965
        timearray.append(time_delta.total_seconds())

        cPot = data["Vars"][0]
        potentialNC[i,:,:] = data[cPot][i]

        cPot = data["Vars"][1]
        efluxNC[i,:,:] = data[cPot][i]

        cPot = data["Vars"][2]
        aveeNC[i,:,:] = data[cPot][i]
        
    timesNC[:] = timearray

    ds.close()
        

        
