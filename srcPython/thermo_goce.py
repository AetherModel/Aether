#!/usr/bin/env python

from datetime import datetime
from datetime import timedelta
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
from gitm_routines import *

def read_goce(file):

    data = {}
    data["times"] = []
    data["alts"] = []
    data["lats"] = []
    data["lons"] = []
    data["lst"] = []
    data["density"] = []
    data["Ue"] = []
    data["Un"] = []
    data["Ur"] = []
    data["densityError"] = []
    data["windError"] = []
    data["FlagOver"] = []
    data["FlagEclipse"] = []
    data["FlagAD"] = []
    data["FlagThuster"] = []

    f = open(file, 'r')

    for line in f:

        if (line.find('#') < 0):
            items = line.split()
            ymd = items[0].split('-')
            hms = items[1].split(':')
            s = float(hms[2])
            data["times"].append(datetime(int(ymd[0]),int(ymd[1]),int(ymd[2]),
                                         int(hms[0]),int(hms[1]),int(s)))
            data["alts"].append(float(items[3])/1000.0)
            data["lons"].append(float(items[4]))
            data["lats"].append(float(items[5]))
            data["lst"].append(float(items[6]))
            data["density"].append(float(items[8]))
            data["Ue"].append(float(items[9]))
            data["Un"].append(float(items[10]))
            data["Ur"].append(float(items[11]))
            data["densityError"].append(float(items[12]))
            data["windError"].append(float(items[13]))
            data["FlagOver"].append(int(items[14]))
            data["FlagEclipse"].append(int(items[15]))
            data["FlagAD"].append(int(items[16]))
            data["FlagThuster"].append(int(items[17]))

    f.close()

    return data

def find_index(time, t):

    if (t < time[0]):
        return 0
    if (t > time[-1]):
        return len(time)-1

    iLow = 0
    iHigh = len(time)
    iMid = int((iHigh + iLow)/2)
    while (iHigh - iLow > 1):
        if (time[iMid] == t):
            iHigh = iMid
            iLow = iMid
        else:
            if (t > time[iMid]):
                iLow = iMid
            else:
                iHigh = iMid
            iMid = int((iHigh + iLow)/2)
    return iMid
            
def smooth_data(time, data, window):

    smoothed = []
    for i, t in enumerate(time):
        iMin = find_index(time, t-window/2)
        iMax = find_index(time, t+window/2)
        s = np.mean(data[iMin : iMax+1])
        smoothed.append(s)

    return smoothed
    

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Main Code!
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

rtod =  180.0 / np.pi
dtor = 1.0 / rtod
headers = read_gitm_headers('3DALL')

nGitmFiles = len(headers["time"])

# get lons, lats, alts:
vars = [0,1,2]
data = read_gitm_one_file(headers["filename"][0], vars)
Alts = data[2][0][0]/1000.0;
Lons = data[0][:,0,0]*rtod;
Lats = data[1][0,:,0]*rtod;
[nLons, nLats, nAlts] = data[0].shape

dLon = Lons[1]-Lons[0]
dLat = Lats[1]-Lats[0]

goceDir_lap = '/Users/ridley/Data/Goce/timeseries_data/'
goceDir_zed = '/data1/Data/Goce/timeseries_data/'

found = False
if (os.path.isdir(goceDir_lap)):
    goceDir = goceDir_lap
    found = True
if (os.path.isdir(goceDir_zed)):
    goceDir = goceDir_zed
    found = True
if (not found):
    print("Can't seem to find GOCE directory...")
    exit()

goceFileFront = 'goce_denswind_ac082_v2_0_'
goceYear = headers["time"][int(nGitmFiles/2)].strftime('%Y-%m.txt')
goceFile = goceDir + goceFileFront + goceYear

if (not os.path.exists(goceFile)):
    print("Can not find GOCE file : ", goceFile)
    exit()

data = read_goce(goceFile)
Vr = np.array(data["Ur"])
Ve = np.array(data["Ue"])
Vn = np.array(data["Un"])
rho = np.array(data["density"])
goceAltsAll = np.array(data["alts"])
goceLatsAll = np.array(data["lats"]) 
goceLstAll = np.array(data["lst"]) 

re = 6378.137
rp = 6356.75
diff = re - rp
rgitm = 6372.0
r = rp + diff * np.cos(goceLatsAll * dtor)
correction = r - rgitm
goceAltsPrime = goceAltsAll + correction

goceRm = (r + goceAltsAll) * 1000.0
meanR = np.mean(goceRm)

G = 6.67259e-11
Me = 5.9722e24
mu = G * Me
v0 = np.sqrt(mu / meanR)
period = 2.0 * np.pi * meanR / v0


iBefore0 = -1
iAfter0 = -1
# here we just assume that we will start with the first set of files.
# this shouldn't matter...
iBefore = 0
iAfter = 1

vars = [3]
nVars = len(vars)

AfterVals = np.zeros(nVars)
BeforeVals = np.zeros(nVars)

gitmRho = []
goceRho = []
goceTime = []
goceLats = []
goceLst = []

for i, time in enumerate(data["times"]):

    while (time > headers["time"][iAfter]):
        iAfter = iAfter+1
        if (iAfter >= nGitmFiles-1):
            break
        
    if (iAfter == nGitmFiles):
        break

    iBefore = iAfter-1

    if (iBefore != iBefore0):
        file = headers["filename"][iBefore]
        BeforeData = read_gitm_one_file(file, vars)
    if (iAfter != iAfter0):
        file = headers["filename"][iAfter]
        AfterData = read_gitm_one_file(file, vars)

    if (time >= headers["time"][iBefore]):
            
        dt = (headers["time"][iAfter] - \
              headers["time"][iBefore]).total_seconds()
        xt = (time - headers["time"][iBefore]).total_seconds() / dt

        lon = data["lons"][i]
        lat = goceLatsAll[i]
        alt = goceAltsPrime[i]
        goceRho.append(rho[i])
        goceTime.append(time)
        goceLats.append(lat)
        goceLst.append(goceLstAll[i])
        
        xLon = (lon-Lons[0])/dLon
        iLon = int(xLon)
        xLon = xLon - iLon
        
        yLat = (lat-Lats[0])/dLat
        jLat = int(yLat)
        yLat = yLat - jLat

        kAlt = 0
        zAlt = 0.0
        if ((alt > Alts[0]) and (nAlts > 1)):
            if (alt > Alts[nAlts-1]):
                # above domain:
                kAlt = nAlts-2
                zAlt = 1.0
            else:
                while (Alts[kAlt] < alt):
                    kAlt = kAlt + 1
                kAlt = kAlt - 1
                zAlt = (alt - Alts[kAlt]) / (Alts[kAlt+1] - Alts[kAlt])
            kAltp1 = kAlt + 1
        else:
            kAltp1 = kAlt

        for i, v in enumerate(vars):
            BeforeVals[i] = \
                (1-xLon)*(1-yLat)*(1-zAlt)*BeforeData[v][iLon][jLat][kAlt]+\
                (  xLon)*(1-yLat)*(1-zAlt)*BeforeData[v][iLon+1][jLat][kAlt]+\
                (1-xLon)*(  yLat)*(1-zAlt)*BeforeData[v][iLon][jLat+1][kAlt]+\
                (  xLon)*(  yLat)*(1-zAlt)*BeforeData[v][iLon+1][jLat+1][kAlt]+\
                (1-xLon)*(1-yLat)*(  zAlt)*BeforeData[v][iLon][jLat][kAltp1]+\
                (  xLon)*(1-yLat)*(  zAlt)*BeforeData[v][iLon+1][jLat][kAltp1]+\
                (1-xLon)*(  yLat)*(  zAlt)*BeforeData[v][iLon][jLat+1][kAltp1]+\
                (  xLon)*(  yLat)*(  zAlt)*BeforeData[v][iLon+1][jLat+1][kAltp1]
            AfterVals[i] = \
                (1-xLon)*(1-yLat)*(1-zAlt)*AfterData[v][iLon][jLat][kAlt]+\
                (  xLon)*(1-yLat)*(1-zAlt)*AfterData[v][iLon+1][jLat][kAlt]+\
                (1-xLon)*(  yLat)*(1-zAlt)*AfterData[v][iLon][jLat+1][kAlt]+\
                (  xLon)*(  yLat)*(1-zAlt)*AfterData[v][iLon+1][jLat+1][kAlt]+\
                (1-xLon)*(1-yLat)*(  zAlt)*AfterData[v][iLon][jLat][kAltp1]+\
                (  xLon)*(1-yLat)*(  zAlt)*AfterData[v][iLon+1][jLat][kAltp1]+\
                (1-xLon)*(  yLat)*(  zAlt)*AfterData[v][iLon][jLat+1][kAltp1]+\
                (  xLon)*(  yLat)*(  zAlt)*AfterData[v][iLon+1][jLat+1][kAltp1]
            
        date = time.strftime('%Y %m %d %H %M %S 00 ')
        pos = '%7.2f %7.2f %8.2f' % (lon, lat, alt)
        vals = ''
        i = 0
        gitmRho.append((1-xt) * BeforeVals[i] + xt * AfterVals[i])
        for i, v in enumerate(vars):
            v = (1-xt) * BeforeVals[i] + xt * AfterVals[i]
            vals = vals + '  %e' % v
        #fpout.write(date+pos+vals+'\n')

    iBefore0 = iBefore
    iAfter0 = iAfter


dt = ((goceTime[-1] - goceTime[0]).total_seconds())/3600.0
if (dt < 2*86400.0):
    StartTime = datetime(goceTime[0].year, \
                         goceTime[0].month, \
                         goceTime[0].day, \
                         goceTime[0].hour)
else:
    StartTime = datetime(goceTime[0].year, \
                         goceTime[0].month, \
                         goceTime[0].day)
    
sTime = StartTime.strftime('%b %d, %Y %H:%M UT') + ' - ' + \
    goceTime[-1].strftime('%b %d, %Y %H:%M UT (Hours)')

tInS = []
for t in goceTime:
    tInS.append((t - StartTime).total_seconds())
tInS = np.array(tInS)

EndTimeS = tInS[-1]
dLat = 1.0
nX = int(np.round(EndTimeS/period))+1
nY = int(180.0/dLat)

gitmRhoMapa = np.zeros([nX, nY]) + np.nan
goceRhoMapa = np.zeros([nX, nY]) + np.nan

gitmRhoMapd = np.zeros([nX, nY]) + np.nan
goceRhoMapd = np.zeros([nX, nY]) + np.nan

mapTime = []
mapTime1D = []
mapLats = np.zeros([nX+1, nY+1])

for i in range(nX+1):
    mapLats[i,:] = np.arange(nY+1) * dLat - 90.0 #+ dLat/2.0
    dt = i * period - period/2
    temp = [ dt ] * (nY+1)
    mapTime.append(temp)
    mapTime1D.append(StartTime + timedelta(seconds = dt))
    
nPts = len(goceTime)

goceLats = np.array(goceLats)
mapTime = np.array(mapTime)/3600.0

ix = np.round(tInS/period)
iy = np.round((goceLats + 90.0)/dLat)

dLat = goceLats[1:] - goceLats[0:-1]


all = [gitmRho, goceRho]
maxi = np.max(all)
fac = 10**(int(np.log10(maxi))-1)
gitmRho = np.array(gitmRho) / fac
goceRho = np.array(goceRho) / fac

asc = 0.0
ascN = 0
des = 0.0
desN = 0
for iPt in range(nPts-1):
    i = int(ix[iPt])
    j = int(iy[iPt])
    if (dLat[iPt] > 0):
        gitmRhoMapa[i][j] = gitmRho[iPt]
        goceRhoMapa[i][j] = goceRho[iPt]
        if (np.abs(goceLats[iPt]) < 30.0):
            asc = asc + goceLst[iPt]
            ascN = ascN + 1
    else:
        gitmRhoMapd[i][j] = gitmRho[iPt]
        goceRhoMapd[i][j] = goceRho[iPt]
        if (np.abs(goceLats[iPt]) < 30.0):
            des = des + goceLst[iPt]
            desN = desN + 1
        
all = [gitmRho, goceRho]
mini = np.min(all)
maxi = np.max(all)

asc = asc / ascN
des = des / desN

# ----------------------------------------------------------------------------
# Ascending Rho

fig = plt.figure(figsize = (10,10))

ax = fig.add_axes([0.075, 0.53, 1.0, 0.43])
cax = ax.pcolor(mapTime, mapLats, goceRhoMapa, vmin = mini, vmax = maxi)
ax.set_xlabel(sTime)
ax.set_ylabel('Latitude (deg)')
ax.set_ylim([-82.0, 82.0])
ax.set_xlim([np.min(mapTime), np.max(mapTime)])
ax.set_title("GOCE (top) vs GITM (bot) - Ascending Node (LT : %4.1f hours)" % asc)
cbar = fig.colorbar(cax, ax=ax, shrink = 0.75, pad=0.02)
cbar.set_label('GOCE Rho (%7.1e kg/m3)' % fac,rotation=90)


ax2 = fig.add_axes([0.075, 0.05, 1.0, 0.43])
cax2 = ax2.pcolor(mapTime, mapLats, gitmRhoMapa, vmin = mini, vmax = maxi)
ax2.set_xlabel(sTime)
ax2.set_ylabel('Latitude (deg)')
ax2.set_ylim([-82.0, 82.0])
ax2.set_xlim([np.min(mapTime), np.max(mapTime)])
cbar2 = fig.colorbar(cax2, ax=ax2, shrink = 0.75, pad=0.02)
cbar2.set_label('GITM Rho (%7.1e kg/m3)' % fac,rotation=90)

outfile = StartTime.strftime('ascending_%Y%m%d.png')
print('Writing file : ' + outfile)
plt.savefig(outfile)
plt.close()

# ----------------------------------------------------------------------------
# Descending Rho

fig = plt.figure(figsize = (10,10))

ax = fig.add_axes([0.075, 0.53, 1.0, 0.43])
cax = ax.pcolor(mapTime, mapLats, goceRhoMapd, vmin = mini, vmax = maxi)
ax.set_xlabel(sTime)
ax.set_ylabel('Latitude (deg)')
ax.set_ylim([-82.0, 82.0])
ax.set_xlim([np.min(mapTime), np.max(mapTime)])
ax.set_title("GOCE (top) vs GITM (bot) - Descending Node (LT : %4.1f hours)" % des)
cbar = fig.colorbar(cax, ax=ax, shrink = 0.75, pad=0.02)
cbar.set_label('GOCE Rho (%7.1e kg/m3)' % fac,rotation=90)


ax2 = fig.add_axes([0.075, 0.05, 1.0, 0.43])
cax2 = ax2.pcolor(mapTime, mapLats, gitmRhoMapd, vmin = mini, vmax = maxi)
ax2.set_xlabel(sTime)
ax2.set_ylabel('Latitude (deg)')
ax2.set_ylim([-82.0, 82.0])
ax2.set_xlim([np.min(mapTime), np.max(mapTime)])
cbar2 = fig.colorbar(cax2, ax=ax2, shrink = 0.75, pad=0.02)
cbar2.set_label('GITM Rho (%7.1e kg/m3)' % fac,rotation=90)

outfile = StartTime.strftime('descending_%Y%m%d.png')
print('Writing file : ' + outfile)
plt.savefig(outfile)
plt.close()

# ----------------------------------------------------------------------------
# Line plots

goceRhoSmoothed = smooth_data(tInS, goceRho, period)
gitmRhoSmoothed = smooth_data(tInS, gitmRho, period)

sTime = goceTime[0].strftime('%b %d, %Y %H:%M UT') + ' - ' + \
    goceTime[-1].strftime('%b %d, %Y %H:%M UT')

fig = plt.figure(figsize = (10,10))
ax = fig.add_axes([0.075, 0.53, 0.88, 0.43])
ax.plot(goceTime, gitmRho, 'b', alpha = 0.1)
ax.plot(goceTime, gitmRhoSmoothed, 'b', label = 'GITM')
ax.plot(goceTime, goceRho, 'r', alpha = 0.1)
ax.plot(goceTime, goceRhoSmoothed, 'r', label = 'GOCE')
ax.set_ylabel('Rho (kg/m3)')
ax.set_xlim(goceTime[0], goceTime[-1])
ax.legend()

ax2 = fig.add_axes([0.075, 0.05, 0.88, 0.43])
rawDiff = (np.array(gitmRho) - np.array(goceRho)) / np.array(goceRho) * 100.0
smoothedDiff = (np.array(gitmRhoSmoothed) - np.array(goceRhoSmoothed)) / np.array(goceRhoSmoothed) * 100.0
ax2.plot(goceTime, rawDiff, 'k', alpha = 0.1)
ax2.plot(goceTime, smoothedDiff, 'k')
ax2.plot(goceTime, rawDiff*0.0, 'g:')
ax2.set_ylabel('GITM - GOCE Rho Diff (%)')
ax2.set_xlim(goceTime[0], goceTime[-1])
ax2.set_xlabel(sTime)

outfile = StartTime.strftime('goce_gitm_rho_line_%Y%m%d.png')
print('Writing file : ' + outfile)
plt.savefig(outfile)
plt.close()

exit()


                                        
