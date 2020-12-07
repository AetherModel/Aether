#!/opt/local/bin/python

# the Scientific Python netCDF 3 interface
# http://dirac.cnrs-orleans.fr/ScientificPython/

from netCDF4 import Dataset
import matplotlib.pyplot as plt

# the 'classic' version of the netCDF4 python interface
# http://code.google.com/p/netcdf4-python/
#from netCDF4_classic import Dataset
import numpy as np
from pylab import cm
import glob

# open netCDF file for reading.

filelist = sorted(glob.glob("3DALL*.nc"))
file = filelist[-1];

ncfile = Dataset(file,'r') 

for var in ncfile.variables.values():
    print(var)

lats = np.array(ncfile.variables['Latitude'])*180.0/3.14159
lons = np.array(ncfile.variables['Longitude'])*180.0/3.14159
alts = np.array(ncfile.variables['Altitude'])/1000.0

nLons = len(lats[:,0,0])
nLats = len(lats[0,:,0])
nAlts = len(lats[0,0,:])

temp = np.array(ncfile.variables['Temperature'])
o = np.array(ncfile.variables['O'])

ncfile.close()

iCut = 0

value = temp # np.log10(o)

# Lat/Lon Cut:
if (iCut == 0):
    AllData2D = value[:,:,40]
    xPos = lons[:,0,0]
    yPos = lats[0,:,0]

# Lon/Alt Cut:
if (iCut == 1):
    AllData2D = value[:,int(nLats/2),:]
    xPos = lons[:,0,0]
    yPos = alts[0,0,:]
    
# Lat/Alt Cut:
if (iCut == 2):
    AllData2D = value[nLons/2,:,:]
    xPos = lats[0,:,0]
    yPos = alts[0,0,:]

fig = plt.figure()
ax = fig.add_subplot(111)

norm = cm.colors.Normalize(vmax=np.max(AllData2D), vmin=np.min(AllData2D))
cmap = cm.plasma

d2d = np.transpose(AllData2D)

# sTime = time.strftime('%y%m%d_%H%M%S')
# outfile = file+'_'+sTime+'.png'
outfile = 'test.png'

# Define plot range:
minX = (xPos[ 1] + xPos[ 2])/2
maxX = (xPos[-2] + xPos[-3])/2
minY = (yPos[ 1] + yPos[ 2])/2
maxY = (yPos[-2] + yPos[-3])/2

#cax = ax.pcolor(xPos, yPos, d2d, vmin=mini, vmax=maxi)#, shading='auto')
cax = ax.pcolor(xPos, yPos, d2d, shading='auto')
ax.set_ylim([minY,maxY])
ax.set_xlim([minX,maxX])

cbar = fig.colorbar(cax)

fig.savefig(outfile)
plt.close()
