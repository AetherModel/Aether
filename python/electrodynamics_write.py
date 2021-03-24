#!/usr/bin/env python

import netCDF4 as nc
import datetime as dt
import numpy as np
from datetime import datetime
from datetime import timedelta

# ------------------------------------------------------------------------
# make the OCFLB
# ------------------------------------------------------------------------

def make_ocflb(mlts, bz):

    mltsR = mlts * np.pi / 12.0

    # This is the shift of the ocflb as you go around in MLT, with
    # midnight being deflected by this much towards the equator,
    # and noon being deflected this much towards the pole:
    ocflb_lat_deflection_amp = 3.0

    # Use studies of cusp latitude to figure out where to put OCFLB
    if (bz > 0):
        CuspLat = 77.2 + 0.11 * bz
    else:
        if (bz > -10):
            CuspLat = 77.2 + 1.1 * bz  # bz is negative!
        else:
            CuspLat = 21.7 * np.exp(0.1*bz) + 58.2

    ocflb0 = CuspLat - ocflb_lat_deflection_amp 
    ocflb = ocflb0 - ocflb_lat_deflection_amp * np.cos(mltsR)

    return ocflb
    
# ------------------------------------------------------------------------
# make the aurora
# ------------------------------------------------------------------------

def make_electron_aurora(mlts, lats, by, bz):

    ocflbBase = make_ocflb(mlts, bz)
    
    amp_eflux = 2.0
    amp_avee = 2.5

    tau_eflux = 1.0
    tau_avee = 3.0

    mltsR = mlts * np.pi / 12.0
    r = 90.0 - lats
    r2d, t2d = np.meshgrid(r, mltsR)
    eflux = amp_eflux * (np.cos(t2d)+2.0)/3.0
    avee = amp_avee * (np.cos(t2d)+5.0)/6.0

    tau_eflux_mlts = tau_eflux * (0.5 + 1.5 * (np.cos(mltsR)+1.0))/3.5
    
    nLats = len(lats)
    tau = np.zeros(nLats)
    for i, ocflb in enumerate(ocflbBase):
        tau_ef = tau_eflux_mlts[i]
        dist = lats - ocflb
        tau[dist >= 0.0] = tau_eflux_mlts[i]
        tau[dist < 0.0] = tau_eflux_mlts[i] * 2.0
        fac = np.exp(-abs(dist**4/tau**4))
        eflux[i,:] = eflux[i,:] * fac

        fac = np.exp(-abs(dist**2/tau_avee**2))
        avee[i,:] = avee[i,:] * fac
        
    return eflux, avee
    
# ------------------------------------------------------------------------
# make potential
# ------------------------------------------------------------------------

def make_potential(mlts, lats, by, bz):

    ocflbBase = make_ocflb(mlts, bz)
    ocflb0 = np.mean(ocflbBase)
    
    # Do everything is kV, then transform to V at end
    amp_bz = -10.0
    amp_by = -5.0
    by_sharpen = (25.0 - np.abs(by)) / 25.0
    tau_potential_base = 7.0 * np.sqrt(by_sharpen)

    mltsR = mlts * np.pi / 12.0

    nMlts = len(mlts)
    nLats = len(lats)
    
    potential = np.zeros([nMlts, nLats])

    # By drives a single cell that is centered at a location that
    # moves as a function of By. So, need to:
    # 1. Redefine the grid in cartesian space:
    r = 90.0 - lats
    r2d, t2d = np.meshgrid(r, mltsR)
    x2d = r2d * np.cos(t2d)
    y2d = r2d * np.sin(t2d)

    # Define the Center of the By-driven cell:
    t0 = (12.0 + by/4.0) * np.pi/12.0
    r0 = (90.0 - np.max(ocflbBase))/2.0
    x0 = r0 * np.cos(t0)
    y0 = r0 * np.sin(t0)

    distance_for_by = np.sqrt((x2d-x0)**2 + (y2d-y0)**2)
    tau_for_by = (90.0 - ocflb0)/1.5
    potential_by = by * amp_by * np.exp(-distance_for_by / tau_for_by)

    potential_bz = bz * amp_bz * np.sin(t2d)
    potential_visc = -0.5 * amp_bz * np.sin(t2d)

    shifted = t2d - 1.0 * np.pi / 12.0
    potential_harang = amp_bz * ((np.cos(shifted)+1.0)/2.0)**3.0
    
    tau = np.zeros(nLats)
    for i, ocflb in enumerate(ocflbBase):
        tau[:] = tau_potential_base
        # equatorward of the OCFLB, fall off faster
        dist = lats - ocflb
        tau[dist < 0.0] = tau_potential_base / 1.5
        fac = np.exp(-abs(dist/tau))

        potential_visc[i,:] = potential_visc[i,:] * fac

        if ((by > 0) & (mlts[i] > 12)):
            fac = fac * by_sharpen
        if ((by < 0) & (mlts[i] < 12)):
            fac = fac * by_sharpen
        potential_bz[i,:] = potential_bz[i,:] * fac

        fac = np.exp(-2.0 * dist**2 / tau**2)
        potential_harang[i,:] = potential_harang[i,:] * fac
        
    potential = (potential_by + \
                 potential_bz + \
                 potential_harang + \
                 potential_visc) * 1000.0

    return potential


# ------------------------------------------------------------------------
# main code
# ------------------------------------------------------------------------

file = 'test.nc'

bz = [-2.0, -2.0, -2.0, 0.0]
by = [-5.0, 0.0,  5.0,  5.0]

time_1965 = datetime(1965, 1, 1, 0, 0, 0)

basetime = datetime(2020, 3, 20, 0, 0, 0)
alltimes = np.array([0.0, 12.0, 12.1667, 24.0]) * 3600.0

nTimes = len(bz)

dLat = 1.0
minLat = 40.0
nLats = (90.0 - minLat)/dLat + 1
dMlt = 0.5
nMlts = (24.0-0.0)/dMlt + 1

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
lat = ds.createDimension(cLat, nLats)
mlt = ds.createDimension(cMlt, nMlts)

timesNC = ds.createVariable(cTime, 'f8', (cTime, ))
timesNC.units = 's'

latsNC = ds.createVariable(cLat, cFloat4, (cLat, ))
latsNC.units = 'deg'
mltsNC = ds.createVariable(cMlt, cFloat4, (cMlt, ))
mltsNC.units = 'hours'

potentialNC = ds.createVariable(cPot, cFloat4, (cTime, cMlt, cLat, ))
potentialNC.units = 'Volts'

efluxNC = ds.createVariable(cEFlux, cFloat4, (cTime, cMlt, cLat, ))
efluxNC.units = 'ergs/cm2/s'
aveeNC = ds.createVariable(cAveE, cFloat4, (cTime, cMlt, cLat, ))
aveeNC.units = 'keV'

lats = np.arange(minLat, 90.0+dLat, dLat)
mlts = np.arange(0.0, 24.0+dMlt, dMlt)

latsNC[:] = lats
mltsNC[:] = mlts

timearray = [] 
for i in np.arange(0,nTimes):
    pot2d = make_potential(mlts, lats, by[i], bz[i])
    eflux2d,avee2d = make_electron_aurora(mlts, lats, by[i], bz[i])
    potentialNC[i,:,:] = pot2d
    efluxNC[i,:,:] = eflux2d
    aveeNC[i,:,:] = avee2d
    time_current = basetime + timedelta(seconds = alltimes[i])
    time_delta = time_current - time_1965
    timearray.append(time_delta.total_seconds())

timesNC[:] = timearray
print(timearray)

ds.close()
