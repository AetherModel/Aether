#!/usr/bin/env python

import datetime as dt
import numpy as np
from datetime import datetime
from datetime import timedelta
from amie_routines import *

# ------------------------------------------------------------------------
# make the OCFLB
# ------------------------------------------------------------------------

def make_ocflb(mlts, by, bz, DoAddVar, phase):

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

    # Add some variability
    by_deflection_amp = 0.0
    
    if (DoAddVar):
        magnitude_by = 1.5 * np.abs(np.sin(np.arctan2(by, bz)))
        if (by > 0.0):
            loc = 8.0
        else:
            loc = 16.0
        HalfWidth = 3.0

        Period = 1000.0 / (111.0 * np.cos(ocflb0*np.pi/180.0) * 15.0)
        
        by_deflection_amp = magnitude_by * np.cos(2*np.pi * mlts/Period + phase)
        by_deflection_amp = by_deflection_amp * np.exp(-np.abs(mlts-loc)/HalfWidth)
        
    ocflb = ocflb0 \
        - ocflb_lat_deflection_amp * np.cos(mltsR) \
        + by_deflection_amp

    return ocflb
    
# ------------------------------------------------------------------------
# make the aurora
# ------------------------------------------------------------------------

def make_electron_aurora(mlts, lats, by, bz, ocflbBase):

    amp_eflux = 2.0
    amp_avee = 2.5

    tau_eflux = 1.0
    tau_avee = 3.0

    mltsR = mlts * np.pi / 12.0
    r = 90.0 - lats
    t2d, r2d = np.meshgrid(mltsR, r)
    eflux = amp_eflux * (np.cos(t2d)+2.0)/3.0
    avee = amp_avee * (np.cos(t2d)+5.0)/6.0

    tau_eflux_mlts = tau_eflux * (0.5 + 1.5 * (np.cos(mltsR)+1.0))/3.5
    
    nLats = len(lats)
    tau = np.zeros(nLats)
    for i, ocflb in enumerate(ocflbBase):
        tau_ef = tau_eflux_mlts[i]
        dist = lats - (ocflb - tau_ef)
        tau[dist >= 0.0] = tau_eflux_mlts[i]
        tau[dist < 0.0] = tau_eflux_mlts[i] * 2.0
        fac = np.exp(-abs(dist**4/tau**4))
        eflux[:,i] = eflux[:,i] * fac

        fac = np.exp(-abs(dist**2/tau_avee**2))
        avee[:,i] = avee[:,i] * fac
        
    return eflux, avee
    
# ------------------------------------------------------------------------
# make potential
# ------------------------------------------------------------------------

def make_potential(mlts, lats, by, bz, ocflbBase):

    ocflb0 = np.mean(ocflbBase)
    
    # Do everything is kV, then transform to V at end
    amp_bz = -10.0
    amp_by = -5.0
    by_sharpen = (25.0 - np.abs(by)) / 25.0
    tau_potential_base = 7.0 * np.sqrt(by_sharpen)

    mltsR = mlts * np.pi / 12.0

    nMlts = len(mlts)
    nLats = len(lats)
    
    #potential = np.zeros([nMlts, nLats])
    potential = np.zeros([nLats, nMlts])

    # By drives a single cell that is centered at a location that
    # moves as a function of By. So, need to:
    # 1. Redefine the grid in cartesian space:
    r = 90.0 - lats
    t2d, r2d = np.meshgrid(mltsR, r)
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

        potential_visc[:,i] = potential_visc[:,i] * fac

        if ((by > 0) & (mlts[i] > 12)):
            fac = fac * by_sharpen
        if ((by < 0) & (mlts[i] < 12)):
            fac = fac * by_sharpen
        potential_bz[:,i] = potential_bz[:,i] * fac

        fac = np.exp(-2.0 * dist**2 / tau**2)
        potential_harang[:,i] = potential_harang[:,i] * fac
        
    potential = (potential_by + \
                 potential_bz + \
                 potential_harang + \
                 potential_visc) * 1000.0

    return potential


# ------------------------------------------------------------------------
# main code
# ------------------------------------------------------------------------

DoAddVarOCFLB = 0

dt = 1.0/60.0
subtimes = np.arange(12,13+dt,dt)
nTimes = len(subtimes)
bysub = np.random.normal(-5.0,2.0,nTimes)
bzsub = bysub * 0.0 - 2.0

bz = np.concatenate([[-2.0], bzsub, [0.0]])
by = np.concatenate([[-5.0], bysub,  [-5.0]])
alltimes = np.concatenate([[0.0], subtimes, [24.0]])

time_1965 = datetime(1965, 1, 1, 0, 0, 0)
basetime = datetime(2020, 3, 20, 0, 0, 0)
alltimes = np.array(alltimes) * 3600.0

nTimes = len(bz)

dLat = 0.25
minLat = 40.0
nLats = (90.0 - minLat)/dLat + 1
dMlt = 0.25
nMlts = (24.0-0.0)/dMlt + 1

lats = np.arange(minLat, 90.0+dLat, dLat)
mlts = np.arange(0.0, 24.0+dMlt, dMlt)


data = {}

data["nLats"] = nLats
data["nMlts"] = nMlts
data["nTimes"] = nTimes

data["lats"] = lats
data["mlts"] = mlts

data["times"] = []
for t in alltimes:
    data["times"].append(basetime + timedelta(seconds = t))

data["imf"] = []
data["ae"] = []
data["dst"] = []
data["hp"] = []
data["cpcp"] = []

data["Vars"] = ['Potential (kV)', \
                'Total Energy Flux (ergs/cm2/s)', \
                'Mean Energy (ergs)']

data["nVars"] = len(data["Vars"])

data["version"] = 1.2

for var in data["Vars"]:
    data[var] = []

vel_km = 0.4 # km/s
vel_rad = vel_km / (111.0 * np.cos(70.0*np.pi/180.0)) * np.pi / 180.0

for i in np.arange(0,nTimes):

    time_current = basetime + timedelta(seconds = alltimes[i])
    data["times"].append(time_current)

    time_delta = time_current - time_1965
    ut = time_delta.total_seconds() % 86400.0
    phase = vel_rad * ut
    
    ocflb = make_ocflb(mlts, by[i], bz[i], DoAddVarOCFLB, phase)
    pot2d = make_potential(mlts, lats, by[i], bz[i], ocflb)
    eflux2d,avee2d = make_electron_aurora(mlts, lats, by[i], bz[i], ocflb)

    #  ; IMF should be (nTimes,4) - V, Bx, By, Bz
    #  ; AE  (nTimes, 4) - AL, AU, AE, AEI?
    #  ; Dst (nTimes, 2) - Dst, Dsti?
    #  ; Hpi (nTimes, 2) - HP, Joule Heating (GW)
    #  ; CPCP (nTimes) - CPCP (kV)

    data["imf"].append([400.0, 0.0, by[i], bz[i]])
    # These are all bad values for now....
    data["ae"].append([0.0, 0.0, 0.0, 0.0])
    data["dst"].append([0.0, 0.0])
    data["hp"].append([0.0, 0.0])
    data["cpcp"].append((np.max(pot2d) - np.min(pot2d))/1000.0)

    cPot = data["Vars"][0]
    data[cPot].append(pot2d)
    
    cEFlux = data["Vars"][1]
    data[cEFlux].append(eflux2d)
    
    cAveE = data["Vars"][2]
    data[cAveE].append(avee2d)

# Write out in new netCDF format:
file = 'test_ed.nc'
amie_write_netcdf(file, data)

# Write out in old AMIE format:
file = 'test_ed.bin'
amie_write_binary(file, data)
