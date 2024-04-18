#!/usr/bin/env python

import datetime as dt
import numpy as np
from datetime import datetime
from datetime import timedelta
from amie_routines import *
from omniweb import *
from kyoto import *
import sys
import matplotlib.pyplot as plt
import matplotlib.dates as dates

# ------------------------------------------------------------------------
# These are stolen from Angeline's code
# ------------------------------------------------------------------------

def bool_string(line):
    """ Determine whether a string should be True or False
    Parameters
    ----------
    line : string
        Line to be tested
    Returns
    -------
    bout : bool
        Boolean output (True/False)
    Raises
    ------
    ValueError
        If the value cannot be interpreted as True or False
    Notes
    -----
    Accepts empty, true, false, t, f, 1, and 0 in any capitalization combo
    """

    line = line.lower()

    bout = None
    if line in ['true', 't', '1', '']:
        bout = True
    elif line in ['false', 'f', '0']:
        bout = False

    if bout is None:
        raise ValueError('input not interpretable as a boolean')

    return bout


def none_string(line):
    """ Determine whether a string should be None
    Parameters
    ----------
    line : string
        Line to be tested
    Returns
    -------
    out_line : string or NoneType
        None if all-lowercase version of line is "none" or line is zero length.
        Otherwise returns original value of line
    """
    out_line = None if line.lower() == "none" or len(line) == 0 else line
    return out_line


def process_command_line_input():
    """ Process command line input, needed to possible ipython use
    Returns
    -------
    input_args : list
        List of input arguements
    """

    input_args = sys.argv
    if input_args[0].find('ipython') >= 0:
        input_args = list()
    else:
        input_args.pop(0)

    return input_args

# ------------------------------------------------------------------------
# This sets the default values of the inputs and parses the inputs
# This is adapted from Angeline's code.
# ------------------------------------------------------------------------

def get_command_line_args(argv):
    """ Parse the arguements and set to a dictionary
    Parameters
    ----------
    argv : list
        List of arguments fed on the command line
    Returns
    -------
    args : dict
        A dictionary containing information about arguements, including:

    """
    # Initialize the arguments to their default values   

    args = {'startdate': '20200101',
            'enddate': '20200102',
            'outfile': 'test.nc',
            'dt': 5,
            'real': True,
            'south': False,
            'tcv': False,
            'substorm': False,
            'ions': False,
            'move': False,
            'cusp': False}

    arg_type = {'startdate': str,
                'enddate': str,
                'outfile': str,
                'dt': float,
                'real': bool,
                'south': bool,
                'tcv': bool,
                'substorm': bool,
                'ions': bool,
                'move': bool,
                'cusp': bool}
    
    # If there is input, set default help to False
    args['help'] = False if len(argv) > 0 else True
    
    # Cycle through all arguments except the first, saving input
    for arg in argv:
        # Treat the file list and formatting seperately
        if arg.find('-') == 0:
            # This is not a filename, remove the dash to get the key
            split_arg = arg.split('=')
            akey = split_arg[0][1:]
            # Get the argument value as the desired type
            if akey not in arg_type.keys():
                raise ValueError(''.join(['unknown command line input, ',
                                          arg, ', try -help for details']))

            if len(split_arg) == 1:
                if arg_type[akey] == bool:
                    arg_val = True
                else:
                    raise ValueError('expected equality after flag {:}'.format(
                        akey))
            else:
                if arg_type[akey] == int:
                    arg_val = int(split_arg[1])
                elif arg_type[akey] == float:
                    arg_val = float(split_arg[1])
                elif arg_type[akey] == str:
                    arg_val = split_arg[1]
                else:
                    # This is boolean input
                    arg_val = bool_string(split_arg[1])

            args[akey] = arg_val
                    
    return args
    
# ------------------------------------------------------------------------
# make the OCFLB
# ------------------------------------------------------------------------

def make_ocflb(mlts, by, bz, DoAddVar, phase, substorm, move):

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
        HalfWidthMLT = 3.0

        Period = 1000.0 / (111.0 * np.cos(ocflb0*np.pi/180.0) * 15.0)
        
        by_deflection_amp = magnitude_by * np.cos(2*np.pi * mlts/Period + phase)
        by_deflection_amp = by_deflection_amp * np.exp(-np.abs(mlts-loc)/HalfWidthMLT)

    if (substorm > -1000.0):
        # let's push the OCFLB poleward a bit during the substorm
        offset_amp = np.abs(substorm) * ocflb_lat_deflection_amp * 2.0
        # The center of the push should move from 23 to 21 MLT during expansion:
        if ((substorm < 0.0) and (move)):
            x = np.abs(substorm)
            center = (x * 22.0 + (1.0-x) * 24.0)
        else:
            center = 22.0
        # print("substorm: ", substorm, center)
        ss_deflection = offset_amp * \
            ((np.cos((mlts-center) * np.pi/12.0) + 1)/2.0)**6
    else:
        ss_deflection = 0.0
        
    ocflb = ocflb0 \
        - ocflb_lat_deflection_amp * np.cos(mltsR) \
        + by_deflection_amp \
        + ss_deflection

    return ocflb
    
# ------------------------------------------------------------------------
# Define characteristics of a substorm
# ------------------------------------------------------------------------

def define_substorm_characteristics(mlts, lats, ae, ocflbBase, substorm, IsPot, move):

    mltsR = mlts * np.pi / 12.0
    r = 90.0 - lats
    t2d, r2d = np.meshgrid(mltsR, r)


#    if (DoAddVar):
#        magnitude = 1.5 * np.abs(np.sin(np.arctan2(by, bz)))
#        if (by > 0.0):
#            loc = 8.0
#        else:
#            loc = 16.0
#        HalfWidthMLT = 3.0
#
#        Period = 1000.0 / (111.0 * np.cos(ocflb0*np.pi/180.0) * 15.0)
#        
#        by_deflection_amp = magnitude_by * np.cos(2*np.pi * mlts/Period + phase)
#        by_deflection_amp = by_deflection_amp * np.exp(-np.abs(mlts-loc)/HalfWidthMLT)

    
    if (substorm > -1000):
    
        tau_eflux = 2.0 + ae/100.0 * 0.35 # how the thickness of the aurora grows

        if ((substorm < 0.0) and (move)):
            x = np.abs(substorm)
            center = (x * 22.0 + (1.0-x) * 24.0) * np.pi / 12.0
        else:
            center = 22.0 * np.pi / 12.0

        diff = (t2d - center) 
        value = ((np.cos(diff) + 1)/2.0)**10
        ind = value > 0.5
        value[ind] = 0.5
        value = value * 2.0

        tau_eflux_mlts = tau_eflux * (0.5 + 1.5 * (np.cos(mltsR)+1.0))/3.5

        # shift equatorward
        shift_mlts = -tau_eflux * 0.75
        nLats = len(lats)
        tau = np.zeros(nLats)
        for i, ocflb in enumerate(ocflbBase):
            tau_ef = tau_eflux_mlts[i]
            shift = shift_mlts
            dist = lats - (ocflb + shift)
            tau[dist >= 0.0] = tau_eflux_mlts[i] * (1.0 + IsPot * 1.00)
            tau[dist < 0.0] = tau_eflux_mlts[i] * (0.75 + IsPot * 0.50)
            fac = np.exp(-abs(dist**4/tau**4))
            value[:,i] = value[:,i] * fac
        
    else:

        value = r2d * 0.0

    return value

        
# ------------------------------------------------------------------------
# make the aurora
# ------------------------------------------------------------------------

def make_electron_aurora(mlts, lats, ae, ocflbBase, substorm, move):

    ssvalue = define_substorm_characteristics(mlts, lats, ae, ocflbBase, substorm, 0, move)
    ssvalue_avee = define_substorm_characteristics(mlts, lats, ae, ocflbBase, substorm, 1, move)
    
    amp_eflux = 3.0 + ae/150.0 * 0.75 # how brightness changes with AE
    amp_avee = 2.5 + ae/1000.0 # Aurora energy increases slightly

    # half e-folding widths:
    tau_eflux = 2.0 + ae/200.0 * 0.25 # how the thickness of the aurora grows
    tau_avee = tau_eflux * 2

    mltsR = mlts * np.pi / 12.0
    r = 90.0 - lats
    t2d, r2d = np.meshgrid(mltsR, r)
    eflux = amp_eflux * (np.cos(t2d)+2.0)/3.0
    avee = amp_avee * (np.cos(t2d)+5.0)/6.0

    tau_eflux_mlts = tau_eflux * (0.5 + 1.5 * (np.cos(mltsR)+1.0))/3.5

    # shift equatorward
    shift_mlts = -tau_eflux * 0.75
    nLats = len(lats)
    tau = np.zeros(nLats)
    for i, ocflb in enumerate(ocflbBase):
        tau_ef = tau_eflux_mlts[i]
        shift = shift_mlts
        dist = lats - (ocflb + shift)
        tau[dist >= 0.0] = tau_eflux_mlts[i] * 1.5
        tau[dist < 0.0] = tau_eflux_mlts[i] * 0.75
        fac = np.exp(-abs(dist**4/tau**4))
        eflux[:,i] = eflux[:,i] * fac

        fac = np.exp(-abs(dist**2/tau_avee**2))
        avee[:,i] = avee[:,i] * fac

    eflux = eflux + ssvalue * ae/25.0

    ssavee = 2.5 + ae/1000.0 * 4.0 # Aurora energy increases slightly
    aveess = ssvalue_avee * ssavee
    avee[aveess > avee] = aveess[aveess > avee]
    
    return eflux, avee
    
# ------------------------------------------------------------------------
# make the ION aurora
# ------------------------------------------------------------------------

def make_ion_aurora(mlts, lats, ae, ocflbBase, dst):

    dstp = -dst
    if (dstp < 0):
        dstp = 0.0
    
    # eflux is 0 mW/m2 + another 10 for every 200 nT change in Dst?
    amp_eflux = 0.0 + dstp/100.0 * 2.0 # how brightness changes with Dst

    # average energy is 40 keV + another 5 keV for every 200 nT change in Dst?
    amp_avee = 30.0 + dstp/200.0 * 5.0 # Aurora energy increases slightly

    # half e-folding width for the electron aurora
    tau_electron_eflux = 2.0 + ae/200.0 * 0.25 # how the thickness of the aurora grows

    # half e-folding width for the ion aurora
    tau_eflux = 3.0 + dstp/200.0 # how the thickness of the aurora grows
    tau_avee = tau_eflux * 2.0
    
    mltsR = mlts* np.pi / 12.0
    r = 90.0 - lats
    t2d, r2d = np.meshgrid(mltsR, r)

    offset = -3.0 * np.pi/12.0

    eflux = amp_eflux * (np.cos(t2d-offset)+2.0)/3.0
    avee = amp_avee * (np.cos(t2d-offset)+5.0)/6.0

    # this is the baseline electron stuff, so don't change
    tau_electron_eflux_mlts = tau_electron_eflux * (0.5 + 1.5 * (np.cos(mltsR)+1.0))/3.5
    tau_ion_eflux_mlts = tau_eflux * (0.5 + 1.5 * (np.cos(mltsR-offset)+1.0))/3.5

    # shift equatorward
    shift_mlts = -tau_electron_eflux_mlts * 1.5 - tau_ion_eflux_mlts * 0.75
    nLats = len(lats)
    tau = np.zeros(nLats)
    for i, ocflb in enumerate(ocflbBase):
        tau_ef = tau_ion_eflux_mlts[i]
        shift = shift_mlts[i]
        dist = lats - (ocflb + shift)
        # poleward of the center:
        tau[dist >= 0.0] = tau_ion_eflux_mlts[i] * 1.0
        # equatorward of the center:
        tau[dist < 0.0] = tau_ion_eflux_mlts[i] * 1.0
        fac = np.exp(-abs(dist**4/tau**4))
        eflux[:,i] = eflux[:,i] * fac

        fac = np.exp(-abs(dist**2/tau_avee**2))
        avee[:,i] = avee[:,i] * fac

    return eflux, avee
    
# ------------------------------------------------------------------------
# make potential
# ------------------------------------------------------------------------

def make_potential(mlts, lats, by, bz, ae, ocflbBase, substorm, move):

    ssvalue = define_substorm_characteristics(mlts, lats, ae, ocflbBase, substorm, 1, move)

    ocflb0 = np.mean(ocflbBase)
    
    # Do everything is kV, then transform to V at end
    amp_bz = -7.5
    amp_by = -5.0
    amp_ae = -30.0 * ae / 500.0 
    by_sharpen = (25.0 - np.abs(by)) / 25.0
    tau_potential_base = 7.0 * np.sqrt(by_sharpen)

    mltsR = mlts * np.pi / 12.0

    nMlts = len(mlts)
    nLats = len(lats)
    
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

    if (bz < 0.0):
        potential_bz = bz * amp_bz * np.sin(t2d)
    else:
        potential_bz = 0.0 * t2d
    potential_visc = -1.0 * amp_bz * np.sin(t2d)

    shifted = t2d - 23.0 * np.pi / 12.0
    potential_harang = 2.0 * amp_bz * ((np.cos(shifted)+1.0)/2.0)**3.0
    
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

        fac = np.exp(-2.0 * (dist + tau_potential_base*0.75)**2 / (tau*1.25)**2)
        potential_harang[:,i] = potential_harang[:,i] * fac
        

    potential_ss = ssvalue * amp_ae
    
    potential = (potential_by + \
                 potential_bz + \
                 potential_ss + \
                 potential_harang + \
                 potential_visc) * 1000.0

    return potential


# ------------------------------------------------------------------------
# This detects substorms based on the Newell and Gjerloev [2011] criteria
# A substorm onset occurs at T0, when:
#  al[T0+1] - al[T0] < -15
#  al[T0+2] - al[T0] < -30
#  al[T0+3] - al[T0] < -45
#  sum(al(T0 + i=4 to 30))/26 - al[T0] < -100
# ------------------------------------------------------------------------

def detect_substorms(al):

    nPts = len(al)
    substorm = np.zeros(nPts) - 1000.0

    IsSubstorm = 0
    HasPeaked = 0
    IsDone = 0
    iTStart = -1
    iTPeak = -1
    iT0 = 0
    while (iT0 < nPts-30):
        al_base = al[iT0]
        ind = np.arange(iT0+4, iT0+31)
        al_mean = np.mean(al[ind])
        if ((al[iT0 + 1] - al_base < -15.0) and
            (al[iT0 + 2] - al_base < -30.0) and
            (al[iT0 + 3] - al_base < -45.0) and
            (al_mean - al_base < -100.0)):
            # This means that there is a substorm
            IsSubstorm = 1
            HasPeaked = 0
            IsDone = 0
            iTStart = iT0
            iTPeak = -1
            iTEnd = -1

            iEnd = iT0 + 60
            if (iEnd > nPts):
                iEnd = nPts
            ind = np.arange(iT0, iEnd)
            mini = np.min(al[ind])
            for i in ind:
                if (al[i] == mini):
                    iTPeak = i

            iEnd = iTPeak + 120
            if (iEnd > nPts):
                iEnd = nPts
            ind = np.arange(iTPeak, iEnd)
            maxi = np.max(al[ind])
            for i in ind:
                if (al[i] == maxi):
                    iTEnd = i
            
            if ((iTPeak - iT0 > 20) and (iTEnd - iTPeak > 20)):

                print("Substorm Start @ ", iT0)
                print("         Peak  @ ", iTPeak)
                print("         End   @ ", iTEnd)
                
                di = (iTPeak - iT0)
                # this is the percentage of the way from the start to the peak:
                # put this as a negative to indicate expansion phase:
                for i in np.arange(iT0, iTPeak+1):
                    substorm[i] = -(i-iT0) / di

                # this is the percentage of the way from the peak to the end:
                di = (iTEnd - iTPeak)
                for i in np.arange(iTPeak, iTEnd+1):
                    substorm[i] = 1.0 - (i-iTPeak) / di

            iT0 = iTEnd

        iT0 = iT0 + 1

    return substorm

# ------------------------------------------------------------------------      
# smooth function                                                               
# ------------------------------------------------------------------------      

def smooth(times, values, dt):

    values = np.array(values)
    smoothed = np.array(values)

    iL = 0
    iH = 0
    n = len(times)
    for i, t in enumerate(times):
        tlow = t - dt
        thigh = t + dt
        while ((times[iL] < tlow) & (iH < n-1)):
            iL = iL+1
        while ((times[iH] < thigh) & (iH < n-1)):
            iH = iH+1
        iM = np.arange(iL,iH)
        v = values[iM]
        smoothed[i] = np.mean(v)

    return smoothed

# ------------------------------------------------------------------------
# main code
# ------------------------------------------------------------------------

# Get the input arguments
args = get_command_line_args(process_command_line_input())

move = args["move"]

if (args["real"]):

    results = download_omni_data(args["startdate"], args["enddate"], "-all")
    omniDirty = parse_omni_data(results)
    omni = clean_omni(omniDirty)
    basetime = dt.datetime.strptime(args["startdate"],"%Y%m%d")
    endtime = dt.datetime.strptime(args["enddate"],"%Y%m%d")

    xp = []
    for t in omni["times"]:
        xp.append((t-basetime).total_seconds())

    fp = omni["bz"]
        
    dt_in_sec = args["dt"]*60.0
    totaltime = (endtime - basetime).total_seconds()
    alltimes = np.arange(0.0, totaltime+dt_in_sec, dt_in_sec)

    dt = 15.0 * 60.0
    omni_smoothed_bz = smooth(xp, omni["bz"], dt)
    
    bx = np.interp(alltimes, xp, omni["bx"])
    by = np.interp(alltimes, xp, omni["by"])
    bz = np.interp(alltimes, xp, omni["bz"])
    vx = np.abs(np.interp(alltimes, xp, omni["vx"]))

    bzs = np.interp(alltimes, xp, omni_smoothed_bz)    
    
    ae = np.interp(alltimes, xp, omni["ae"])
    al = np.interp(alltimes, xp, omni["al"])
    au = np.interp(alltimes, xp, omni["au"])

    # --------------------------------------
    # If we have to get Dst, do the same thing:
    # --------------------------------------

    ions = args["ions"]
    if (ions):
        # This blindly assumes that the kyoto Dst exists....
        dstfile = download_kyoto_dst(args["startdate"])
        kyoto = parse_kyoto_dst(dstfile, args["startdate"])

        xp = []
        for t in kyoto["times"]:
            xp.append((t-basetime).total_seconds())
        dst = np.interp(alltimes, xp, kyoto["dst"])
    
    if (args["substorm"]):
        dt_omni = (omni["times"][1] - omni["times"][0]).total_seconds()
        if (dt_omni != 60.0):
            print("In order to run with a substorm, dt for omni needs to be 60s")
            exit()
        substorm_omni = detect_substorms(np.array(omni["al"]))
        substorm = np.interp(alltimes, xp, substorm_omni)

    else:
        substorm = np.zeros(len(al)) - 1000.0
            
else:

    dt = 1.0/60.0
    subtimes = np.arange(12,13+dt,dt)
    nTimes = len(subtimes)
    bysub = np.random.normal(-5.0,2.0,nTimes)
    bzsub = bysub * 0.0 - 2.0

    bx = np.zeros(nTimes)
    by = np.concatenate([[-5.0], bysub,  [-5.0]])
    bz = np.concatenate([[-2.0], bzsub, [0.0]])
    alltimes = np.concatenate([[0.0], subtimes, [24.0]])

    basetime = datetime(2020, 3, 20, 0, 0, 0)
    alltimes = np.array(alltimes) * 3600.0
        
DoAddVarOCFLB = 0

time_1965 = datetime(1965, 1, 1, 0, 0, 0)

nTimes = len(bz)

dLat = 0.5
minLat = 50.0
nLats = (90.0 - minLat)/dLat + 1
dMlt = 1.0/3.0
nMlts = (24.0-0.0)/dMlt + 1

lats = np.arange(minLat, 90.0+dLat, dLat)
mlts = np.arange(0.0, 24.0+dMlt, dMlt)

theta2d, r2d = np.meshgrid(mlts * np.pi/12.0 - np.pi/2.0, 90.0 - lats)
area = 111.0 * 111.0 * 1000.0 * 1000.0 * np.sin(r2d*np.pi/180.0)

data = {}

data["nLats"] = nLats
data["nMlts"] = nMlts
data["nTimes"] = nTimes

data["lats"] = lats
data["mlts"] = mlts

data["times"] = []
#for t in alltimes:
#    data["times"].append(basetime + timedelta(seconds = t))

data["imf"] = []
data["ae"] = []
data["dst"] = []
data["hp"] = []
data["cpcp"] = []

hparray = []
hpiarray = []
bzarray = []
byarray = []
aearray = []

data["Vars"] = ['Potential (kV)', \
                'Electron Energy Flux (ergs/cm2/s)', \
                'Electron Mean Energy (keV)']
if (ions):
    data["Vars"].append('Ion Energy Flux (ergs/cm2/s)')
    data["Vars"].append('Ion Mean Energy (keV)')

data["nVars"] = len(data["Vars"])

data["version"] = 1.3

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

    byNow = by[i]
    if (args["south"]):
        byNow = -byNow

    aeCurrent = ae[i]
    aesat = 1000.0 * (1-np.exp(-aeCurrent/500.0))
    if (aesat < aeCurrent):
        aeCurrent = aesat
    
    ocflb = make_ocflb(mlts, byNow, bzs[i], DoAddVarOCFLB, phase, substorm[i], move)
    pot2d = make_potential(mlts, lats, byNow, bz[i], aeCurrent, ocflb, substorm[i], move)
    eflux2d, avee2d = make_electron_aurora(mlts, lats, aeCurrent, ocflb, substorm[i], move)
    power = eflux2d/1000.0 * area # In Watts
    hp = np.sum(power)/1.0e9 # In GW

    if (ions):
        ionEflux2d, ionAvee2d = make_ion_aurora(mlts, lats, aeCurrent, ocflb, dst[i])
        ionPower = ionEflux2d/1000.0 * area
        ionHp = np.sum(ionPower)/1.0e9 # in GW
    else:
        ionHp = hp * 0.0

    # print(ut, int(by[i]), int(bz[i]), int(ae[i]), int(hp))
        
    #  ; IMF should be (nTimes,4) - V, Bx, By, Bz
    #  ; AE  (nTimes, 4) - AL, AU, AE, AEI?
    #  ; Dst (nTimes, 2) - Dst, Dsti?
    #  ; Hpi (nTimes, 2) - HP, Joule Heating (GW)
    #  ; CPCP (nTimes) - CPCP (kV)

    data["imf"].append([vx[i], bx[i], by[i], bz[i]])
    bzarray.append(bz[i])
    byarray.append(by[i])
    # These are all bad values for now....
    data["ae"].append([al[i], au[i], ae[i], 0.0])
    aearray.append(ae[i])
    data["dst"].append([0.0, 0.0])
    data["hp"].append([hp, ionHp])
    hparray.append(hp)
    hpiarray.append(ionHp)
    data["cpcp"].append((np.max(pot2d) - np.min(pot2d))/1000.0)

    cPot = data["Vars"][0]
    data[cPot].append(pot2d)
    
    cEFlux = data["Vars"][1]
    data[cEFlux].append(eflux2d)
    
    cAveE = data["Vars"][2]
    data[cAveE].append(avee2d)

    if (ions):
        cIonEFlux = data["Vars"][3]
        data[cIonEFlux].append(ionEflux2d)
        cIonAveE = data["Vars"][4]
        data[cIonAveE].append(ionAvee2d)

if (args["outfile"].find(".nc") > 0):
    # Write out in new netCDF format:
    amie_write_netcdf(args["outfile"], data)
else:
    # Write out in old AMIE format:
    amie_write_binary(args["outfile"], data)


fig = plt.figure(figsize = (10,10))
zeros = np.array(bz) * 0.0

ax = fig.add_subplot(511)
ax.plot(data["times"], bzarray)
ax.plot(data["times"], zeros, 'k:')
ax.set_ylabel('IMF Bz (nT)')
ax.set_xlim(data["times"][0],data["times"][-1])

ax = fig.add_subplot(512)
ax.plot(data["times"], byarray)
ax.plot(data["times"], zeros, 'k:')
ax.set_ylabel('IMF By (nT)')
ax.set_xlim(data["times"][0],data["times"][-1])

ax = fig.add_subplot(514)
ax.plot(data["times"], hparray)
if (ions):
    ax.plot(data["times"], hpiarray)
#ax.plot(times, zeros, 'k:')
ax.set_ylabel('Hemispheric Power (GW)')
ax.set_xlim(data["times"][0],data["times"][-1])

ax = fig.add_subplot(513)
ax.plot(data["times"], aearray)
#ax.plot(times, zeros, 'k:')
ax.set_ylabel('AE (nT)')
ax.set_xlim(data["times"][0],data["times"][-1])

ax = fig.add_subplot(515)
ax.plot(data["times"], data["cpcp"])
#ax.plot(times, zeros, 'k:')
ax.set_ylabel('CPCP (kV)')
ax.set_xlim(data["times"][0],data["times"][-1])


plotfile = 'electrodynamics.png'
fig.savefig(plotfile)
plt.close()
