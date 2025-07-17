#!/usr/bin/env python

# Authors of this code:
# Daniel A. Brandt, Ph.D., Michigan Tech Research Institute, daabrand@mtu.edu

# This file contains a suite of tools that do the following:
# 1 - Obtain F10.7 data between any dates of the user's choosing.
# 2 - Generate NEUVAC irradiances between any two dates of the user's choosing.
# 3 - Output the NEUVAC irradiances to a .csv file to be used by Aether, either in the b37 or b59 bins.

# Top-level imports:
import argparse
import numpy as np
from datetime import datetime
import pathlib
from pathlib import Path
import pandas as pd
from scipy.interpolate import CubicSpline
import urllib.request, pickle
from fism import saveFism
import os, sys
import pooch

# Directory management:
here = pathlib.Path(__file__).parent.resolve()
euvDir = here.parent.joinpath('share/run/UA/inputs')

# Physical constants:
h = 6.62607015e-34 # Planck's constant in SI units of J s
c = 299792458 # Speed of light in m s^-1

# Global variable(s):
# Waves Table (coefficients for the old NEUVAC Model):
# Format: 0-Min, 1-Max, 2-S_1i, 3-S_Ai, 4-S_Di, 5-I_i, 6-Pi, 7-Ai
waveTable = np.array([
    [1700.00, 1750.00, 1.31491e-06, 6.71054e-06, 5.78034e-07, 0.00355128, 1.05517, 0.901612],
    [1650.00, 1700.00, 5.19285e-07, 2.62376e-06, 3.08447e-07, 0.00218156, 1.06245, 0.964892],
    [1600.00, 1650.00, 3.85348e-07, 1.73851e-06, 3.34911e-07, 0.00115310, 1.07246, 0.959562],
    [1550.00, 1600.00, 2.96220e-07, 1.29250e-06, 2.61812e-07, 0.000814814, 1.04567, 0.967804],
    [1500.00, 1550.00, 2.35326e-07, 1.21123e-06, 2.27793e-07, 0.000566574, 1.13520, 0.970257],
    [1450.00, 1500.00, 1.86793e-07, 5.96399e-07, 1.48283e-07, 0.000331058, 1.01564, 0.940506],
    [1400.00, 1450.00, 1.96396e-07, 5.84154e-07, 1.82438e-07, 0.000207013, 1.67546, 0.945697],
    [1350.00, 1400.00, 1.04362e-07, 5.02422e-07, 1.45100e-07, 0.000153277, 1.04246, 0.992749],
    [1300.00, 1350.00, 1.74403e-07, 6.32214e-07, 4.03009e-07, 0.000311075, 1.00964, 1.09381],
    [1250.00, 1300.00, 7.12738e-08, 2.44220e-07, 9.56532e-08, 9.68823e-05, 1.15737, 1.01121],
    [1200.00, 1250.00, 8.74335e-06, 5.02272e-05, 1.32536e-05, 0.00263307, 1.46273, 0.987493],
    [1215.67, 1215.67, 6.43713e-06, 5.16823e-05, 1.11399e-05, 0.00247063, 1.26340, 0.998295],
    [1150.00, 1200.00, 1.15468e-07, 2.74916e-07, 1.65125e-07, 0.000105178, 1.66887, 1.00997],
    [1100.00, 1150.00, 7.71861e-08, 2.15061e-07, 1.44227e-07, 5.16157e-05, 0.971988, 1.05634],
    [1050.00, 1100.00, 5.84127e-08, 3.08808e-07, 1.25160e-07, 4.65227e-05, 1.58808, 1.05327],
    [1000.00, 1050.00, 2.23073e-07, 6.92710e-07, 5.19444e-07, 5.44992e-05, 0.449052, 1.10271],
    [1031.91, 1031.91, 6.18723e-08, 1.21679e-07, 2.28527e-07, 3.14905e-05, 1.42684, 1.17863],
    [1025.72, 1025.72, 1.61504e-07, 4.38856e-07, 2.79663e-07, 1.06365e-05, 1.09262, 1.05186],
    [950.00, 1000.00, 1.70358e-07, 5.20531e-07, 3.86006e-07, 3.34989e-05, 0.491283, 1.09676],
    [977.02, 977.02, 1.51857e-07, 5.60743e-07, 2.74541e-07, 6.71100e-06, 1.44918, 1.04869],
    [900.00, 950.00, 7.27646e-08, 4.53511e-07, 1.91513e-07, 3.93851e-05, 1.21476, 1.06473],
    [850.00, 900.00, 1.45264e-07, 2.82927e-07, 4.22856e-07, 4.83494e-05, 1.15579, 1.14948],
    [800.00, 850.00, 6.69560e-08, 1.26613e-07, 1.76066e-07, 3.69687e-05, 1.14722, 1.12832],
    [750.00, 800.00, 3.22816e-08, 7.81757e-08, 6.32959e-08, 4.42679e-05, 0.969748, 1.06692],
    [789.36, 789.36, 1.19733e-08, 2.53334e-08, 1.58546e-08, 1.25539e-05, 1.48302, 1.00982],
    [770.41, 770.41, 7.33597e-09, 2.10650e-08, 1.63125e-08, 8.88041e-06, 1.18634, 1.06584],
    [765.15, 765.15, 4.85967e-09, 1.05567e-08, 5.42104e-09, 1.15262e-05, 1.17912, 1.03352],
    [700.00, 750.00, 1.85139e-08, 3.63837e-08, 3.29576e-08, 1.72134e-05, 1.25328, 1.06364],
    [703.36, 703.36, 5.34708e-09, 9.65120e-09, 4.54419e-09, 8.80278e-06, 1.51207, 0.972520],
    [650.00, 700.00, 1.79851e-08, 6.39605e-08, 1.86000e-08, 1.41950e-05, 1.11181, 0.945801],
    [600.00, 650.00, 1.52595e-07, 5.29641e-07, 1.41837e-07, 3.96165e-05, 1.00554, 0.949913],
    [629.73, 629.73, 4.96048e-08, 2.46454e-07, 3.12902e-08, 1.59200e-05, 1.01611, 0.846628],
    [609.76, 609.76, 2.80641e-08, 3.24530e-07, 1.81554e-08, 1.68460e-06, 0.973085, 0.793355],
    [550.00, 600.00, 1.12234e-07, 6.29889e-07, 1.56092e-07, 2.79143e-05, 0.961457, 0.970150],
    [584.33, 584.33, 7.91646e-08, 3.05430e-07, 5.14430e-08, 1.70372e-05, 0.844250, 0.881026],
    [554.31, 554.31, 2.47485e-08, 2.68042e-07, 5.40951e-08, 1.16226e-06, 1.08699, 1.01483],
    [500.00, 550.00, 1.12037e-07, 7.84515e-07, 6.32364e-08, 4.55230e-06, 1.13480, 0.816868],
    [450.00, 500.00, 1.10016e-07, 3.96192e-07, 7.37101e-08, 2.62692e-05, 1.15344, 0.865234],
    [465.22, 465.22, 9.60010e-09, 1.75358e-08, 6.91440e-11, 1.45142e-05, 1.62256, -0.203971],
    [400.00, 450.00, 5.15555e-08, 2.89821e-07, 3.85807e-08, 1.64207e-05, 1.36652, 0.893190],
    [350.00, 400.00, 3.91955e-07, 1.43942e-06, 3.16713e-07, -2.36108e-06, 1.05819, 0.910235],
    [368.07, 368.07, 1.38855e-07, 7.21254e-07, 1.01814e-07, 8.71098e-07, 1.26707, 0.890513],
    [300.00, 350.00, 1.35439e-06, 1.09238e-05, 8.24308e-07, 4.35250e-05, 1.22619, 0.816515],
    [303.78, 303.78, 7.43959e-07, 5.94012e-06, 4.05188e-07, 9.23799e-05, 1.32976, 0.796970],
    [303.31, 303.31, 5.25977e-07, 7.87164e-06, 3.07932e-07, 7.87468e-05, 0.945961, 0.759694],
    [250.00, 300.00, 9.10710e-07, 3.91586e-06, 1.20177e-06, -9.64301e-06, 1.07360, 0.958369],
    [284.15, 284.15, 8.67633e-07, 6.00671e-06, 3.97664e-07, -0.000107230, 1.20608, 0.773950],
    [256.30, 256.30, 6.44996e-08, 4.12637e-07, 1.05193e-07, 6.61853e-06, 1.48670, 1.03265],
    [200.00, 250.00, 4.83013e-07, 1.18898e-06, 8.94772e-07, 5.34779e-05, 1.04532, 1.07888],
    [150.00, 200.00, 7.13305e-07, 2.47623e-06, 9.78936e-07, 0.000261230, 1.47374, 1.01156],
    [100.00, 150.00, 4.03676e-08, 2.28270e-07, 4.43965e-08, 2.16162e-05, 1.09062, 0.970310],
    [50.00, 100.00, 1.69769e-07, 6.93618e-07, 2.89457e-07, 2.03013e-05, 1.07887, 1.06022],
    [32.00, 50.00, 1.23478e-07, 4.43644e-07, 1.75749e-07, -1.34567e-05, 1.27409, 1.01254],
    [23.00, 32.00, 6.10174e-08, 2.34313e-07, 1.10591e-07, -1.22729e-05, 0.699812, 1.04841],
    [16.00, 23.00, 2.23866e-07, 7.97533e-07, 3.03563e-07, -5.62012e-05, 0.706360, 0.987835],
    [8.00, 16.00, 3.10773e-07, 1.22767e-06, 3.74797e-07, -8.41459e-05, 1.39529, 0.963859],
    [4.00, 8.00, 1.17378e-08, 7.13970e-08, 1.38839e-08, -3.63146e-06, 0.811119, 0.920702],
    [2.00, 4.00, 3.97985e-09, 4.12085e-08, 4.71914e-09, -1.86099e-06, 1.15214, 0.916686],
    [1.00, 2.00, 3.52498e-09, 1.57342e-08, 4.03741e-09, -8.84488e-07, 0.951714, 0.943490]
    ])

# Helper Functions:
def rollingAverage(myData, window_length=1, impute_edges=True, center=True):
    """
    Using pandas, compute a rolling average of over 'data' using a window length of 'windowlength'. Sets the leading and
    trailing windows to the values of the original data.
    :param myData: arraylike
        The data over which to compute the rolling average.
    :param window_length: int
        The size of the window over which to average.
    :param impute_edges: bool
        A boolean determining whether the edges will be interpolated. Default is True.
    :param center: bool
        A boolean determining whether the centered average will be used.
    :return: rolled, arraylike
        The rolling average data.
    """
    myDataframe = pd.DataFrame(data=myData, columns=['Var'])
    myDataframe['Rolling'] = myDataframe['Var'].rolling(window=window_length, center=center).mean()
    firstValidIndex = myDataframe['Rolling'].first_valid_index()
    lastValidIndex = myDataframe['Rolling'].last_valid_index()
    if impute_edges == True:
        # Sample x-axis:
        sampleXaxis = np.linspace(0, window_length, window_length)
        middleIndex = int(0.5*window_length)
        # Use cubic interpolation to fill the gaps on the edges:
        leadingEdgeStartingVal = myDataframe['Var'][:window_length].values[0]
        leadingEndingVal = myDataframe['Rolling'][firstValidIndex]
        leadingEdgeMiddleVal = np.mean([leadingEdgeStartingVal, leadingEndingVal])
        leadingSpline = CubicSpline([sampleXaxis[0], sampleXaxis[middleIndex], sampleXaxis[-1]],
                                    [leadingEdgeStartingVal, leadingEdgeMiddleVal, leadingEndingVal])
        leadingImputedValues = leadingSpline(sampleXaxis)

        trailingEdgeStartingVal = myDataframe['Rolling'][lastValidIndex]
        trailingEndingVal = myDataframe['Var'].values[-1]
        trailingEdgeMiddleVal = np.mean([trailingEdgeStartingVal, trailingEndingVal])
        trailingSpline = CubicSpline([sampleXaxis[0], sampleXaxis[middleIndex], sampleXaxis[-1]],
                                    [trailingEdgeStartingVal, trailingEdgeMiddleVal, trailingEndingVal])
        trailingImputedValues = trailingSpline(sampleXaxis)
        # Ingest the imputed values:
        myDataframe['Rolling'][:window_length] = leadingImputedValues
        myDataframe['Rolling'][-window_length:] = trailingImputedValues
    else:
        myDataframe['Rolling'][:window_length] = myDataframe['Var'][:window_length]
        myDataframe['Rolling'][-window_length:] = myDataframe['Var'][-window_length:]
    rolled = myDataframe['Rolling'].values
    return rolled

def readCLS(filename):
    """
    Load in flare-corrected, Sun-Earth distance adjusted flux values recorded by the Collecte Localisation Satellites
    (CLS).
    :param filename: str
        The location of the data file.
    :return times: list
        The datetimes for each data value.
    :return data: ndarray
        The solar flux data for F30, F15, F10.7, F8, and F3.2.
    """
    times = []
    precisionVals = []
    with open(filename, 'r') as myFile:
        allLines = myFile.readlines()
        data = np.zeros((len(allLines)-25, 5))
        i = 0
        j = 0
        for line in allLines:
            if i >= 25:
                elements = line.split()
                data[j, :] = np.array([float(elements[5]), float(elements[9]), float(elements[13]), float(elements[17]), float(elements[21])])
                times.append( datetime(int(elements[0]), int(elements[1]), int(elements[2]), 12) )
                precisionVals.append( [float(elements[6]), float(elements[10]), float(elements[14]), float(elements[18]), float(elements[22])] )
                j += 1
            i += 1
    # Print the precision:
    # print('Mean precision values...')
    # print('F30: '+str(np.nanmean([element[0] for element in precisionVals]))+' sfu') # 6
    # print('F15: ' + str(np.nanmean([element[1] for element in precisionVals])) + ' sfu') # 8
    # print('F10.7: ' + str(np.nanmean([element[2] for element in precisionVals])) + ' sfu') # 13
    # print('F8: ' + str(np.nanmean([element[3] for element in precisionVals])) + ' sfu') # 12
    # print('F3.2: ' + str(np.nanmean([element[4] for element in precisionVals])) + ' sfu') # 11
    return times, data

def getCLSF107(dateStart, dateEnd, truncate=True):
    """
    Obtains flare-corrected F10.7 data from Collecte Localisation Satellites. The "adjusted" here means that it has been
    adjusted from measurements from Earth to 1AU. (Aether/GITM need measurements at 1AU, so they can adjust to the 
    proper sun-planet distance, where planet can be Earth, Venus, Mars, etc.) A description of the data is
    provided here: https://spaceweather.cls.fr/services/radioflux/. 
    Downloads the most recent measurements to a file. Reads the file and extracts the F10.7 values between two dates. 
    Note that if the ending date is less than or equal to the last date in the version of the file that has already 
    been downloaded, the file IS NOT re-downloaded, but simply parsed. Otherwise, the file is redownloaded.
    :param dateStart: str
        The starting date in YYYYMMDD format.
    :param dateEnd: str
        The ending date in YYYYMMDD format.
    :param truncate: bool
        Controls whether to truncate the data to exclude the most recent 81 days. Defaults is True.
    """
    dateStart = dateStart[:4]+'-'+dateStart[4:6]+'-'+dateStart[6:]
    dateEnd = dateEnd[:4]+'-'+dateEnd[4:6]+'-'+dateEnd[6:]
    dateTimeStart = datetime.strptime(dateStart, '%Y-%m-%d')
    dateTimeEnd = datetime.strptime(dateEnd, '%Y-%m-%d')
    fname = euvDir.joinpath("radio_flux_adjusted_observation.txt")
    if fname.exists():
        # Read in the file:
        times, data = readCLS(fname)
        # Check if the ending date exceeds the ending date in the file. If so, redownloading the file:
        if times[-1] > dateTimeEnd:
            out = urllib.request.urlretrieve(
                'ftp://ftpsedr.cls.fr/pub/previsol/solarflux/observation/radio_flux_adjusted_observation.txt', fname)
        times, data = readCLS(fname)
    else:
        # Download the file:
        out = urllib.request.urlretrieve('ftp://ftpsedr.cls.fr/pub/previsol/solarflux/observation/radio_flux_adjusted_observation.txt', fname)
        times, data = readCLS(fname)

    # Compute the 81-day (centered) averaged F10.7 and 54-day averaged (:
    F107 = data[:, 2]
    F107A = rollingAverage(F107, window_length=81, impute_edges=True)
    F107B = rollingAverage(F107, window_length=54, impute_edges=True, center=False)
    print(f'\n\nlen(F107) = {len(F107)}, len(F107B) = {len(F107B)}, len(F107B) = {len(F107B)}')
    # Extract the values in the desired time range:
    goodInds = np.where((np.asarray(times) >= dateTimeStart) & (np.asarray(times) <= dateTimeEnd))[0]
    # Truncation:
    if truncate and len(goodInds) >= 2*81:
        goodInds = goodInds[:-81]
    return np.asarray(times)[goodInds], np.asarray(F107)[goodInds], np.asarray(F107A)[goodInds], np.asarray(F107B)[goodInds]

def mycorrelate2d(df, normalized=False):
    """
    Compute the correlation matrix from 2D data, where each row is cross correlated with the others.
    This function handles NaN values by ignoring them.
    :param df: ndarray
        A 2D array of dimensions n x m.
    :param normalized: bool
        Determines whether the resulting correlation matrix is normalized. Default is False.
    :returns ccm: ndarray
        The [normalized] cross-correlation matrix.
    Source: https://stackoverflow.com/questions/54292947/basics-of-normalizing-cross-correlation-with-a-view-to-comparing-signals
    """
    # Initialize cross correlation matrix with zeros
    ccm = np.zeros((df.shape[1], df.shape[1]))
    # Fill in each entry of the matrix one-by-one:
    for i in range(df.shape[1]):
        outer_row = df[:, i]
        for j in range(df.shape[1]):
            inner_row = df[:, j]
            goodInds = np.logical_and(~np.isnan(outer_row), ~np.isnan(inner_row))
            if (not normalized):
                x = np.correlate(outer_row[goodInds], inner_row[goodInds])
            else:
                x = get_cc(outer_row[goodInds], inner_row[goodInds])
                # a = (inner_row - np.mean(inner_row)) / (np.std(inner_row) * len(inner_row))
                # b = (outer_row - np.mean(outer_row)) / (np.std(outer_row) )
                # x = np.correlate(a, b)
            ccm[i, j] = x
    return ccm

def get_cc(array1, array2, normalize=True):
    """
    Compute the cross-correlation of two 1D arrays of the same length.
    :param array1: ndarray
        A 1D array of length n.
    :param array2: ndarray
        A 1D array of length n.
    :return c: float
        The normalized correlation of the two arrays.
    """
    if normalize:
        a = (array1 - np.mean(array1)) / (np.std(array1) * len(array1))
        b = (array2 - np.mean(array2)) / (np.std(array2))
        c = np.correlate(a, b)
    else:
        c = np.correlate(array1, array2)
    return c

def loadPickle(pickleFilename):
    """
    Given the name of a (pre-existing) pickle file, load its contents.
    :param: pickleFilename, str
        A string with the location/name of the filename.
    :return: var
        The loaded data.
    """
    with open(pickleFilename, 'rb') as pickleFile:
        var = pickle.load(pickleFile)
    return var

# CORE NEUVAC FUNCTIONS:
def irrFunc(F107input, A, B, C, D, E, F):
    F107, F107A = F107input
    return A * (F107 ** B) + C * (F107A ** D) + E * (F107A - F107) + F

def neuvacEUV(f107, f107b, bands=None, tableFile=None, statsFiles=None):
    """
    Use a parametric model to compute solar flux in the 59 conventional wavelength bands used by Aether/GITM. Capable
    of returning perturbed irradiance values that are perturbed according to the variations of in the intensity of each
    bin
    :param f107: ndarray
        F10.7 values.
    :param f107b: ndarray
        81-day center-averaged F10.7 values; must be the same length as f107.
    :param bands: str
        If None or 'NEUVAC', returns irradiances in the GITM Bands. If 'EUVAC', returns them in the 37 bands used by
        EUVAC. If 'SOLOMON', returns them in the 22 bands used by Solomon and Qian.
    :param tableFile: str
        Corresponds to the .txt file holding the NEUVAC coefficients most recently-generated by fitNeuvac.py. If
        not given, simply uses the table file corresponding to the selected bin structure. Default is None.
    :param statsFiles: Bool
        Determines whether data for uncertainty quantification is exploited. Involves usage of a list containing
        2 elements where the first element is a file containing the 59x59 correlation matrix and the second
        element is a file containing the 1x59 standard deviation values for NEUVAC. NOT REQUIRED.
    :return euvIrradiance: ndarray
        A nxm ndarray where n is the number of EUV irradiance values and m is the number of wavelength bands.
    :return perturbedEuvIrradiance: ndarray
        A nxm ndarray where n is the number of EUV irradiance values perturbed due to inherent uncertainty and m is the
        number of wavelength bands.
    :return savedPerts: ndarray
        A nxm ndarray of the perturbations (time series of the NEUVAC+Perturbation - NEUVAC)
    :return cc2: ndarray
        A mxm ndarray of the correlation matrix between each wavelength's time-series of the NEUVAC+Perturbation -
        NEUVAC.
    """
    if type(f107) != np.ndarray:
        f107 = np.asarray([f107])
        f107b = np.asarray([f107b])
        if bands == 'SOLOMON':
            solarFlux = np.zeros((1, 22))
        else:
            solarFlux = np.zeros((1, waveTable.shape[0]))
    else:
        if bands == 'SOLOMON':
            solarFlux = np.zeros((len(f107), 22))
        else:
            solarFlux = np.zeros((len(f107), waveTable.shape[0]))
    euvIrradiance = np.zeros_like(solarFlux)
    perturbedEuvIrradiance = np.zeros_like(solarFlux)
    # Gather the model parameters:
    if tableFile is None:
        if bands == 'SOLOMON':
            tableFile = euvDir.joinpath('neuvac_table_stan_bands.txt') #'../data/neuvac_table_stan_bands.txt'
        else:
            tableFile = euvDir.joinpath('neuvac_table.txt') #'../data/neuvac_table.txt'
    neuvacTable = []
    with open(here.parent.joinpath(tableFile)) as neuvacFile: # open(tableFile)
        contents = neuvacFile.readlines()
        i = 0
        for line in contents:
            if i > 17:
                neuvacTable.append([float(element) for element in line.split(' ')])
            i+=1
    neuvacTable = np.asarray(neuvacTable)

    # If no stats file is provided, simply return the base model output (using the required table file):
    if not statsFiles:
        # Loop across the F10.7 (and F10.7A) values:
        for i in range(len(f107)):
            k = 0
            for j in (range(solarFlux.shape[1])):
                irrRes = irrFunc([f107[i], f107b[i]], *neuvacTable[j, 2:])
                if irrRes < 0:
                    irrRes = 0
                    euvIrradiance[i, k] = irrRes
                else:
                    euvIrradiance[i, k] = irrRes
                k += 1
        if bands == 'EUVAC':  # Returns values ONLY for those corresponding to the wavelengths used by EUVAC
            return euvIrradiance[:, 7:44], None, None, None
        else:
            return euvIrradiance, None, None, None
    else:
        # Include statistical data for calculating uncertainties via perturbations:
        if bands == 'NEUVAC':
            statsFiles = [euvDir.joinpath('corMat.pkl'), euvDir.joinpath('sigma_NEUVAC.pkl')]
        elif bands == 'EUVAC':
            statsFiles = [euvDir.joinpath('corMatEUVAC.pkl'), euvDir.joinpath('sigma_EUVAC.pkl')]
        else:
            statsFiles = [euvDir.joinpath('corMatStanBands.pkl'), euvDir.joinpath('sigma_NEUVAC_StanBands.pkl')]
        corMatFile = statsFiles[0]
        corMat = loadPickle(corMatFile)
        sigmaFile = statsFiles[1]
        STDNeuvacResids = loadPickle(sigmaFile)
        # Loop across the F10.7 (and F10.7A) values:
        nTimes = len(f107)
        nWaves = solarFlux.shape[1]
        savedPerts = np.zeros((nTimes, nWaves))
        for i in range(len(f107)):
            # Loop across the wavelengths (59 conventional wavelengths):
            k = 0
            P_n = []
            for j in (range(solarFlux.shape[1])):
                # Percentage perturbation:
                P_j = np.random.normal(0, 1.0)
                P_n.append(P_j)
                P_1 = P_n[0]
                # Normalized Correlated Perturbation:
                if bands == 'SOLOMON':
                    if j < 5:
                        C_j1 = corMat[0, j] # 3 # Only consider correlation with the third wavelength bin of the SOLOMON bins!
                    else:
                        C_j1 = corMat[5, j] # Only consider correlation with the fifth wavelength bin of the SOLOMON bins!
                else:
                    if j < 7:
                        # Only consider correlation with the third wavelength bin (of the NEUVAC bins!) when bands are below 8.
                        C_j1 = corMat[0, j] # 2
                    else:
                        # Only consider correlation with the first wavelength bin (of the EUVAC bins!) when bands are above 8.
                        C_j1 = corMat[7, j]
                N_j = C_j1 * P_1 + (1.0 - C_j1) * P_j
                # Actual Normalized Correlated Perturbation:
                A_j = STDNeuvacResids[j] * N_j
                irrRes = irrFunc([f107[i], f107b[i]], *neuvacTable[j, 2:])
                if irrRes < 0:
                    irrRes = 0
                    euvIrradiance[i, k] = irrRes
                    if irrRes + A_j < 0:
                        perturbedEuvIrradiance[i, k] = 0
                    else:
                        perturbedEuvIrradiance[i, k] = irrRes + A_j
                else:
                    euvIrradiance[i, k] = irrRes
                    perturbedEuvIrradiance[i, k] = irrRes + A_j
                savedPerts[i, j] = A_j
                k += 1

        # Generate a correlation matrix of the perturbations (to compare to the input correlation matrix as a sanity check):
        cc2 = mycorrelate2d(savedPerts, normalized=True)

        if bands == 'EUVAC':  # Returns values ONLY for those corresponding to the wavelengths used by EUVAC
            return euvIrradiance[:, 7:44], perturbedEuvIrradiance[:, 7:44], savedPerts, cc2
        else:
            return euvIrradiance, perturbedEuvIrradiance, savedPerts, cc2

# -----------------------------------------------------------------------------------------------------------------------------------------
# Argument Parsing Function:
def get_args():

    parser = argparse.ArgumentParser(description = 'Create NEUVAC input data')
    parser.add_argument('start',
                        help='Start date (format YYYYMMDD)',
                        type=str)
    parser.add_argument('end',
                        help='End date (format YYYYMMDD)',
                        type=str)
    parser.add_argument('-b', '--binning',
                        help="Binning scheme to use. Can be [solomon,neuvac,euvac] "
                        "(case insensitive)",
                        type=str, default="neuvac")
    args = parser.parse_args()

    return args
# -----------------------------------------------------------------------------------------------------------------------------------------

# Example Execution:

# python neuvac.py 20110319 20110321 -b euvac

args = get_args()
dateStart = args.start
dateEnd = args.end
binning_scheme = args.binning

# Load F10.7 data (from Collecte Localisation Satellites):
times, F107, F107A, F107B = getCLSF107(dateStart, dateEnd)

# Generate NEUVAC Irradiance from that F10.7 data:
if binning_scheme == 'HFG' or binning_scheme == 'SOLOMON' or binning_scheme == 'Solomon' or binning_scheme == 'solomon':
    # SOLOMON (STAN BANDS; b23)
    irradiance, _, _, _ = neuvacEUV(F107, F107B, bands='SOLOMON')
else:
    if binning_scheme == 'NEUVAC' or binning_scheme == 'Neuvac' or binning_scheme == 'neuvac':
        # NEUVAC BINS (b59)
        irradiance, _, _, _ = neuvacEUV(F107, F107B, bands='NEUVAC')
    else:
        # EUVAC BINS (b37)
        irradiance, _, _, _ = neuvacEUV(F107, F107B, bands='EUVAC')

# Save the NEUVAC Irradiance to a file that Aether can use:
tag = str(irradiance.shape[1])
fname = os.getcwd() + '/neuvac_file_'+tag+'.txt'
saveFism(irradiance, times, fname)

