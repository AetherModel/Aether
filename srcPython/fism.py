#!/usr/bin/env python

# Authors of this code:
# Daniel A. Brandt, Ph.D., Michigan Tech Research Institute, daabrand@mtu.edu
# Aaron J. Ridley, Ph.D., University of Michigan, ridley@umich.edu

# This file contains a suite of tools that do the following:
# 1. Download FISM2 data for a time period the user desires.
# 2. Rebin that data into the binning scheme the user desires (i.e. EUVAC-37, NEUVAC-59, or SOLOMON).
# 3. Outputs a FISM2 file with the rebinnined irradiances in the desired bins (for use by euv.cpp)

# Top-level imports:
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib
matplotlib.use('Qt5Agg')
import pathlib
import os
from urllib.request import urlretrieve
from netCDF4 import Dataset
import scipy.integrate as integ

# Directory management:
here = pathlib.Path(__file__).parent.resolve()

# Physical constants:
h = 6.62607015e-34 # Planck's constant in SI units of J s
c = 299792458 # Speed of light in m s^-1

# Helper Functions:
def getFism(dateStart, dateEnd, stanBands=False):
    """
    Download FISM2 daily data, either in the Standard Bands or not. Default is to download the data NOT in the
    Standard Bands. Note that this function downloads the data to a pre-determined location. The download is also
    'smart'. If the data already exists (over the exact same time period). It is NOT overwitten. If the requested data
    covers times OUTSIDE the time period of existing data, the NEW data is simply appended or prepended to a file
    containing the existing data.
    Args:
        dateStart: str
            The start date in 'YYYY-MM-DD' format.
        dateEnd: str
            The end date in 'YYYY-MM-DD' format.
        stanBands: bool
            If True, downloads FISM2 data in the STAN BANDS. Default is False.
    Returns:
        fismFile: str
            The location of the downloaded FISM data.
    """
    dateStartDatetime = datetime.strptime(dateStart, "%Y-%m-%d")
    dateEndDatetime = datetime.strptime(dateEnd, "%Y-%m-%d")

    # Helper function to manage the obtaining of data, given a URL:
    def urlObtain(url, fname):
        if os.path.isfile(fname) == True:
            print('File already exists (loading in data) '+str(fname))
        else:
            urlretrieve(url, fname)

    if stanBands:
        URL = 'https://lasp.colorado.edu/eve/data_access/eve_data/fism/daily_bands/daily_bands.nc'
        fname = 'FISM2_daily_stan_bands.nc'
        saveLoc = here.joinpath('../../tmp/')
        if os.path.exists(saveLoc) == False:
            os.makedirs(saveLoc)
        urlObtain(URL, saveLoc.joinpath(fname))
        datetimes, wavelengths, irradiance, uncertainties = readFism(saveLoc.joinpath(fname), stanBands=True)
    else:
        URL = 'https://lasp.colorado.edu/eve/data_access/eve_data/fism/daily_hr_data/daily_data.nc'
        fname = 'FISM2_daily_bands.nc'
        saveLoc = here.joinpath('../../tmp/')
        if os.path.exists(saveLoc) == False:
            os.makedirs(saveLoc)
        urlObtain(URL, saveLoc.joinpath(fname))
        datetimes, wavelengths, irradiance, uncertainties = readFism(saveLoc.joinpath(fname))

    # Subset the data in time, and save the subset data to a relative path:
    subset_inds = np.where((datetimes >= dateStartDatetime) & (datetimes <= dateEndDatetime))[0]
    subset_times = datetimes[subset_inds]
    subset_irradiance = irradiance[subset_inds, :]
    subset_uncertainty = uncertainties[subset_inds, :]

    return subset_times, wavelengths, subset_irradiance, subset_uncertainty

def rebin(fism_out, binning_scheme='EUVAC', zero=True):
    """
    Takes the output of getFism and rebins the data into whatever format the user desires.
    Args:
        fism_out: arraylike
            The output of getFism. Contains 4 elements: (1) datetime values for the FISM2 spectra, (2) the wavelengths
            of the spectrum, (3) the actual FISM2 irradiance spectra, and (4) uncertainties on the FISM2 irradiances (is
            non-NaN only for FISM2 data not in the Standard Bands).
        binning_scheme: str
            Determines the binning scheme to be used. Valid arguments include the following:
            'EUVAC': Uses the 37 wavelength band scheme described in Richards, et al. 1994; doi.org/10.1029/94JA00518
            'NEUVAC': Uses the 59 wavelength band scheme described in Brandt and Ridley, 2024; doi.org/10.1029/2024SW004043
            'HFG': Uses the 23 wavelength band scheme described in Solomon and Qian, 2005; https://doi.org/10.1029/2005JA011160
            'SOLOMON': Same situation as for argument 'HFG'.
            NOTE: If 'HFG' or 'SOLOMON' is chosen, the values of fism_out must correspond to getFism being run with the
            argument stanBands = True. If this IS NOT the case, an error will be thrown.
        zero: bool
            Controls whether singular (bright) wavelength lines are set to a value of zero after they are extracted.
            Default is True.
    Returns:
        fism2_file: str
            The location of the rebinned FISM data.
        fism2_data: arraylike
            Contains 3 elements: (a) a list of datetimes for the data and (b) the rebinned FISM2 data.
    """
    # Unpack the contents of fism_out:
    datetimes, wavelengths, irradiance, uncertainties = fism_out

    # Get the native wavelength resolution of the input data:
    # nativeResolution = np.concatenate((np.diff(wavelengths), np.array([np.diff(wavelengths)[-1]])), axis=0)
    nativeWavelengths = wavelengths.copy()

    if binning_scheme != 'HFG' and binning_scheme != 'SOLOMON':
        if binning_scheme == 'EUVAC':
            # Grab the euv_37.csv file:
            fileStr = str(here.joinpath('euv_files/euv_37.csv'))
            bin_bounds = read_euv_csv_file(fileStr)
            tag = '_37'
        elif binning_scheme == 'NEUVAC':
            # Grab the euv_59.csv file:
            fileStr = str(here.joinpath('euv_files/euv_59.csv'))
            bin_bounds = read_euv_csv_file(fileStr)
            tag = '_59'
        else:
            raise FileNotFoundError('The .csv files for specifying bin boundaries cannot be found!!')

        # Perform the rebinning!
        shorts = bin_bounds['short'] / 10.
        longs = bin_bounds['long'] / 10.
        newWaves = 0.5 * (shorts + longs)

        # Instantiate the new data array:
        if len(irradiance.shape) < 2:
            fism2_data = np.zeros((1, newWaves.shape[0]))
        else:
            fism2_data = np.zeros((irradiance.shape[0], newWaves.shape[0]))

        # First go through all the wavelengths that are singular
        myData = irradiance
        for iWave, short in enumerate(shorts):
            long = longs[iWave]
            if (long == short):
                i = np.argmin(np.abs(wavelengths - short))
                i2 = np.argmin(np.abs(nativeWavelengths - short))
                try:
                    fism2_data[:, iWave] = myData[:, i] * (nativeWavelengths[i2 + 1] - nativeWavelengths[i2])
                except:
                    fism2_data[:, iWave] = myData[i] * (nativeWavelengths[i2 + 1] - nativeWavelengths[i2])
                if zero == True:
                    # Zero out bin so we don't double count it.
                    try:
                        myData[:, i] = np.zeros_like(myData[:, i])
                    except:
                        myData[i] = 0.0

        # Then go through the ranges
        for iWave, short in enumerate(shorts):
            long = longs[iWave]
            if (long != short):
                d1 = np.abs(wavelengths - short)
                iStart = np.argmin(d1)
                d2 = np.abs(wavelengths - long)
                iEnd = np.argmin(d2)
                wave_int = 0.0
                # For wavelengths at or below 0.2 nm, just compute the sum:
                if long <= 0.2:
                    for i in range(iStart + 1, iEnd + 1):
                        fism2_data[:, iWave] += myData[:, i] * \
                                             (wavelengths[i + 1] - wavelengths[i])
                        wave_int += (wavelengths[i + 1] - wavelengths[i])
                else:
                    # For issues computing the sum, integrate instead:
                    try:
                        fism2_data[:, iWave] = integ.trapezoid(myData[:, iStart:iEnd], wavelengths[iStart:iEnd], axis=1)
                    except:
                        fism2_data[:, iWave] = integ.trapezoid(myData[iStart:iEnd], wavelengths[iStart:iEnd])

                    # # Plotting for a sanity check:
                    # plt.figure()
                    # plt.plot(wavelengths[iStart:iEnd], myData[0, iStart:iEnd], marker='o')
                    # plt.scatter(newWaves, fism2_data[0, :])
                    # plt.show()

    elif binning_scheme == 'HFG' or binning_scheme == 'SOLOMON':
        # Determine whether the supplied data already conforms to the Solomon and Qian binning scheme.
        tag = '_solomon'
        if fism_out[2].shape[1] != 23:
            raise ValueError("Incorrect dimensions for element 3 of argument 'fism_out'. Dimensions must be (n,23), "
                             "resulting from running function 'getFism' with argument 'stanBands'=True. ")
        # Should the data confirm to the proper dimensions, there is no rebinning step that needs to be done. Simply
        # continue.
        fism2_data = fism_out[2]
    else:
        # If the input irradiance data DOES NOT conform to the Solomon and Qian binning scheme, throw an error.
        raise ValueError("Invalid value for argument 'binning_scheme'. Must be 'EUVAC', 'NEUVAC', 'HFG', or 'SOLOMON'.")

    # Save the rebinned data to a relative path (outside the package directory) in the form of a .txt file:
    fism2_file = here.joinpath('../../tmp/fism2_file'+tag+'.txt')
    saveFism(fism2_data, datetimes, fism2_file)

    return fism2_file, fism2_data

def saveFism(data, times, filename):
    """
    Takes (rebinned) FISM2 data and saves it a .txt file at a user-defined location.
    Args:
        data: numpy.ndarray
            Irradiance data in a nxm array where the first dimension corresponds to observations (the spectrum number)
            and the second dimension corresponds to wavelengths.
        times: numpy.ndarray
            The time values at which each spectrum is recorded.
        filename: str
            The desired location where the data will be saved.
    Returns:
        Nothing. Simply saves a file.
    """
    # A helper function for working with integers:
    def numStr(num):
        if int(num) < 10:
            return ' ' + str(int(num))
        else:
            return str(int(num))

    # Define a helper function for opening a file to write the data, in such a way as to include parent directories if
    # needed:
    def safe_open_w(path):
        ''' Open "path" for writing, creating any parent directories as needed.
        (https://stackoverflow.com/questions/23793987/write-a-file-to-a-directory-that-doesnt-exist)
        '''
        os.makedirs(os.path.dirname(path), exist_ok=True)
        return open(path, 'w')

    # Open the new file and begin writing, line by line:
    with safe_open_w(str(filename)) as output:
        # Write the header information:
        output.write("#START\n")
        # Write the irradiances themselves:
        firstLine = ['%.6g' % (element) for element in data[0, :]]
        firstLine_joined = ' '.join(firstLine)
        # The first line should always be a duplicate of the first line of data, but starting at UTC=00:00 of the first date:
        output.write(' ' + str(times[0].year) + ' ' + numStr(times[0].month) + ' ' + numStr(
            times[0].day) + '  0  0  0 ' + firstLine_joined + '\n')
        # The rest of the lines can be straight from the data:
        for i in range(data.shape[0]):
            currentLine_joined = ' '.join(['%.6g' % (element) for element in data[i, :]])
            output.writelines(' ' + str(times[i].year) + ' ' + numStr(times[i].month) + ' ' + numStr(
                times[i].day) + ' ' + numStr(times[i].hour) + '  0  0 ' + currentLine_joined + '\n')
        # The last line should occur 12 hours from the last datapoint, but have duplicate values there:
        lastLine_joined = ' '.join(['%.6g' % (element) for element in data[-1, :]])
        lastTime = times[-1] + timedelta(hours=12)
        output.write(' ' + str(lastTime.year) + ' ' + numStr(lastTime.month) + ' ' + numStr(
            lastTime.day) + '  0  0  0 ' + lastLine_joined + '\n')

    print('FISM2 data saved to: ')
    os.system('readlink -f '+str(filename))
    return

def readFism(fism_file, stanBands=False):
    """
    Load in spectrum data from a FISM2 file.
    Args:
        fism_file: str
            The location of a FISM2 NETCDF4 file.
        stanBands: bool
            If True, expects the data to be read in to be in the STAN BANDS. Default is False.
    Returns:
        datetimes: numpy.ndarray
            An array of datetimes for reach FISM2 spectrum.
        wavelengths: numpy.ndarray
            A one-dimensional array of wavelengths at which there are irradiance values.
        irradiances: numpy.ndarray
            A two-dimensional array of irradiance values in each wavelength bin.
        uncertainties: numpy.ndarray
            A two-dimensional array of irradiance uncertainty values at each bin. Only returned if stanBands is False.
    """
    fism2Data = Dataset(fism_file)
    wavelengths = np.asarray(fism2Data.variables['wavelength'])
    if stanBands:
        flux = np.asarray(fism2Data.variables['ssi']) # photons/cm^2/second
        pFlux = flux * 1e4 # photons/m^2/second
        # Convert fluxes to irradiances:
        irradiance = np.zeros_like(flux)
        for i in range(flux.shape[1]):
            irradiance[:, i] = spectralIrradiance(pFlux[:, i], wavelengths[i] * 10.)# W/m^2
        uncertainties = np.full_like(irradiance, np.nan)
    else:
        irradiance = np.asarray(fism2Data.variables['irradiance'])
        uncertainties = np.asarray(fism2Data.variables['uncertainty'])
    dates = fism2Data.variables['date']
    datetimes = []
    for j in range(len(dates)):
        year = dates[j][:4]
        day = dates[j][4:]
        currentDatetime = datetime(int(year), 1, 1) + timedelta(int(day) -1 ) + timedelta(hours=12)
        datetimes.append(currentDatetime)
    datetimes = np.asarray(datetimes)

    return datetimes, wavelengths, irradiance, uncertainties

def spectralIrradiance(photonFlux, wavelength):
    """
    Convert the photon flux to the corresponding spectral irradiance, given a specific wavelength.
    Args:
        photonFlux: numpy.ndarray, float, or int
            Photon flux in units of photons s^-1 m^-2. For a singular wavelength, units are in photons m^-2.
        wavelength: wavelength: float
            A specific wavelength in Angstroms.
    Returns:
        irradiance: numpy.ndarray or float
            The corresponding spectral irradiance in units of W/m^2/nm.
    """
    photonEnergy = (h*c) / (wavelength*1e-10) # Convert the wavelength in the denominator to meters.
    irradiance= photonFlux * photonEnergy
    return irradiance

def read_euv_csv_file(file):
    """
    Originally written by Aaron J. Ridley, within the file 'fism2_process.py':
    https://github.com/aaronjridley/EUV/blob/main/fism2_process.py

    This file reads in binning data from a CSV file that specifies bin boundaries and cross sections for either the
    EUVAC model or the NEUVAC model.
    Args:
        file: str
            The location of the .csv file to be read.
    Returns:
        wavelengths: numpy.ndarray
            The wavelength bin boundaries for either the EUVAC model or the NEUVAC model.
    """
    fpin = open(file, 'r')

    iFound = 0
    afac = []
    f74113 = []
    for line in fpin:
        aline = line.split(',')
        s = aline[-1].strip().split('.')[0]
        if (aline[0].strip() == "Short"):
            if (s.isnumeric()):
                short = np.asarray(aline[5:], dtype=float)
            else:
                short = np.asarray(aline[5:-1], dtype=float)
            iFound += 1
        if (aline[0].strip() == "Long"):
            if (s.isnumeric()):
                long = np.asarray(aline[5:], dtype=float)
            else:
                long = np.asarray(aline[5:-1], dtype=float)
        if (aline[0].strip() == "F74113"):
            if (s.isnumeric()):
                f74113 = np.asarray(aline[5:], dtype=float)
            else:
                f74113 = np.asarray(aline[5:-1], dtype=float)
            iFound += 1
        if (aline[0].strip() == "AFAC"):
            if (s.isnumeric()):
                afac = np.asarray(aline[5:], dtype=float)
            else:
                afac = np.asarray(aline[5:-1], dtype=float)
            iFound += 1
    # Save and convert from Angstroms to nm (FISM is in nm)
    wavelengths = {'short': short / 10.0,
                   'long': long / 10.0,
                   'afac': afac,
                   'f74113': f74113}
    return wavelengths

# Execution (testing):
if __name__ == '__main__':
    # Download some FISM2 data for the time period stated by the user.
    dateStart = '2015-08-13'
    dateEnd = '2015-08-19'
    # fismOut = getFism(dateStart, dateEnd)
    # rebinnedFismFile, rebinnedFismData = rebin(fismOut, binning_scheme='EUVAC')
    # rebinnedFismFile_N, rebinnedFismData_N = rebin(fismOut, binning_scheme='NEUVAC')

    fismOut_S = getFism(dateStart, dateEnd, stanBands=True)
    rebinnedFismFile_S, rebinnedFismData_S = rebin(fismOut_S, binning_scheme='SOLOMON')

