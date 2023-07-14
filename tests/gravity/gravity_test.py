from glob import glob
import os
import re
import sys
from xml.dom import minicompat
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.colors as colors
import numpy as np
import cartopy.crs as ccrs

from aetherpy import logger
import read_routines
from aetherpy.plot import data_prep
from aetherpy.plot import movie_routines
from aetherpy.utils import inputs
from aetherpy.utils import time_conversion


np.set_printoptions(threshold=sys.maxsize)

def read_nc_files(filename):
    header = read_routines.read_blocked_netcdf_header(filename)
    data = read_routines.read_blocked_netcdf_file(filename)
    return header,data




def compare_gravity(argv): 
    pi = 3.141592653589793
    dtor =  3.141592653589793/180
    header,data = read_nc_files(argv[1])

    galtitude = data['Gvertical'][0]
    glongitude = data['Geast'][0]
    glatitude = data['Gnorth'][0]

    gpotential = data['Gpotential'][0]

    longitude = data['lon'][0] * dtor
    latitude = data['lat'][0] * dtor


    radius = data['radius'][0]
    z = data['z'][0]
    

    polar_radius = 6356.8 * 1000 
    equator_radius = 6378.1 * 1000
 

    delta_radius = equator_radius - polar_radius

    radius_at_45 = polar_radius + (delta_radius * np.cos( latitude[0][5][0]))

    equator_radius = polar_radius + (delta_radius * np.cos(latitude[0][11][0]))

    polar_radius = polar_radius + (delta_radius * np.cos(latitude[0][0][0]))


    polar_radius += 90000
    equator_radius += 90000
    radius_at_45 += 90000 

    cg = 6.674080039e-11 
    mass = 5.97e+24
    mu = mass * cg


    gravity_potential_polar = ( -mu/ polar_radius)  + (((3*(0.00108262668 * mu) ) / (( 2 * (polar_radius**3)))) * ((np.sin(latitude[0][0])**2 ) - 1.0))

    gravity_potential_equator = ( -mu/ equator_radius)  + (((3*(0.00108262668 * mu) ) / (( 2 * (equator_radius**3))))     * ((np.sin(latitude[0][11])**2 ) - 1.0))

    gravity_potential_45 = ( -mu/ radius_at_45)  + (((3*(0.00108262668 * mu) ) / (( 2 * (radius_at_45**3)))) * ((np.sin(latitude[0][5])**2 ) - 1.0))
    
    print("Aether potential polar:  ", gpotential[0][0][0])
    print("\n")
    print("Calculated polar potential",gravity_potential_polar[0])
    print("\n")
    print("Aether potential equator:  ", gpotential[0][11][0])
    print("\n")
    print("Calculated equator potential",gravity_potential_equator[0])
    print("\n")
    print("Aether potential 45:  ", gpotential[0][5][0])
    print("\n")
    print("Calculated 45 potential",gravity_potential_45[0])


if __name__ == '__main__':
    compare_gravity(sys.argv)

