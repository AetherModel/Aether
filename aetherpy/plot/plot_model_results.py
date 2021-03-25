#!/usr/bin/env python
""" Standard model visualization routines
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os

from aetherpy.io import read_routines
from aetherpy.utils import inputs
from aetherpy.plot import movie_routines 


#-----------------------------------------------------------------------------
# Define the support routines

def get_help(file_vars=None):
    """ Provide string explaining how to run the command line interface

    Parameters
    ----------
    file_vars : list or NoneType
        List of file variables or None to exclude this output (default=None)

    Returns
    -------
    help_str : str
        String with formatted help statement

    """

    help_str = 'Usage:\n{:s} -[flags] [filenames]\n'.format(__name__)
    help_str += 'Flags:\n'
    help_str += '       -help : print this message, include filename for '
    help_str += 'variable names and indices\n'
    help_str += '       -var=number : index of variable to plot\n'
    help_str += '       -cut=alt, lat, or lon : which cut you would like\n'
    help_str += '       -alt=number : alt in km or grid number (closest)\n'
    help_str += '       -lat=number : latitude in degrees (closest)\n'
    help_str += '       -lon=number: longitude in degrees (closest)\n'
    help_str += '       -log : plot the log of the variable\n'
    help_str += '       -winds : overplot winds\n'
    help_str += '       -tec : plot the TEC variable\n'
    help_str += '       -movie=number : provide a positive frame rate to '
    help_str += 'create a movie\n'
    help_str += '       -ext=str : figure or movie extension\n'
    help_str += 'At end, list the files you want to plot. This code should '
    help_str += 'work with either GITM files (*.bin) or Aether netCDF files '
    help_str += '(*.nc)'

    if file_vars is not None:
        help_str += "File Variables (index, name):\n"
        for ivar, var in enumerate(file_vars):
            help_str += "               ({:d}, {:s})\n".format(ivar, var)

    return help_str


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
        filelist (list of filenames), gitm (flag that is true for GITM input,
        determined by examining filelist naming convention),
        var (variable index to plot), cut (coordinate to hold constant),
        diff (difference with other plots),
        movie (framerate for movie, which is > 0 if a movie is desired),
        ext (output extension), winds (flag to plot with winds),
        alt (to plot), lat (to plot), lon (to plot),
        log (flag to use log scale), and help (flag to display help)

    """
    # Initialize the arguments to their default values
    args = {'filelist': [], 'log': False, 'var': 15, 'alt': 400,
            'lon': np.nan, 'lat': np.nan, 'cut': 'alt', 'help': False,
            'winds': False, 'diff': False, 'gitm': False, 'movie': 0,
            'ext': 'png'}

    arg_type = {'filelist': list, 'log': bool, 'var': int, 'alt': int,
                'lon': float, 'lat': float, 'cut': str, 'help': bool,
                'winds': bool, 'diff': bool, 'gitm': bool,
                'movie': int, 'ext': str}

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

            if len(split_arg) == 0:
                if arg_type[akey] == bool:
                    arg_val = True
                else:
                    raise ValueError('expected equality after flag {:}'.format(
                        akey))
            else:
                if arg_type[akey] == int:
                    arg_val = int(split_arg[1])
                elif arg_type[akey] == float:
                    arg_val] = float(split_arg[1])
                elif arg_type[akey] == str:
                    arg_val = split_arg[1]
                else:
                    # This is boolean input
                    arg_val = inputs.bool_string(split_arg[1])

            # Assign the output
            if akey.find('tec') == 0:
                args['var'] = 34
            else:
                args[akey] = arg_val
        else:
            # Save the filenames
            args['filelist'].append(arg)

            gitm = arg.find('bin') == len(arg) - 3

            if len(args['filelist']) == 1:
                args['gitm'] = gitm
            elif gitm != args['gitm']:
                raise ValueError('input files are of a mixed type')

    # Update default movie extention for POSIX systems
    if args['movie'] > 0 and args['ext'] == 'png':
        if (os.name == "posix"):
            args['ext'] = "mkv"
        else:
            args['ext'] = "mp4"

    return args


#-----------------------------------------------------------------------------
# Define the main routine

def main():
    # Get the input arguments
    args = get_command_line_args(inputs.process_command_line_input())

    if args['help'] and len(args['filelist']) == 0:
        help_str = get_help()
        print(help_str)
        return

    # Read the file header
    if args['gitm']:
        header = read_routines.read_gitm_header(args["filelist"])
    else:
        header = aether_routines.read_aether_headers(args["filelist"])

    # If help is requested for a specific file, return it here
    if args['help']:
        help_str = get_help(header['vars'])
        print(help_str)
        return

    if args["var"] >= len(header["vars"]):
        raise ValueError("requested variable doesn't exist: {:d}>{:d}".format(
            args["var"], len(header["vars"])))

    # Define the plotting inputs
    plot_vars = [0, 1, 2, args["var"]]

if (args["winds"]):
    if (args['cut']=='alt'):
        iUx_ = 16
        iUy_ = 17
    if (args['cut']=='lat'):
        iUx_ = 16
        iUy_ = 18
    if (args['cut']=='lon'):
        iUx_ = 17
        iUy_ = 18
    plot_vars.append(iUx_)
    plot_vars.append(iUy_)
    AllWindsX = []
    AllWindsY = []

Var = header["vars"][args["var"]]

iVar_ = args["var"]

AllData2D = []
AllAlts = []
AllTimes = []

j = 0
iCut = -1
for file in args['filelist']:

    if args['gitm']:
        data = read_gitm_one_file(file, plot_vars)
        iVar_ = args["var"]
    else:
        if (j == 0):
            VarList = []
            for v in plot_vars:
                VarList.append(header["vars"][v])
        data = read_aether_one_file(file, VarList)
        iVar_ = 3
    if (j == 0):
        [nLons, nLats, nAlts] = data[0].shape
        Alts = data[2][0][0]/1000.0;
        Lons = np.degreed(data[0][:,0,0])
        Lats = np.degreed(data[1][0,:,0])
        if (args['cut'] == 'alt'):
            xPos = Lons
            yPos = Lats
            if (len(Alts) > 1):
                if (args["alt"] < 50):
                    iAlt = args["alt"]
                else:
                    if (args["alt"] > Alts[nAlts-3]):
                        iAlt = nAlts-3
                    else:
                        iAlt = 2
                        while (Alts[iAlt] < args["alt"]):
                            iAlt=iAlt+1
            else:
                iAlt = 0
            Alt = Alts[iAlt]
            iCut = iAlt
            
        if (args['cut'] == 'lat'):
            xPos = Lons
            yPos = Alts
            if (args["lat"] < Lats[1]):
                iLat = int(nLats/2)
            else:
                if (args["lat"] > Lats[nLats-2]):
                    iLat = int(nLats/2)
                else:
                    iLat = 2
                    while (Lats[iLat] < args["lat"]):
                        iLat=iLat+1
            Lat = Lats[iLat]
            iCut = iLat
            
        if (args['cut'] == 'lon'):
            xPos = Lats
            yPos = Alts
            if (args["lon"] < Lons[1]):
                iLon = int(nLons/2)
            else:
                if (args["lon"] > Lons[nLons-2]):
                    iLon = int(nLons/2)
                else:
                    iLon = 2
                    while (Lons[iLon] < args["lon"]):
                        iLon=iLon+1
            Lon = Lons[iLon]
            iCut = iLon
                        
    AllTimes.append(data["time"])
    
    if (args["tec"]):
        iAlt = 2
        tec = np.zeros((nLons, nLats))
        for Alt in Alts:
            if (iAlt > 0 and iAlt < nAlts-3):
                tec = tec + data[iVar_][:,:,iAlt] * (Alts[iAlt+1]-Alts[iAlt-1])/2 * 1000.0
            iAlt=iAlt+1
        AllData2D.append(tec/1e16)
    else:
        if (args['cut'] == 'alt'):
            AllData2D.append(data[iVar_][:,:,iAlt])
        if (args['cut'] == 'lat'):
            AllData2D.append(data[iVar_][:,iLat,:])
        if (args['cut'] == 'lon'):
            AllData2D.append(data[iVar_][iLon,:,:])
        if (args["winds"]):
            if (args['cut'] == 'alt'):
                AllWindsX.append(data[iUx_][:,:,iAlt])
                AllWindsY.append(data[iUy_][:,:,iAlt])
            if (args['cut'] == 'lat'):
                AllWindsX.append(data[iUx_][:,iLat,:])
                AllWindsY.append(data[iUy_][:,iLat,:])
            if (args['cut'] == 'lon'):
                AllWindsX.append(data[iUx_][iLon,:,:])
                AllWindsY.append(data[iUy_][iLon,:,:])
    j=j+1
    

AllData2D = np.array(AllData2D)
if (args['IsLog']):
    AllData2D = np.log10(AllData2D)
if (args["winds"]):
    AllWindsX = np.array(AllWindsX)
    AllWindsY = np.array(AllWindsY)

Negative = 0

maxi  = np.max(AllData2D)*1.01
mini  = np.min(AllData2D)*0.99

if (mini < 0):
    Negative = 1

if (Negative):
    maxi = np.max(abs(AllData2D))*1.05
    mini = -maxi

if (args['cut'] == 'alt'):
    maskNorth = ((yPos>45) & (yPos<90.0))
    maskSouth = ((yPos<-45) & (yPos>-90.0))
    DoPlotNorth = np.max(maskNorth)
    DoPlotSouth = np.max(maskSouth)
    if (DoPlotNorth):
        maxiN = np.max(abs(AllData2D[:,:,maskNorth]))*1.05
        if (Negative):
            miniN = -maxiN
        else:
            miniN = np.min(AllData2D[:,:,maskNorth])*0.95
    if (DoPlotSouth):
        maxiS = np.max(abs(AllData2D[:,:,maskSouth]))*1.05
        if (Negative):
            miniS = -maxiS
        else:
            miniS = np.min(AllData2D[:,:,maskSouth])*0.95
    
dr = (maxi-mini)/31
levels = np.arange(mini, maxi, dr)

i = 0

# Define plot range:
minX = (xPos[ 1] + xPos[ 2])/2
maxX = (xPos[-2] + xPos[-3])/2
minY = (yPos[ 1] + yPos[ 2])/2
maxY = (yPos[-2] + yPos[-3])/2

file = "var%2.2d_" % args["var"]
file = file+args['cut']+"%3.3d" % iCut

if args['movies'] > 0:
    img_file_fmt = movie_routines.setup_movie_dir(file)
else:
    img_file_fmt = "".join(file, '_', "{:}.", args['ext'])

iter = 0
for time in AllTimes:

    ut = time.hour + time.minute/60.0 + time.second/3600.0
    shift = ut * 15.0
    
    fig = plt.figure(constrained_layout=False,
                     tight_layout=True, figsize=(10, 8.5))

    gs1 = mpl.gridspec.GridSpec(nrows=2, ncols=2, wspace=0.0, hspace=0)
    gs = mpl.gridspec.GridSpec(nrows=2, ncols=2, wspace=0.0, left=0.0,
                               right=0.9)

    norm = mpl.cm.colors.Normalize(vmax=mini, vmin=maxi)
    if (mini >= 0):
        cmap = mpl.cm.plasma
    else:
        cmap = mpl.cm.bwr

    d2d = np.transpose(AllData2D[iter])
    if (args["winds"]):
        Ux2d = np.transpose(AllWindsX[iter])
        Uy2d = np.transpose(AllWindsY[iter])

    fmt_input = iter if args['movie'] > 0 else time.strftime('%y%m%d_%H%M%S')
    outfile = img_file_fmt.format(fmt_input)

    ax = fig.add_subplot(gs1[1, :2])

    cax = ax.pcolor(xPos, yPos, d2d, vmin=mini, vmax=maxi, cmap=cmap)

    if (args["winds"]):
        ax.quiver(xPos,yPos,Ux2d,Uy2d)
    ax.set_ylim([minY,maxY])
    ax.set_xlim([minX,maxX])

    if (args['cut'] == 'alt'):
        ax.set_ylabel('Latitude (deg)')
        ax.set_xlabel('Longitude (deg)')
        title = time.strftime('%b %d, %Y %H:%M:%S')+'; Alt : '+"%.2f" % Alt + ' km'
        ax.set_aspect(1.0)

    if (args['cut'] == 'lat'):
        ax.set_xlabel('Longitude (deg)')
        ax.set_ylabel('Altitude (km)')
        title = time.strftime('%b %d, %Y %H:%M:%S')+'; Lat : '+"%.2f" % Lat + ' km'

    if (args['cut'] == 'lon'):
        ax.set_xlabel('Latitude (deg)')
        ax.set_ylabel('Altitude (km)')
        title = time.strftime('%b %d, %Y %H:%M:%S')+'; Lon : '+"%.2f" % Lon + ' km'

    ax.set_title(title)
    cbar = fig.colorbar(cax, ax=ax, shrink = 0.75, pad=0.02)
    cbar.set_label(Var,rotation=90)

    if (args['cut'] == 'alt'):
        
        if (DoPlotNorth):
            # Top Left Graph Northern Hemisphere
            ax2 = fig.add_subplot(gs[0, 0],projection='polar')
            r, theta = np.meshgrid(90.0-yPos[maskNorth], \
                                   (xPos+shift-90.0)*3.14159/180.0)
            cax2 = ax2.pcolor(theta, r, AllData2D[iter][:,maskNorth], \
                              vmin=miniN, vmax=maxiN, cmap=cmap)
            xlabels = ['', '12', '18', '00']
            ylabels = ['80', '70', '60', '50']
            ax2.set_xticklabels(xlabels)
            ax2.set_yticklabels(ylabels)
            cbar2 = fig.colorbar(cax2, ax=ax2, shrink = 0.5, pad=0.01)
            ax2.grid(linestyle=':', color='black')
            pi = 3.14159
            ax2.set_xticks(np.arange(0,2*pi,pi/2))
            ax2.set_yticks(np.arange(10,50,10))
            
        if (DoPlotSouth):
            # Top Right Graph Southern Hemisphere
            r, theta = np.meshgrid(90.0+yPos[maskSouth], \
                                   (xPos+shift-90.0)*3.14159/180.0)
            ax3 = fig.add_subplot(gs[0, 1],projection='polar')
            cax3 = ax3.pcolor(theta, r, AllData2D[iter][:,maskSouth], \
                              vmin=miniS, vmax=maxiS, cmap=cmap)
            xlabels = ['', '12', '18', '00']
            ylabels = ['80', '70', '60', '50']
            ax3.set_xticklabels(xlabels)
            ax3.set_yticklabels(ylabels)
            cbar3 = fig.colorbar(cax3, ax=ax3, shrink = 0.5, pad=0.01)
            ax3.grid(linestyle=':', color='black')
            pi = 3.14159
            ax3.set_xticks(np.arange(0,2*pi,pi/2))
            ax3.set_yticks(np.arange(10,50,10))

    print("Writing file : "+outfile)
    fig.savefig(outfile)
    plt.close()

    iter = iter + 1

if args['movie'] > 0:
    movie_routines.save_movie(file, ext=args['ext'], rate=args['movie'])
