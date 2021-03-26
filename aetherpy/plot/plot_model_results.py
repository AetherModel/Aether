#!/usr/bin/env python
""" Standard model visualization routines
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os

from aetherpy.io import read_routines
from aetherpy.utils import inputs
from aetherpy.plot import data_prep, movie_routines


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

    # Update plotting variables to include the wind, if desired
    if args["winds"]:
        plot_vars.append(16 if args['cut'] in ['alt', 'lat'] else 17)
        plot_vars.append(18 if args['cut'] in ['lat', 'lon'] else 17)
        all_winds_x = []
        all_winds_y = []

    # Prepare to load the desired file data
    var = header["vars"][args["var"]]
    ivar = args["var"]

    all_2dim_data = []
    all_times = []

    for j, filename in enumerate(args['filelist']):
        # Read in the data file
        if args['gitm']:
            data = read_routines.read_gitm_file(filename, plot_vars)
            ivar = args["var"]
        else:
            if j == 0:
                var_list = []
                for pvar in plot_vars:
                    var_list.append(header["vars"][pvar])
            data = read_routines.read_aether_file(filename, var_list)
            ivar = 3

        # For the first file, initialize the necessary plotting data
        if j == 0:
            # Get 1D arrays for the coordinates
            alts = data[2][0][0] / 1000.0  # Convert from m to km
            lons = np.degrees(data[0][:, 0, 0])  # Convert from rad to deg
            lats = np.degrees(data[1][0, :, 0])  # Convert from rad to deg

            # Find the desired index to cut along to get a 2D slice
            icut, cut_data, x_pos, y_pos, z_val = data_prep.get_cut_index(
                lons, lats, alts, args[args['cut']], args['cut'])


        # Save the time data
        all_times.append(data["time"])

        # Save the z-axis data
        if args["tec"]:
            all_2dim_data.append(data_prep.calc_tec(alts, data[ivar], 2, -4))
        else:
            all_2dim_data.append(data[ivar][cut_data])

            if (args["winds"]):
                all_winds_x.append(data[plot_vars[-1]][cut_data])
                all_winds_y.append(data[plot_vars[-1]][cut_data])

    # Convert data list to a numpy array
    all_2dim_data = np.array(all_2dim_data)

    if args["winds"]:
        all_winds_x = np.array(all_winds_x)
        all_winds_y = np.array(all_winds_y)

    # If desired, take the log of the data
    if args['log']:
        all_2dim_data = np.log10(all_2dim_data)

# HERE
Negative = 0

maxi  = np.max(all_2dim_data)*1.01
mini  = np.min(all_2dim_data)*0.99

if (mini < 0):
    Negative = 1

if (Negative):
    maxi = np.max(abs(all_2dim_data))*1.05
    mini = -maxi

if (args['cut'] == 'alt'):
    maskNorth = ((y_pos>45) & (y_pos<90.0))
    maskSouth = ((y_pos<-45) & (y_pos>-90.0))
    DoPlotNorth = np.max(maskNorth)
    DoPlotSouth = np.max(maskSouth)
    if (DoPlotNorth):
        maxiN = np.max(abs(all_2dim_data[:,:,maskNorth]))*1.05
        if (Negative):
            miniN = -maxiN
        else:
            miniN = np.min(all_2dim_data[:,:,maskNorth])*0.95
    if (DoPlotSouth):
        maxiS = np.max(abs(all_2dim_data[:,:,maskSouth]))*1.05
        if (Negative):
            miniS = -maxiS
        else:
            miniS = np.min(all_2dim_data[:,:,maskSouth])*0.95
    
dr = (maxi-mini)/31
levels = np.arange(mini, maxi, dr)

i = 0

# Define plot range:
minX = (x_pos[ 1] + x_pos[ 2])/2
maxX = (x_pos[-2] + x_pos[-3])/2
minY = (y_pos[ 1] + y_pos[ 2])/2
maxY = (y_pos[-2] + y_pos[-3])/2

filename = "var%2.2d_" % args["var"]
filename = filename+args['cut']+"%3.3d" % icut

if args['movies'] > 0:
    img_file_fmt = movie_routines.setup_movie_dir(filename)
else:
    img_file_fmt = "".join(filename, '_', "{:}.", args['ext'])

iter = 0
for time in all_times:

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

    d2d = np.transpose(all_2dim_data[iter])
    if (args["winds"]):
        Ux2d = np.transpose(all_winds_x[iter])
        Uy2d = np.transpose(all_winds_y[iter])

    fmt_input = iter if args['movie'] > 0 else time.strftime('%y%m%d_%H%M%S')
    outfile = img_file_fmt.format(fmt_input)

    ax = fig.add_subplot(gs1[1, :2])

    cax = ax.pcolor(x_pos, y_pos, d2d, vmin=mini, vmax=maxi, cmap=cmap)

    if (args["winds"]):
        ax.quiver(x_pos,y_pos,Ux2d,Uy2d)
    ax.set_ylim([minY,maxY])
    ax.set_xlim([minX,maxX])

    if (args['cut'] == 'alt'):
        ax.set_ylabel('Latitude (deg)')
        ax.set_xlabel('Longitude (deg)')
        title = time.strftime('%b %d, %Y %H:%M:%S')+'; alt : '+"%.2f" % z_val + ' km'
        ax.set_aspect(1.0)

    if (args['cut'] == 'lat'):
        ax.set_xlabel('Longitude (deg)')
        ax.set_ylabel('Altitude (km)')
        title = time.strftime('%b %d, %Y %H:%M:%S')+'; Lat : '+"%.2f" % z_val + ' km'

    if (args['cut'] == 'lon'):
        ax.set_xlabel('Latitude (deg)')
        ax.set_ylabel('Altitude (km)')
        title = time.strftime('%b %d, %Y %H:%M:%S')+'; Lon : '+"%.2f" % z_val + ' km'

    ax.set_title(title)
    cbar = fig.colorbar(cax, ax=ax, shrink = 0.75, pad=0.02)
    cbar.set_label(var,rotation=90)

    if (args['cut'] == 'alt'):
        
        if (DoPlotNorth):
            # Top Left Graph Northern Hemisphere
            ax2 = fig.add_subplot(gs[0, 0],projection='polar')
            r, theta = np.meshgrid(90.0-y_pos[maskNorth], \
                                   (x_pos+shift-90.0)*3.14159/180.0)
            cax2 = ax2.pcolor(theta, r, all_2dim_data[iter][:,maskNorth], \
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
            r, theta = np.meshgrid(90.0+y_pos[maskSouth], \
                                   (x_pos+shift-90.0)*3.14159/180.0)
            ax3 = fig.add_subplot(gs[0, 1],projection='polar')
            cax3 = ax3.pcolor(theta, r, all_2dim_data[iter][:,maskSouth], \
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
    movie_routines.save_movie(filename, ext=args['ext'], rate=args['movie'])
