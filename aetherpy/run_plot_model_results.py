#!/usr/bin/env python
""" Standard model visualization routines
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os

from aetherpy.io import read_routines
from aetherpy.utils import inputs, time_conversion
from aetherpy.plot import data_prep, movie_routines


# ----------------------------------------------------------------------------
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

    mname = os.path.join(
        os.path.commonpath([inputs.__file__, data_prep.__file__]),
        'run_plot_model_results.py') if __name__ == '__main__' else __name__

    help_str = 'Usage:\n{:s} -[flags] [filenames]\n'.format(mname)
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
    args = {'filelist': [], 'log': False, 'var': 15, 'alt': 400, 'tec': False,
            'lon': np.nan, 'lat': np.nan, 'cut': 'alt', 'winds': False,
            'diff': False, 'gitm': False, 'movie': 0, 'ext': 'png'}

    arg_type = {'filelist': list, 'log': bool, 'var': int, 'alt': int,
                'tec': bool,
                'lon': float, 'lat': float, 'cut': str, 'help': bool,
                'winds': bool, 'diff': bool, 'gitm': bool, 'tec': bool,
                'movie': int, 'ext': str}

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


# ----------------------------------------------------------------------------
# Define the main plotting routine

def plot_model_results():
    # Get the input arguments
    args = get_command_line_args(inputs.process_command_line_input())

    if len(args['filelist']) == 0:
        help_str = get_help()
        print(help_str)
        return

    # Read the file header
    if args['gitm']:
        header = read_routines.read_gitm_headers(args["filelist"])
    else:
        files = args["filelist"]
        header = read_routines.read_aether_headers(files)

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

    # Define plotting limits
    symmetric = False
    cmap = mpl.cm.plasma

    maxi = all_2dim_data.max() * 1.01
    mini = all_2dim_data.min() * 0.99

    if mini < 0.0:
        symmetric = True
        cmap = mpl.cm.bwr
        maxi = abs(all_2dim_data).max() * 1.05
        mini = -maxi

    if args['cut'] == 'alt':
        mask_north = ((y_pos > 45) & (y_pos < 90.0))
        mask_south = ((y_pos < -45) & (y_pos > -90.0))
        plot_north = mask_north.max()
        plot_south = mask_south.max()

        if plot_north:
            maxi_north = abs(all_2dim_data[:, :, mask_north]).max() * 1.05

            if symmetric:
                mini_north = -maxi_north
            else:
                mini_north = all_2dim_data[:, :, mask_north].min() * 0.95

        if plot_south:
            maxi_south = abs(all_2dim_data[:, :, mask_south]).max() * 1.05

            if symmetric:
                mini_south = -maxi_south
            else:
                mini_south = all_2dim_data[:, :, mask_south].min() * 0.95

    # Define plot range
    minx = (x_pos[1] + x_pos[2]) / 2.0
    maxx = (x_pos[-2] + x_pos[-3]) / 2.0
    miny = (y_pos[1] + y_pos[2]) / 2.0
    maxy = (y_pos[-2] + y_pos[-3]) / 2.0

    # Prepare the output filename
    filename = "var{:02d}_{:s}{:03d}".format(args["var"], args['cut'], icut)

    if args['movie'] > 0:
        img_file_fmt = movie_routines.setup_movie_dir(filename)
    else:
        img_file_fmt = filename+'_{:}.'+args['ext']

    # Create a plot for each time
    for itime, utime in enumerate(all_times):
        # Initialize the figure
        fig = plt.figure(constrained_layout=False,
                         tight_layout=True, figsize=(10, 8.5))

        gs1 = mpl.gridspec.GridSpec(nrows=2, ncols=2, wspace=0.0, hspace=0)
        gs = mpl.gridspec.GridSpec(nrows=2, ncols=2, wspace=0.0, left=0.0,
                                   right=0.9)
        ax = fig.add_subplot(gs1[1, :2])

        # Plot the global data set (square plot at bottom if three plots):
        con = ax.pcolor(x_pos, y_pos, all_2dim_data[itime].transpose(),
                        vmin=mini, vmax=maxi, cmap=cmap)

        # Add the winds, if desired
        if args["winds"]:
            ax.quiver(x_pos, y_pos, all_winds_x[itime].transpose(),
                      all_winds_y[itime].transpose())
            ax.set_ylim([miny, maxy])
            ax.set_xlim([minx, maxx])

        # Set the labels and aspect ratio
        ax.set_title("{:s}; {:s}: {:.2f} {:s}".format(
            utime.strftime("%d %b %Y %H:%M:%S UT"), args['cut'], z_val,
            'km' if args['cut'] == 'alt' else r'$^\circ$'))
        ax.set_xlabel(r'Latitude ($^\circ$)' if args['cut'] == 'lon'
                      else r'Longitude ($^\circ$)')
        ax.set_ylabel(r'Latitude ($^\circ$)' if args['cut'] == 'alt'
                      else r'Altitude (km)')
        if args['cut'] == 'alt':
            ax.set_aspect(1.0)

        # Set the colorbar
        cbar = fig.colorbar(con, ax=ax, shrink=0.75, pad=0.02)
        cbar.set_label(header["vars"][args["var"]], rotation=90)

        # If this is an altitude slice, add polar dials
        if args['cut'] == 'alt' and (plot_north or plot_south):
            # Set the common inputs
            shift = time_conversion.calc_time_shift(utime)
            xlabels = ['', '12', '18', '00']
            ylabels = [r'80$^\circ$', r'70$^\circ$', r'60$^\circ$',
                       r'50$^\circ$']
            xticks = np.arange(0, 2 * np.pi, np.pi / 2.0)
            yticks = np.arange(10, 50, 10)

            if plot_north:
                # Top Left Graph Northern Hemisphere
                ax2 = fig.add_subplot(gs[0, 0], projection='polar')
                rad, theta = np.meshgrid(90.0 - y_pos[mask_north],
                                         np.radians(x_pos + shift - 90.0))
                conn = ax2.pcolor(theta, rad,
                                  all_2dim_data[itime][:, mask_north],
                                  vmin=mini_north, vmax=maxi_north, cmap=cmap)
                ax2.set_xticklabels(xlabels)
                ax2.set_yticklabels(ylabels)
                ax2.grid(linestyle=':', color='black')
                ax2.set_xticks(xticks)
                ax2.set_yticks(yticks)
                fig.colorbar(conn, ax=ax2, shrink=0.5, pad=0.01)

            if plot_south:
                # Top Right Graph Southern Hemisphere
                rad, theta = np.meshgrid(90.0 + y_pos[mask_south],
                                         np.radians(x_pos + shift - 90.0))
                ax3 = fig.add_subplot(gs[0, 1], projection='polar')
                cons = ax3.pcolor(theta, rad,
                                  all_2dim_data[itime][:, mask_south],
                                  vmin=mini_south, vmax=maxi_south, cmap=cmap)
                ax3.set_xticklabels(xlabels)
                ax3.set_yticklabels(ylabels)
                ax3.grid(linestyle=':', color='black')
                ax3.set_xticks(xticks)
                ax3.set_yticks(yticks)
                fig.colorbar(cons, ax=ax3, shrink=0.5, pad=0.01)

        # Format the output filename
        fmt_input = itime if args['movie'] > 0 else utime.strftime(
            '%y%m%d_%H%M%S')
        outfile = img_file_fmt.format(fmt_input)

        # Save the output file
        fig.savefig(outfile)
        plt.close(fig)

    # Create a movie, if desired
    if args['movie'] > 0:
        movie_routines.save_movie(filename, ext=args['ext'],
                                  rate=args['movie'])

    return


# Needed to run main script as the default executable from the command line
if __name__ == '__main__':
    plot_model_results()
