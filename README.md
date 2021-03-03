# Aethers
This is the home of the Aether model of the thermosphere and ionosphere

The Aether model has been developed using gnu c++ (version 9.3.0). If
you are using this, hopefully it will just work out of the box.

## Dependencies:

1. Aether uses the netcdf library (netcdf-cxx4). We will eventually
make a configuration file that will check to see if you have this
installed, but right now it is hardcoded to be in
/opt/local/lib. Sorry. Also, the python code provided for
visualization uses the netcdf library (netCDF4).

## Quick Start:

These are unix commands, assuming that you have access to a unix/linux
terminal. This has been tested on a MacBook Pro.

1. git clone https://github.com/AetherModel/Aether

2. cd Aether

3. git checkout develop

4. Modify src/Makefile to point to the correct netcdf library (in LINK
command near the bottom - may need to change the location and the name
of the library, as well as the FLAGS variable at the top). Also, if
you are using a different c++ compiler, you will need to change that
at this time (top two variables).  We will fix this so that there is a
configuration script at some point.

5. make (may throw warning about ../lib directory not found.)

6. make rundir

7. cd run

8. ./aether.exe (should run 10 minutes with no issues).

9. python ../python/read_python.py

10. compare test.png to ../inputs/test.png to see if they are similar.


