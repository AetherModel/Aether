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

2. The armadillo include files need to be placed somewhere that
can be accessed by the Makefile.

## Quick Start:

These are unix commands, assuming that you have access to a unix/linux
terminal. This has been tested on a MacBook Pro.

1. git clone https://github.com/AetherModel/Aether

2. cd Aether

3. git checkout develop

4. Copy example Makefile.OSX/ubuntu and modify src/Makefile.OS to
adjust compiler options and to point to the correct netcdf library and
armadillo include file locations.  We will fix this so that there is a
configuration script at some point.

5. make 

6. make rundir

7. cd run

8. ./aether.exe (should run 1 hour of simulation with no issues).

9. Make some plots:

../python/plot_model_results.py -var=24 3DALL_20110320_010000.nc -alt=110

../python/plot_model_results.py -var=14 3DALL_20110320_010000.nc -alt=300

10. compare png files to ../inputs/*.png to see if they are similar.

## Code Manual:

To create the code documentation manual, download Doxygen for your operating
system and run:

```bash
cd doc
doxygen Doxyfile
```
