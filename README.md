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

```bash
git clone https://github.com/AetherModel/Aether
```

```bash
cd Aether
```

```bash
git checkout develop
```

Make sure you have [CMake](https://cmake.org/) installed. If you don't:

For MacOS [homebrew](https://formulae.brew.sh/formula/cmake):
```bash
brew install cmake
```

For Ubuntu/Debian Linux:
```bash
sudo apt install cmake
```

To compile Aether:
```bash
mkdir build
cd build
cmake ..
make -j
```

You should have a file called `aether` by now. To run an example you must set
up your directory like this:

```bash
cmake -DMAKE_RUNDIR=Y ..
cd ../run
```

Make some plots:

```bash
../python/plot_model_results.py -var=24 3DALL_20110320_010000.nc -alt=110

../python/plot_model_results.py -var=14 3DALL_20110320_010000.nc -alt=300
```

Compare png files to ../inputs/*.png to see if they are similar.


