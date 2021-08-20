# Aethers
This is the home of the Aether model of the thermosphere and ionosphere

The Aether model has been developed using gnu c++ (version 9.3.0). If
you are using this, hopefully it will just work out of the box.

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

To compile Aether with NetCDF:
```bash
mkdir build
cd build
cmake -DUSE_NETCDF=Y ..
make -j
```

Once you have compiled you can install Aether with an example run directory
structure like this:

```bash
make install
cd ../run
```

If you want a different run directory you could run `cmake` like this:
```bash
cmake -DRUN_DIR=/my/run/dir ..
```

Make some plots:

```bash
../python/plot_model_results.py -var=24 3DALL_20110320_010000.nc -alt=110

../python/plot_model_results.py -var=14 3DALL_20110320_010000.nc -alt=300
```

Compare png files to ../inputs/*.png to see if they are similar.


