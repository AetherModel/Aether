# Aethers
This is the home of the Aether model of the thermosphere and ionosphere

The Aether model has been developed using gnu c++ (version 9.3.0). If
you are using this, hopefully it will just work out of the box.

## Dependencies:

1. Aether uses the netcdf library (netcdf-cxx4), but we wrote a
binary output file also, so we can choose which one we want.

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
cd ..
cp -R share/run ./run.test
cd run.test
./aether
```

Output files are in UA/output.

We are working on aetherpy to make plots.

Compare png files to tests/outputs_pngs/*.png to see if they are similar.

## Code Manual:

To create the code documentation manual, download Doxygen for your operating
system and run:

```bash
cd doc
doxygen Doxyfile
```
