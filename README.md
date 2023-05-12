# Aether
This is the home of the Aether model of the thermosphere and ionosphere

The Aether model has been developed using gnu c++ (version 10). If
you are using this, hopefully it will just work out of the box. We have 
been doing development of Aether on Mac OSX, and in Ubuntu Linux.  We have 
also used the Windows Subsystem for Linux, Ubuntu distribution, which 
works similarly to the native Linux distribution.

If you are a student and don't know how to work with a large code base
(i.e., multiple source codes in multiple directories), you may consider
starting with the doc/student.md file.

## Dependencies:

1. Aether uses [CMake](https://cmake.org/) instead of make. If you don't have
   it installed,

For MacOS [homebrew](https://formulae.brew.sh/formula/cmake):
```bash
sudo brew install cmake
```

This can also be installed using [macports](https://www.macports.org/)
```bash
sudo port install cmake
```

For Ubuntu/Debian Linux:
```bash
sudo apt install cmake
```

2. Aether uses the netcdf library (netcdf-cxx4). As above, netCDF can be
   installed using a package manager

On Mac, this is can be awkward, depending on which c++ compiler you are using.
Since there is one that essentially comes with Mac OSX, called clang, the
default compiler is often this.  Much of the other software is not built with
this, so you need to switch compilers, which can be challenging.  

On Mac, if you want the clang compiled version of netcdf, then:
```bash
sudo port install netcdf-cxx4
```

If you want the gcc version of netcdf (recommended), then:
```bash
sudo port install netcdf-cxx4 +gcc10
```

On Ubuntu, gcc is the default compiler, it seems like you can probably just do:
```bash
sudo apt-get install libnetcdf-dev
sudo apt install libnetcdf-c++4-dev
```

3. Aether uses the nlohman json package for reading and writing json files.

On Ubuntu:

```bash
sudo apt-get install -y nlohmann-json-dev
```

On Mac:

```bash
sudo port install nlohmann-json 
```

4. The armadillo headers need to be installed. Simplistically, Armadillo is a
system that allows matrix math to be done in C++ easily. We mostly use it for
doing math with matrices (like multiplication, addition, etc.), but it is much
more powerful than this.  You will notice that there are not many 3D loops in
Aether, which is due to Armadillo.  To make this all fast, it is best to install
the lapack abd blas libraries too.

On Ubuntu:

```bash
sudo apt-get install liblapack-dev
sudo apt-get install libblas-dev
sudo apt-get install libboost-dev
sudo apt-get install libarmadillo-dev
```

On Mac:

```bash
sudo port install lapack
sudo port install OpenBLAS
sudo port install boost
sudo port install armadillo
 ```

## Quick Start:

These are unix commands, assuming that you have access to a unix/linux
terminal. This has been tested on a MacBook Pro and Ubuntu.

```bash
git clone https://github.com/AetherModel/Aether
```

```bash
cd Aether
```

```bash
git checkout develop
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

To compile Aether with double precision:
```bash
mkdir build
cd build
cmake -DUSE_DOUBLE_PRECISION=Y ..
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

There are essentially two input files that specify the settings in the code.
When you are in a run directory, they are:

1. UA/inputs/defaults.json.  These set the default inputs for the run
and should not be modified at all.  You can look at these and copy the
settings that you want to change to this file:

2. aether.json.  This file can and should be modified to direct the
code to run the way that you would like.  You can copy settings from
the default.json file and then modify them here. This will be covered
in the reference manual, once we have written one.

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
