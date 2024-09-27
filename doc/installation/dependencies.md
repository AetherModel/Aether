# Installing Dependencies

Instructions for MacOS are written for use with MacPorts, as the development
team has mostly used MacPorts and not Homebrew. Commands for the two package
managers are nearly identical, although Homebrew may be preferred if you do not
have root access on your machine. For example, the following two commands will
both install `cmake`. Translating the following guide from MacPorts to Homebrew
is left to the user.

```bash
sudo port install cmake
brew install cmake
```

> To check if any of the following dependencies are already installed on your
system, you will need access to the command line (terminal). To check if (and
where) `cmake` is installed, for example, the command `which cmake` can be run.
If a path is printed, `cmake` is installed. To check the version, run `cmake
--version`.

The layout of this page is as follows:

- [Installing Dependencies](#installing-dependencies)
  - [Install gcc](#install-gcc)
  - [Install cmake](#install-cmake)
  - [Install JSON libraries](#install-json-libraries)
  - [Install Armadillo (and boost)](#install-armadillo-and-boost)
  - [Install NetCDF (optional)](#install-netcdf-optional)

## Install gcc

This comes installed by default on Ubuntu. On MacOS this can be installed, for
example, using:

```bash
sudo port install gcc11
```

> As development began, gcc11 was the latest version; there are newer versions
> of `gcc` available now (latest version is gcc14), which have not yet been
> validated.

## Install cmake

Aether uses [CMake](https://cmake.org/) instead of `make`. If you don't have it
installed, you need it.

For MacOS, this can be installed with:

```bash
sudo port install cmake
```

For Ubuntu/Debian Linux:

```bash
sudo apt install cmake
```

This can be done on RedHat using yum also.

## Install JSON libraries

Aether uses the nlohman json package for reading and writing json files.

On Ubuntu:

```bash
sudo apt-get install -y nlohmann-json3-dev
```

On Mac:

```bash
sudo port install nlohmann-json 
```

## Install Armadillo (and boost)

Simplistically, Armadillo is a system that allows matrix math to be done in C++
easily. We mostly use it for doing math with matrices (like multiplication,
addition, etc.), but it is much more powerful than this.  You will notice that
there are not many 3D loops in Aether, which is due to Armadillo.  To make this
all fast, it is best to install the `lapack` and `blas` libraries too.

On Ubuntu:

```bash
sudo apt-get install liblapack-dev
sudo apt-get install libblas-dev
sudo apt-get install libboost-dev
sudo apt-get install libarmadillo-dev
sudo apt-get install openmpi-bin libopenmpi-dev
```

On Mac:

```bash
sudo port install lapack
sudo port install OpenBLAS
sudo port install boost
sudo port install armadillo
sudo port install openmpi-bin libopenmpi-dev
 ```

## Install NetCDF (optional)

We have removed the strict dependency for NetCDF, but a lot of codes used
NetCDF, so it doesn't hurt to try to install the libraries. Aether uses the
NetCDF library (netcdf-cxx4). As above, NetCDF can be installed using a package
manager.

On Mac, if you want the clang compiled version of netcdf, then:

```bash
sudo port install netcdf-cxx4
```

If you want the gcc version of netcdf, then:

```bash
sudo port install netcdf-cxx4 +gcc10
```

On Ubuntu, gcc is the default compiler, it seems like you can probably just do:

```bash
sudo apt-get install libnetcdf-dev
sudo apt install libnetcdf-c++4-dev
```
