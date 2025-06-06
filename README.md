# Aether

This is the home of the Aether model of the thermosphere and ionosphere.

The Aether model has been developed using GNU C++ (versions 9, 10, 11). If you
are using this, hopefully it will just work out of the box. We have been doing
development of Aether on MacOS and Ubuntu Linux.  We have also used the Windows
Subsystem for Linux, Ubuntu distribution, which works similarly to the native
Linux distribution.

> If you are a student and don't know how to work with a large code base (i.e.,
multiple source codes in multiple directories), you may consider starting within
the [Students](doc/student.md) page.

All other users may wish to continue with the short installation guide below. If
you run into issues or have questions, first consult the in-progress
[documentation](doc/README.md). More complete documentation is actively being
developed.

## Contents

- [Aether](#aether)
  - [Contents](#contents)
  - [Quick Start](#quick-start)
    - [Dependencies](#dependencies)
    - [Getting the Code](#getting-the-code)
    - [Compiling \& Running](#compiling--running)
  - [Code Manual](#code-manual)
  - [Further Documentation](#further-documentation)

## Quick Start

The following guide serves to demonstrate how to download, install, and then run
the Aether model. More details on each topic will be linked, and some more
detail is available [here](doc/README.md).

### Dependencies

On MacOS, installing some of the dependencies can be awkward, depending on which
C++ compiler you are using. Since there is one that essentially comes with
MacOS, called `clang`, the default compiler is often this.  Much of the other
software is not built with this, so you need to switch compilers, which can be
challenging.

While the development team has tried to remove as many of the dependencies as
possible to reduce this issues, it is still a good idea to try to get a
different compiler to work.  If you want to use NetCDF files or the Fortran
compiled options (like MSIS), you will probably have to do this.

A package manager is recommended to install the dependencies of the Aether
model. More details on the use of package managers can be found on [the
dependencies page](doc/installation/dependencies.md). For MacOS and Ubuntu Linux, the validated configurations are:

<details>
  <summary> MacOS </summary>
  
  | Dependency    | Tested version |
  |---------------|----------------|
  | armadillo     | 11.4           |
  | boost         | 1.76           |
  | cmake         | 2.24           |
  | gcc           | 10, 11         |
  | netcdf        | 4.9            |
  | netcdf-cxx4   | 4.3            |
  | nlohmann-json | 3.11           |
  | OpenBLAS      | 0.3            |
  | mpich         | 4.1            |

  MacOS has two predominant package managers: [Homebrew](https://brew.sh) and
  [MacPorts](https://www.macports.org/). Either will work.
  
</details>

<details>
  <summary> Linux </summary>
  
  | Dependency         | Tested version |
  |--------------------|----------------|
  | cmake              | 2.24           |
  | gcc                | 10, 11, 12     |
  | libarbadillo-dev   |                |
  | libblas-dev        |                |
  | libboost-dev       |                |
  | liblapack-dev      |                |
  | libnetcdf-dev      |                |
  | libnetcdf-c++4-dev |                |
  | libopenmpi-dev     |                |
  | nlohmann-json3-dev |                |
  | openmpi-bin        | same as gcc    |

  The specific package manager to use depends on which distribution of Linux
  Aether is installed on. More details on Linux package managers can be found
  [here, for
  example](https://www.linode.com/docs/guides/linux-package-management-overview/#comparison-of-package-managers).

</details>

### Getting the Code

This has been tested on a MacBook Pro and Ubuntu.

```bash
git clone git@github.com:AetherModel/Aether.git
cd Aether
git checkout develop
```

> ***NOTES:***
>
> - (This is assuming that you are installing the root version of Aether and not
a forked version.  If you are using a forked version, replace the
"`AetherModel`" with the appropriate location of the fork.)
> - GitHub now recommends cloning repositories via `ssh` instead of `https`. For
development work, this is the way to go. If the above `git clone` command fails
and you are not planning on contributing to Aether development, the first line
can be replaced with "`git clone https://github.com/AetherModel/Aether`". For
help setting up your machine to use `ssh` to connect to Github, see [this
page](https://docs.github.com/en/authentication/connecting-to-github-with-ssh)

### Compiling & Running

To compile Aether:

```bash
mkdir build
cd build
cmake ..
make [-j4]
```

The `-j4` is optional and uses 4 processors, but you could use 2, 4, 8, or
whatever.

More compilation options can be found on the [Compiling
Aether](doc/installation/build_opts.md) page.

Once you have compiled, you can install Aether with an example run directory
structure like this:

```bash
cd ..
cp -R share/run ./run.test
cd run.test
./aether
```

There are essentially two input files that specify the settings in the code.
When you are in a run directory, they are:

1. UA/inputs/defaults.json.  These set the default inputs for the run and should
not be modified at all.  You can look at these and copy the settings that you
want to change to this file.

2. aether.json.  This file can and should be modified to direct the code to run
the way that you would like.  You can copy settings from the default.json file
and then modify them here. This will be covered in the reference manual, once we
have written one.

You can check to make sure that these are valid json files (not checking the
content, though) with:

```bash
cd run.test
python -m json.tool aether.json
python -m json.tool UA/inputs/defaults.json
```

At this time, there is no checker to see if all of the settings in each of the
inputs files are actually valid and Aether understands them.

Output files are in UA/output.

We are working on aetherpy to make plots.

Compare png files to tests/outputs_pngs/*.png to see if they are similar.

## Code Manual

To create the code documentation manual, download Doxygen for your operating
system and run:

```bash
cd doc
doxygen Doxyfile
```

## Further Documentation

A more complete documentation is in
[development](https://github.com/AetherModel/AetherDocumentation), which can be
browsed [here](https://aetherdocumentation.rtfd.io/).

The collection of Markdown files within the [doc](doc/README.md) folder of this
repository will provide more information and, at the time of writing, is more up
to date than the official documentation. It is recommended to start there.
