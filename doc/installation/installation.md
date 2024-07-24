# Installation instructions

Before installing the Aether model, it is recommended to ensure the required
[dependencies](#dependencies) for your system are met. The next step is to
[download](#get-aether) and then [build](#build-aether) the software. Your
install can then be [tested](#test-the-executable).

## Dependencies

The Aether development team has tested several possible configurations of the
following dependencies. Note that other versions of these programs may possibly
work, however the listed versions are the recommended starting point.

If you are working on an HPC cluster, it is likely that all of these
dependencies will already be installed. Consult the documentation of your
specific system to find out how to load the software, which will likely be done
with `module`.

### MacOS
  
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
  
### Linux
  
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
[here](https://www.linode.com/docs/guides/linux-package-management-overview/#comparison-of-package-managers).

> Programs such as `conda` and `snap` are **not** recommended to install
> dependencies.

## Get Aether

Clone Aether from GitHub into a local directory of your choosing.  We currently
recommend checking out the `develop` branch, as the model is still awaiting its
first release.

```bash
cd [some/directory]
git clone git@github.com:AetherModel/Aether.git
cd Aether
git checkout develop 
```

## Build Aether

Start by creating a build directory for the necessary files.  In this example,
we show how to do so from the Aether project directory and assume you start
there.

```bash
mkdir build
cd build
```

Next, use `cmake` to create the necessary make files in the build directory

```bash
cmake ..
```

This is also the point where compilation options can be chosen. These take
the structure:

```bash
cmake -DFLAG=VALUE ..
```

Here `FLAG` is a flag name and `VALUE` is the desired value (note the `-D`).  A
more complete discussion of the available compilation flags can be found on the
[Compilation Options](build_opts.md) page.

If your default compiler isn't a GCC compiler, you will likely need to specify
the desired GCC compiler at this step using:

```bash
cmake -DCMAKE_CXX_COMPILER=<gcc or mpi executable with full path>
```

Finally, run the make command in the build directory:

```bash
make [-jN]
```

The `-j` flag tells the computer the number of processers to use when running
the code (where N is that number).  To use a single processor, just do not use
the flag.

If the make command fails, try using the `VERBOSE=1` flag.  This will include
additional output that may make it easier to debug the process.

## Test the executable

To test the model executable, exit the build diretory and run the test script:

```bash
cd ..
cp -R share/run ./run.test
cd run.test
./aether
```

The output should look like this:

```bash
> Need to NOT adjust F10.7, but that isn't included yet!!!
> Writing file : 3DALL_20110320_000000
> Writing file : 3DBFI_20110320_000000
> Wall Time : 4s (left : 1111h); Current Time : 2011 3 20 0 0 0 0 
> Wall Time : 4s (left : 23m); Current Time : 2011 3 20 0 0 10 0 
> Wall Time : 5s (left : 14m); Current Time : 2011 3 20 0 0 20 0 
> Wall Time : 5s (left : 9m); Current Time : 2011 3 20 0 0 30 0 
> Wall Time : 5s (left : 7m); Current Time : 2011 3 20 0 0 40 0 
> Wall Time : 6s (left : 7m); Current Time : 2011 3 20 0 0 50 0 
> Wall Time : 6s (left : 5m); Current Time : 2011 3 20 0 1 0 0 
> Wall Time : 6s (left : 5m); Current Time : 2011 3 20 0 1 10 0 
> Wall Time : 7s (left : 5m); Current Time : 2011 3 20 0 1 20 0 
> Wall Time : 7s (left : 4m); Current Time : 2011 3 20 0 1 30 0 
etc.
```

For more details about running Aether, see the documentation about creating and
modifying input files in the document running_aether.md.
