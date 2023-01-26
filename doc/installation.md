# Installation instructions

## Dependencies

Before beginning Aether installation, make sure the following dependencies are
available on your system.  We recommend using a package manager, such as
MacPorts, Homebrew (Mac OS X), apt-get (Ubuntu), or Chocolatey (Windows).

| Mac OS X                       |
|--------------------------------|
| Dependency    | Tested version |
|---------------|----------------|
| armadillo     | 11.4           |
| boost         | 1.76           |
| cmake         | 2.24           |
| gcc           | 10, 11, 12     |
| netcdf-cxx4   | 4.9            |
| nlohmann-json | 3.11           |
| OpenBLAS      | 0.3            |
| openMPI       | same as gcc    |



| Ubuntu                              |
|-------------------------------------|
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
| nlohmann-json-dev  |                |
| openMPI            | same as gcc    |

# Get Aether

Clone Aether from GitHub into a local directory of your choosing.  We currently
recommend checking out the `develop` branch, as the model is still awaiting its
first release.

# Build Aether

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
cmake .. -DFLAG=VALUE
```

where `FLAG` is a flag name and `VALUE` is the desired value.  A table of model
options is included below. Depending on your system, you may also need to
either declare local environment variables or define them using this flag
system.  These will show up as cmake warnings.

| Aether Flag          | Value | Description                       |
|----------------------|-------|-----------------------------------|
| USE_NETDCDF          | Y     | Allow NetCDF input/output files   |
| USE_DOUBLE_PRECISION | Y     | Run Aether with higher precission |

Finally, run the make command in the build directory:

```bash
make -j
```

The `-j` flag DOES THIS.

# Test the executable

To test the model executable, exit the build diretory and run the test script:

```bash
cd ..
cp -R share/run ./run.test
cd run.test
./aether
```

The output should LOOK LIKE THIS

For more details about running Aether, see the documentation about creating and
modifying input files HERE.