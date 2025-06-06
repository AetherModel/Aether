# Compilation Options

When compiling the code, you are able to set variables to pass to the compiler.
This can enable things like NetCDF optputs, for example. This page serves to
give a rough description of the available compilation flags and their purpose.

Depending on your system, you may also need to either declare local environment
variables or define them using this flag system.  These will show up as `cmake`
warnings.

## Table of Compilation Options

The following table lists many of the optional compilation flags. For boolean
flags, the value can be set to `ON` or `Y` to enable it. All boolean options are
disabled by default.

| Compilation Flag       | Type    | Description                              |
| ---------------------- | ------- | ---------------------------------------- |
| `USE_DOUBLE_PRECISION` | boolean | Run Aether with higher precision         |
| `CMAKE_CXX_COMPILER`   | path    | Path to gcc or MPI executable            |
| `USE_FORTRAN`          | boolean | Enables the compilation of Fortran code. |
| `USE_NETCDF`           | boolean | Write outputs to NetCDF format?          |

Recall, the flags are specified with `cmake -DFLAG=VALUE ..`, where  `FLAG` is a
flag name and `VALUE` is the desired value.

### Additional Test flags

The following flags are used for testing vaious aspects of the Aether model.
These are not recommended to be enabled unless you are involved in development.

- `TEST_INTERPOLATION`
- `TEST_COORD`
- `TEST_EXCHANGE`
  - Enables `USE_DOUBLE_PRECISION` as well.
- `TEST_GRADIENT`

> These flags are mutually exclusive; multiple *cannot* be used at the same time.
