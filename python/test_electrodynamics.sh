#!/bin/sh

# This produces 2 electrodynamics files (one in amie binary and one in netcdf):
./electrodynamics_write.py

# this tests the reading of the AMIE binary (producing amie_test.png)
./amie_read_binary.py test_ed.bin

# this tests the reading of the netcdf file (producing netcdf_test.png)
./electrodynamics_read.py test_ed.nc


