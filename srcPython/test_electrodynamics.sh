#!/bin/sh

## test the electrodynamics.
#./electrodynamics_write.py -startdate=20110320 -enddate=20110321 -dt=60 -real -outfile=b20110320n_omni.bin
#./amie_read_binary.py -step=1 b20110320n_omni.bin
#
#./electrodynamics_write.py -startdate=20110320 -enddate=20110321 -dt=60 -real -south -outfile=b20110320s_omni.bin
#./amie_read_binary.py -step=1 b20110320s_omni.bin

# make electrodynamics files
./electrodynamics_write.py -startdate=20130317 -enddate=20130318 -dt=1 -real -ions -outfile=b20130317n_omni_ions.bin
./electrodynamics_write.py -startdate=20130317 -enddate=20130318 -dt=1 -real -ions -outfile=b20130317s_omni_ions.bin -south
omni_download_write_swmf.py 20130317 20130318

# plot:
./amie_read_binary.py -step=60 b20130317n_omni_ions.bin
./amie_read_binary.py -step=60 b20130317s_omni_ions.bin

#./electrodynamics_write.py -startdate=20130318 -enddate=20130319 -dt=1 -real -outfile=b20130318n_omni.bin
#./amie_read_binary.py b20130318n_omni.bin
#
#./electrodynamics_write.py -startdate=20130318 -enddate=20130319 -dt=1 -real -outfile=b20130318s_omni.bin -south
#./amie_read_binary.py b20130318s_omni.bin


## This produces 2 electrodynamics files (one in amie binary and one in netcdf):
#./electrodynamics_write.py
#
## this tests the reading of the AMIE binary (producing amie_test.png)
#./amie_read_binary.py test_ed.bin
#
## this tests the reading of the netcdf file (producing netcdf_test.png)
#./electrodynamics_read.py test_ed.nc


