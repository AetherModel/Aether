#!/bin/sh

# remove old directories:
rm -rf build run.*

# make code WITH Fortran:
mkdir build
cd build
cmake -DUSE_FORTRAN=Y ../../..
make -j4
cd ..
cp -R ../../share/run ./run.ie
cd run.ie

# run IE with spherical
cp ../aether.json.sphere.ie_test ./aether.json
mpirun -np 4 ./aether
cd UA/output
../../../../../srcPython/postAether.py -rm
# This assumes aetherpy is installed and the plotter is in the bin directory:
~/bin/run_plot_block_model_results.py -var=Temperature -alt=20 3DALL_20110320_001500.nc
~/bin/run_plot_block_model_results.py -var=Potential -alt=20 3DALL_20110320_001500.nc
# into UA directory
cd .. 
mv output output.ie_sphere
mkdir output
# into run directory
cd .. 

# run with MSIS (this should be different from the other three runs!):
cp ../aether.json.cube.ie_test ./aether.json
mpirun -np 6 ./aether
cd UA/output
../../../../../srcPython/postAether.py -rm
# This assumes aetherpy is installed and the plotter is in the bin directory:
~/bin/run_plot_block_model_results.py -var=Temperature -alt=20 3DALL_20110320_001500.nc
~/bin/run_plot_block_model_results.py -var=Potential -alt=20 3DALL_20110320_001500.nc
# into UA directory
cd .. 
mv output output.ie_sphere
mkdir output
# into run directory
cd .. 



