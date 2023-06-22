#!/bin/sh

# remove old directories:
rm -rf build run.*

# make code WITH Fortran:
mkdir build
cd build
cmake -DUSE_FORTRAN=Y ../../..
make -j4
cd ..
cp -R ../../share/run ./run.w_fortran
cd run.w_fortran

# run without MSIS:
cp ../aether.json.wo_msis ./aether.json
mpirun -np 6 ./aether
cd UA/output
../../../../../srcPython/postAether.py -rm
# This assumes aetherpy is installed and the plotter is in the bin directory:
~/bin/run_plot_block_model_results.py -var=Temperature -alt=20 3DALL_20110320_001500.nc
~/bin/run_plot_block_model_results.py -var=O -alt=20 3DALL_20110320_001500.nc
# into UA directory
cd .. 
mv output output.wo_msis
mkdir output
# into run directory
cd .. 

# run with MSIS (this should be different from the other three runs!):
cp ../aether.json.w_msis ./aether.json
mpirun -np 6 ./aether
cd UA/output
../../../../../srcPython/postAether.py -rm
# This assumes aetherpy is installed and the plotter is in the bin directory:
~/bin/run_plot_block_model_results.py -var=Temperature -alt=20 3DALL_20110320_001500.nc
~/bin/run_plot_block_model_results.py -var=O -alt=20 3DALL_20110320_001500.nc
# into UA directory
cd .. 
mv output output.w_msis
mkdir output
# into run directory
cd .. 

# back into the test directory
cd .. 

# make code WITHOUT Fortran:
rm -rf build
mkdir build
cd build
cmake ../../..
make -j4
# back into the test directory
cd ..  
cp -R ../../share/run ./run.wo_fortran
cd run.wo_fortran

# run without MSIS:
cp ../aether.json.wo_msis ./aether.json
mpirun -np 6 ./aether
cd UA/output
../../../../../srcPython/postAether.py -rm
# This assumes aetherpy is installed and the plotter is in the bin directory:
~/bin/run_plot_block_model_results.py -var=Temperature -alt=20 3DALL_20110320_001500.nc
~/bin/run_plot_block_model_results.py -var=O -alt=20 3DALL_20110320_001500.nc
# into UA directory
cd .. 
mv output output.wo_msis
mkdir output
# into run directory
cd .. 

# run with MSIS (this should not work, but should be identical to above):
cp ../aether.json.w_msis ./aether.json
mpirun -np 6 ./aether
cd UA/output
../../../../../srcPython/postAether.py -rm
# This assumes aetherpy is installed and the plotter is in the bin directory:
~/bin/run_plot_block_model_results.py -var=Temperature -alt=20 3DALL_20110320_001500.nc
~/bin/run_plot_block_model_results.py -var=O -alt=20 3DALL_20110320_001500.nc
cd ..
mv output output.w_msis
mkdir output
cd ..


