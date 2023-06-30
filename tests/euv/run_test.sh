#!/bin/sh

# remove old directories:
rm -rf run.*

# run 59 wavelengths, neuvac
cp -R ../../share/run ./run.59_neuvac_nopei
cd run.59_neuvac_nopei
cp ../aether.json.59_neuvac_nopei ./aether.json
mpirun -np 4 ./aether
cd UA/output
../../../../../srcPython/postAether.py -rm
# This assumes aetherpy is installed and the plotter is in the bin directory:
~/bin/run_plot_block_model_results.py -var=Temperature -alt=20 3DALL_20110320_003000.nc
~/bin/run_plot_block_model_results.py -var=O2+ -alt=3 3DALL_20110320_003000.nc
# into test directory
cd ../../..

# run 59 wavelengths, neuvac
cp -R ../../share/run ./run.59_neuvac
cd run.59_neuvac
cp ../aether.json.59_neuvac ./aether.json
mpirun -np 4 ./aether
cd UA/output
../../../../../srcPython/postAether.py -rm
# This assumes aetherpy is installed and the plotter is in the bin directory:
~/bin/run_plot_block_model_results.py -var=Temperature -alt=20 3DALL_20110320_003000.nc
~/bin/run_plot_block_model_results.py -var=O2+ -alt=3 3DALL_20110320_003000.nc
# into test directory
cd ../../..

# run 59 wavelengths, euvac
cp -R ../../share/run ./run.59_euvac
cd run.59_euvac
cp ../aether.json.59_euvac ./aether.json
mpirun -np 4 ./aether
cd UA/output
../../../../../srcPython/postAether.py -rm
# This assumes aetherpy is installed and the plotter is in the bin directory:
~/bin/run_plot_block_model_results.py -var=Temperature -alt=20 3DALL_20110320_003000.nc
~/bin/run_plot_block_model_results.py -var=O2+ -alt=3 3DALL_20110320_003000.nc
# into test directory
cd ../../..

# run 37 wavelengths, euvac
cp -R ../../share/run ./run.37_euvac
cd run.37_euvac
cp ../aether.json.37_euvac ./aether.json
mpirun -np 4 ./aether
cd UA/output
../../../../../srcPython/postAether.py -rm
# This assumes aetherpy is installed and the plotter is in the bin directory:
~/bin/run_plot_block_model_results.py -var=Temperature -alt=20 3DALL_20110320_003000.nc
~/bin/run_plot_block_model_results.py -var=O2+ -alt=3 3DALL_20110320_003000.nc
# into test directory
cd ../../.. 



