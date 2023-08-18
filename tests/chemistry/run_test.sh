#!/bin/sh

# remove old directories:
rm -rf run.*

# run without perturbing:
cp -R ../../share/run ./run.no_perturb
cd run.no_perturb

cp ../aether.json.no_perturb ./aether.json
mpirun -np 1 ./aether
cd UA/output
../../../../../srcPython/postAether.py -rm
# This assumes aetherpy is installed and the plotter is in the bin directory:
../../../../../aetherpy/bin/aether_plot_simple.py -var=Temperature -alt=20 3DALL_20110320_003000.nc
../../../../../aetherpy/bin/aether_plot_simple.py -var=O2+ -alt=10 3DALL_20110320_003000.nc
# into test directory
cd ../../.. 

# run with perturbing (4 ensemble members)
cp -R ../../share/run ./run.w_perturb
cd run.w_perturb

cp ../aether.json.w_perturb ./aether.json
mpirun -np 4 ./aether
cd UA/output
../../../../../srcPython/postAether.py
# This assumes aetherpy is installed and the plotter is in the bin directory:
../../../../../aetherpy/bin/aether_plot_simple.py -var=Temperature -alt=20 3DALL_20110320_003000_std.nc
../../../../../aetherpy/bin/aether_plot_simple.py -var=Temperature -alt=20 3DALL_20110320_003000_mean.nc
../../../../../aetherpy/bin/aether_plot_simple.py -var=O2+ -alt=10 3DALL_20110320_003000_std.nc
../../../../../aetherpy/bin/aether_plot_simple.py -var=O2+ -alt=10 3DALL_20110320_003000_mean.nc
# into test directory
cd ../../..


