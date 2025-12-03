#!/bin/sh

PLOTTER=~/Software/PyITM/bin/plot_alt_wpolar.py

rm -rf ./run.test_*

# ----------------------------------------------------------------------
# run the idealized acheron test:

cp -R ../../share/run ./run.test_acheron_ideal
cd run.test_acheron_ideal
cp ../aether.json.ideal ./aether.json
mpirun -np 4 ./aether

# post process and plot:
cd UA/output
~/bin/postAether.py -rm

# This assumes pyitm is installed and the plotter is in the bin directory:
${PLOTTER} -var=Tn -alt=1300 3DALG_20110320_001000.nc
${PLOTTER} -var=O+ -alt=1100 3DALM_20110320_001000.nc
# into test directory
cd ../../..

# ----------------------------------------------------------------------
# run the acheron test:

cp -R ../../share/run ./run.test_acheron
cd run.test_acheron
cp ../aether.json.acheron ./aether.json
mpirun -np 4 ./aether

# post process and plot:
cd UA/output
~/bin/postAether.py -rm

# This assumes pyitm is installed and the plotter is in the bin directory:
${PLOTTER} -var=Tn -alt=1300 3DALG_20110320_001000.nc
${PLOTTER} -var=O+ -alt=1100 3DALM_20110320_001000.nc
# into test directory
cd ../../..

