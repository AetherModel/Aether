#!/bin/bash

# Stop on errors
# See https://vaneyckt.io/posts/safer_bash_scripts_with_set_euxo_pipefail/
set -Eeuo pipefail

echo "This test takes about 8 minutes to run"

# Initialize the environment
rm -rf run.test

# Test sphere
cp -r ../../share/run run.sphere
cd run.sphere
# Put satellite file and json file
cp ../sat_20110320.csv UA/inputs/
cp ../aether_sat_sphere.json aether.json
# Run aether
mpirun -np 4 ./aether
# Generate graph
../../../srcPython/satellite_test.py
# Put the graph to parent dir
mv Satellite_log.png ../Satellite_sphere_log.png
cd ..

# Test cubesphere
cp -r ../../share/run run.cube
cd run.cube
cp ../aether_sat_cube.json aether.json
cp ../sat_20110320.csv UA/inputs/
mpirun -np 6 ./aether
../../../srcPython/satellite_test.py
mv Satellite_log.png ../Satellite_cubesphere_log.png
cd ..

# Clear environment
# rm -rf run.test
echo "The difference between the satellite file and log file is shown in two png files."
