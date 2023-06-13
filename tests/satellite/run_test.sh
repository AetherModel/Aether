#!/bin/bash

# Stop on errors
# See https://vaneyckt.io/posts/safer_bash_scripts_with_set_euxo_pipefail/
set -Eeuo pipefail

echo "This test takes about 8 minutes to run"

# Initialize the environment
rm -rf run.test
cp -r ../../share/run run.sphere

# Test sphere
cd run.sphere
cp ../sat_20110320.csv UA/inputs/
cp ../aether_sat_sphere.json aether.json
mpirun -np 4 ./aether
../../../srcPython/satellite_combine.py
../../../srcPython/satellite_test.py
mv Satellite_log.png ../Satellite_sphere_log.png
cd ..

# Test cubesphere
rm -rf run.test
cp -r ../../share/run run.cube
cd run.cube
cp ../aether_sat_cube.json aether.json
cp ../sat_20110320.csv UA/inputs/
mpirun -np 6 ./aether
../../../srcPython/satellite_combine.py
../../../srcPython/satellite_test.py
mv Satellite_log.png ../Satellite_cubesphere_log.png
cd ..

# Clear environment
#rm -rf run.test
echo "The difference between the satellite file and log file is shown in two png files."
