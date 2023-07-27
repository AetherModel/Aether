#!/bin/sh

# remove old directories:
rm -rf run.*

cp -R ../../share/run ./run.tests
cd run.tests

# ------------------------------------------------
# run default:
cp ../aether.json.v0 ./aether.json
mpirun -np 4 ./aether
cd UA/output
../../../../../srcPython/postAether.py -rm
# This assumes aetherpy is installed and the plotter is in the bin directory:
~/bin/run_plot_block_model_results.py -var=Temperature -alt=10 3DNEU_20110320_001000.nc
~/bin/run_plot_block_model_results.py -var=O -alt=20 3DNEU_20110320_001000.nc
~/bin/run_plot_block_model_results.py -var="O Vertical Wind" -alt=5 3DNEU_20110320_001000.nc
# into UA directory
cd .. 
mv output output.v0
mkdir output
# into run directory
cd ..

# ------------------------------------------------
# run momentum = F:
cp ../aether.json.v1.momentum_f ./aether.json
mpirun -np 4 ./aether
cd UA/output
../../../../../srcPython/postAether.py -rm
# This assumes aetherpy is installed and the plotter is in the bin directory:
~/bin/run_plot_block_model_results.py -var=Temperature -alt=10 3DNEU_20110320_001000.nc
~/bin/run_plot_block_model_results.py -var=O -alt=20 3DNEU_20110320_001000.nc
~/bin/run_plot_block_model_results.py -var="O Vertical Wind" -alt=5 3DNEU_20110320_001000.nc
# into UA directory
cd .. 
mv output output.v1
mkdir output
# into run directory
cd ..

# ------------------------------------------------
# run energy = F:
cp ../aether.json.v2.energy_f ./aether.json
mpirun -np 4 ./aether
cd UA/output
../../../../../srcPython/postAether.py -rm
# This assumes aetherpy is installed and the plotter is in the bin directory:
~/bin/run_plot_block_model_results.py -var=Temperature -alt=10 3DNEU_20110320_001000.nc
~/bin/run_plot_block_model_results.py -var=O -alt=20 3DNEU_20110320_001000.nc
~/bin/run_plot_block_model_results.py -var="O Vertical Wind" -alt=5 3DNEU_20110320_001000.nc
# into UA directory
cd .. 
mv output output.v2
mkdir output
# into run directory
cd ..

# ------------------------------------------------
# run Uniform = true:
cp ../aether.json.v3.uniform ./aether.json
mpirun -np 4 ./aether
cd UA/output
../../../../../srcPython/postAether.py -rm
# This assumes aetherpy is installed and the plotter is in the bin directory:
~/bin/run_plot_block_model_results.py -var=Temperature -alt=10 3DNEU_20110320_001000.nc
~/bin/run_plot_block_model_results.py -var=O -alt=20 3DNEU_20110320_001000.nc
~/bin/run_plot_block_model_results.py -var="O Vertical Wind" -alt=5 3DNEU_20110320_001000.nc
# into UA directory
cd .. 
mv output output.v3
mkdir output
# into run directory
cd ..


