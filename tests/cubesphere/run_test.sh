#!/bin/sh

rm -rf run.cube

cp -R ../../share/run ./run.cube
cd run.cube
cp ../../../share/examples/aether_cubesphere.json ./aether.json
mpirun -np 6 ./aether
# plot the output
cd UA/output
../../../../../srcPython/postAether.py -rm
# [O]:
run_plot_block_model_results.py -var=O -alt=30 3DALL_20110320_001000.nc
# Tn:
run_plot_block_model_results.py -var=Temperature -alt=30 3DALL_20110320_001000.nc
# [e-]
run_plot_block_model_results.py -var=e- -alt=5 3DALL_20110320_001000.nc

