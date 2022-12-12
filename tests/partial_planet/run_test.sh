#!/bin/sh

export TEST=partial
export DIR=./run.$TEST

rm -rf $DIR

cp -R ../../share/run $DIR
cd $DIR
rm -f ./aether.json
cp ../../../share/examples/aether_$TEST.json ./aether.json
mpirun -np 4 ./aether
# plot the output
cd UA/output
../../../../../srcPython/postAether.py -rm

# not going to make any other plots, due to z - alt incompatability
# [O]:
#run_plot_block_model_results.py -var=O -alt=30 3DALL_20110320_001000.nc
# Tn:
#run_plot_block_model_results.py -var=Temperature -alt=30 3DALL_20110320_001000.nc
# [e-]
#run_plot_block_model_results.py -var=e- -alt=5 3DALL_20110320_001000.nc
