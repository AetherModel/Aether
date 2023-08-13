#!/bin/sh

# simple test to make sure that the hdf5 output and plotting works ok!

rm -rf run.hdf5
 
cp -r ../../share/run ./run.hdf5
cd run.hdf5
mpirun -np 4 ./aether

cd UA/output

../../../../../srcPython/postAether.py -rm -hdf5
~/bin/run_plot_block_model_results.py -var=Temperature -alt=30 3DALL_20110320_001000.hdf5

ls

