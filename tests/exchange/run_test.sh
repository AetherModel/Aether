#!/bin/sh

rm -rf build
mkdir build
cd build
cmake -DTEST_EXCHANGE:BOOL=TRUE ../../..
make -j4
cd ..

rm -rf run.w_cubesphere
cp -R ../../share/run ./run.w_cubesphere
cd run.w_cubesphere
cp ../aether.json.cubesphere ./aether.json
mpirun -np 6 ./aether
cd UA/output
../../../../../srcPython/postAether.py -rm

~/bin/run_plot_block_model_results.py -var=lat 3DALL_20110320_000000.nc
~/bin/run_plot_block_model_results.py -var=lon 3DALL_20110320_000000.nc




