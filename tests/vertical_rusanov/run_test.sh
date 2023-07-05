#!/bin/sh


rm -rm ./run.test_*

cp -R ../../share/run ./run.test_rusanov
cd run.test_rusanov
cp ../aether.json.rusanov ./aether.json
mpirun -np 4 ./aether

# post process and plot:
cd UA/output
../../../../../srcPython/postAether.py -rm
# not sure where plotting code is located....
~/Software/aetherpy/vAaron/aetherpy/bin/run_plot_block_model_results.py 3DALL_20110320_*.nc -var="Vertical Wind" -alt=30

# back into main test directory:
cd ../../..
