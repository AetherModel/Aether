#!/bin/sh

rm -rf ./run.test_*

# ----------------------------------------------------------------------
# run the rusanov test:

cp -R ../../share/run ./run.test_rusanov
cd run.test_rusanov
cp ../aether.json.rusanov ./aether.json
mpirun -np 4 ./aether

# post process and plot:
cd UA/output
~/bin/postAether.py -rm
# not sure where plotting code is located....
~/bin/run_plot_block_model_results.py 3DNEU_20110320_*.nc -var="N2 Vertical Wind" -alt=30
~/bin/run_plot_block_model_results.py 3DNEU_20110320_*.nc -var="N2" -alt=30

# back into main test directory:
cd ../../..

# ----------------------------------------------------------------------
# run the hydrostatic test:

cp -R ../../share/run ./run.test_hydro
cd run.test_hydro
cp ../aether.json.hydro ./aether.json
mpirun -np 4 ./aether

# post process and plot:
cd UA/output
~/bin/postAether.py -rm
# not sure where plotting code is located....
~/bin/run_plot_block_model_results.py 3DNEU_20110320_*.nc -var="Vertical Wind" -alt=30
~/bin/run_plot_block_model_results.py 3DNEU_20110320_*.nc -var="N2" -alt=30

# back into main test directory:
cd ../../..
