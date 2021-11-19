#!/bin/sh

rm -rf run.halves run.whole

cp -R ../../share/run ./run.whole
cd run.whole
cp ../aether.whole.json ./aether.json
./aether
# plot the output
cd UA/output
# [O]:
run_plot_model_results.py -var=3 -alt=250 3DALL_20110320_010000.nc
# Tn:
run_plot_model_results.py -var=14 -alt=250 3DALL_20110320_010000.nc
# [e-]
run_plot_model_results.py -var=23 -alt=250 3DALL_20110320_010000.nc

cd ../../..

cp -R ../../share/run ./run.halves
cd run.halves
# first part of the run
cp ../aether.first.json ./aether.json
./aether
cd UA ; rm -f restartIn ; cp -R restartOut restartIn ; cd ..
# second part of the run
cp ../aether.second.json ./aether.json
./aether
# plot the output
cd UA/output
# [O]:
run_plot_model_results.py -var=3 -alt=250 3DALL_20110320_010000.nc
# Tn:
run_plot_model_results.py -var=14 -alt=250 3DALL_20110320_010000.nc
# [e-]
run_plot_model_results.py -var=23 -alt=250 3DALL_20110320_010000.nc

cd ../../..

# you can then look at:
# restarted run results: run.halves/UA/output/*.png
# whole run results: run.whole/UA/output/*.png
# using whatever viewer you have
