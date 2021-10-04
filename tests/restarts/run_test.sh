#!/bin/sh

cp -R ../../share/run ./run.whole
cd run.whole
cp ../aether.whole.json ./aether.json
./aether
# plot the output
cd UA/output
# [O]:
run_plot_model_results.py -var=3 -alt=350 3DALL_20110320_010000.bin
# Tn:
run_plot_model_results.py -var=14 -alt=350 3DALL_20110320_010000.bin
# [e-]
run_plot_model_results.py -var=23 -alt=350 3DALL_20110320_010000.bin

cp -R ../../share/run ./run.halves
cd run.halves
# first part of the run
cp ../aether.first.json ./aether.json
./aether
# second part of the run
cp ../aether.second.json ./aether.json
./aether
# plot the output
cd UA/output
# [O]:
run_plot_model_results.py -var=3 -alt=350 3DALL_20110320_010000.bin
# Tn:
run_plot_model_results.py -var=14 -alt=350 3DALL_20110320_010000.bin
# [e-]
run_plot_model_results.py -var=23 -alt=350 3DALL_20110320_010000.bin

