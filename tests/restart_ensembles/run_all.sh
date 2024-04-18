#!/bin/sh

rm -rf run.test
cp -R ../../share/run ./run.test
cd run.test

rm -f aether.json UA/output/*
rm -rf UA/restart* ; mkdir UA/restartOut
rm -rf UA/output_whole
rm -rf UA/output_nonperturbed
rm -rf UA/output_perturbed
rm -rf UA/output_save

cp ../aether.json.whole aether.json
mpirun -np 20 ./aether
cd UA/output ; ../../../../../srcPython/plot_logfiles.py -vars=10 log_m00*.txt ; cd -
cd UA ; mv output output_whole ; mkdir output ; cd -
rm -f aether.json

rm -f aether.json UA/output/* UA/restartOut/*
cp ../aether.json.start aether.json
mpirun -np 20 ./aether
cd UA ; mv restartOut restartIn ; mkdir restartOut ; cd -
cd UA ; cp -R output output_save ; cd ..

rm -f aether.json
cp ../aether.json.restart1 aether.json
mpirun -np 20 ./aether
cd UA/output ; ../../../../../srcPython/plot_logfiles.py -vars=10 log_m00*.txt ; cd -
cd UA ; mv output output_nonperturbed ; cp -R output_save output ; cd -

rm -f aether.json
cp ../aether.json.restart2 aether.json
mpirun -np 20 ./aether
cd UA/output ; ../../../../../srcPython/plot_logfiles.py -vars=10 log_m00*.txt ; cd -
cd UA ; mv output output_perturbed ; mkdir output ; cd -


