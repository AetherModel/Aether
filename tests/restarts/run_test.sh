#!/bin/sh

rm -rf run.halves run.whole

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
../../../../../srcPython/postAether.py -alt=-1 -rm

# [O]:
aether_plot_simple.py -var=density_O -alt=250 3DALG_20110320_00*.nc
# Tn:
aether_plot_simple.py -var=Temperature_neutral -alt=250 3DALG_20110320_00*.nc
# [e-]
aether_plot_simple.py -var=density_e- -alt=120 3DALG_20110320_00*.nc
# East Velocity
aether_plot_simple.py -var=velocity_east_neutral -alt=250 3DALG_20110320_00*.nc

cd ../../..

cp -R ../../share/run ./run.whole
cd run.whole
cp ../aether.whole.json ./aether.json
./aether
# plot the output
cd UA/output
../../../../../srcPython/postAether.py -alt=-1 -rm
# [O]:
aether_plot_simple.py -var=density_O -alt=250 3DALG_20110320_00*.nc
# Tn:
aether_plot_simple.py -var=Temperature_neutral -alt=250 3DALG_20110320_00*.nc
# [e-]
aether_plot_simple.py -var=density_e- -alt=120 3DALG_20110320_00*.nc
# East Velocity
aether_plot_simple.py -var=velocity_east_neutral -alt=250 3DALG_20110320_00*.nc

cd ../../..

# you can then look at:
# restarted run results: run.halves/UA/output/*.png
# whole run results: run.whole/UA/output/*.png
# using whatever viewer you have
