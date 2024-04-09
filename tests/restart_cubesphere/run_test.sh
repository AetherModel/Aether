#!/bin/sh

rm -rf run.halves run.whole

cp -R ../../share/run ./run.halves
cd run.halves
# first part of the run
cp ../aether.first.json ./aether.json
mpirun -np 6 ./aether
cd UA ; rm -f restartIn ; cp -R restartOut restartIn ; cd ..
# second part of the run
cp ../aether.second.json ./aether.json
mpirun -np 6 ./aether
# plot the output
cd UA/output
../../../../../srcPython/postAether.py -alt=-1 -rm

# [O]:
aether_plot_simple.py -var=density_O -alt=250 3DALL_20110320_003000.nc -polar
# Eastward Ion Velocity
aether_plot_simple.py -var=velocity_east_ion -alt=250 3DALL_20110320_003000.nc -polar
# [e-]
aether_plot_simple.py -var=density_e- -alt=250 3DALL_20110320_003000.nc -polar

cd ../../..

cp -R ../../share/run ./run.whole
cd run.whole
cp ../aether.whole.json ./aether.json
mpirun -np 6 ./aether
# plot the output
cd UA/output
../../../../../srcPython/postAether.py -alt=-1 -rm
# [O]:
aether_plot_simple.py -var=density_O -alt=250 3DALL_20110320_003000.nc -polar
# Eastward Ion Velocity 
aether_plot_simple.py -var=velocity_east_ion -alt=250 3DALL_20110320_003000.nc -polar
# [e-]
aether_plot_simple.py -var=density_e- -alt=250 3DALL_20110320_003000.nc -polar

cd ../../..

# you can then look at:
# restarted run results: run.halves/UA/output/*.png
# whole run results: run.whole/UA/output/*.png
# using whatever viewer you have
