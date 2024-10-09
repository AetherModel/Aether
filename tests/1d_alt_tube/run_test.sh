#!/bin/sh


rm -rf run.1d
 
cp -r ../../share/run ./run.1d
cd run.1d
cp ../aether.json.1d_alt_tube ./aether.json

./aether

cd UA/output

~/bin/postAether.py -rm -alt=-1
~/bin/aether_plot_simple.py 3DALG_20110320_*.nc
