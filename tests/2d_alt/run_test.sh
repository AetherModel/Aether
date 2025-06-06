#!/bin/sh


rm -rf run.2d
 
cp -r ../../share/run ./run.2d
cd run.2d
cp ../aether.json.2d ./aether.json

./aether

cd UA/output

~/bin/postAether.py -rm
~/bin/aether_plot_simple.py -var=Temperature_neutral -alt=2 3DALG_20110320_130000.nc
