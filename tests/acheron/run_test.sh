#!/bin/sh


rm -rf run.acheron
 
cp -r ../../share/run ./run.acheron
cd run.acheron
cp ../aether.json.acheron ./aether.json

./aether

cd UA/output

~/bin/postAether.py -rm
~/bin/aether_plot.py -var=Temperature -alt=30 3DALL_20110320_001000.nc
