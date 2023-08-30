#!/bin/sh

rm -rf run.no_advect run.some_advect run.n2_advect

# Run with n2 as the only advected neutral:
cp -R ../../share/run ./run.n2_advect
cd run.n2_advect
cp ../n2_advect.json ./aether.json
./aether
cd UA/output
~/bin/postAether.py -rm

~/bin/aether_plot_simple.py -var=O -alt=300 3DNEU_20110320_001000.nc
~/bin/aether_plot_simple.py -var="O Vertical Wind" -alt=105 3DNEU_20110320_001000.nc
~/bin/aether_plot_simple.py -var=N2 -alt=300 3DNEU_20110320_001000.nc
~/bin/aether_plot_simple.py -var="N2 Vertical Wind" -alt=105 3DNEU_20110320_001000.nc
cd ../../..

# Run with more advected species (specified in earth.in)
cp -R ../../share/run ./run.some_advect
cd run.some_advect
cp ../some_advect.json ./aether.json
./aether
cd UA/output
~/bin/postAether.py -rm

~/bin/aether_plot_simple.py -var=O -alt=300 3DNEU_20110320_001000.nc
~/bin/aether_plot_simple.py -var="O Vertical Wind" -alt=105 3DNEU_20110320_001000.nc
~/bin/aether_plot_simple.py -var=N2 -alt=300 3DNEU_20110320_001000.nc
~/bin/aether_plot_simple.py -var="N2 Vertical Wind" -alt=105 3DNEU_20110320_001000.nc
cd ../../..

