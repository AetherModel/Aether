#!/bin/sh

PLOTTER=~/Software/PyITM/bin/plot_alt_wpolar.py

# remove old directories:
rm -rf run.*

# run 59 wavelengths, FISM
RUN=59_fism
cp -R ../../share/run ./run.${RUN}
cd run.${RUN}
cp ../aether.json.${RUN} ./aether.json
mpirun -np 4 ./aether
cd UA/output
../../../../../srcPython/postAether.py -rm
# This assumes pyitm is installed and the plotter is in the bin directory:
${PLOTTER} -var=Tn -alt=300 3DALG_20110320_001000.nc
${PLOTTER} -var=O2+ -alt=120 3DALM_20110320_001000.nc
# into test directory
cd ../../..

# run 59 wavelengths, FISM
RUN=59_fism_nopei
cp -R ../../share/run ./run.${RUN}
cd run.${RUN}
cp ../aether.json.${RUN} ./aether.json
mpirun -np 4 ./aether
cd UA/output
../../../../../srcPython/postAether.py -rm
# This assumes pyitm is installed and the plotter is in the bin directory:
${PLOTTER} -var=Tn -alt=300 3DALG_20110320_001000.nc
${PLOTTER} -var=O2+ -alt=120 3DALM_20110320_001000.nc
# into test directory
cd ../../..

# run 59 wavelengths, neuvac
RUN=59_neuvac
cp -R ../../share/run ./run.${RUN}
cd run.${RUN}
cp ../aether.json.${RUN} ./aether.json
mpirun -np 4 ./aether
cd UA/output
../../../../../srcPython/postAether.py -rm
# This assumes pyitm is installed and the plotter is in the bin directory:
${PLOTTER} -var=Tn -alt=300 3DALG_20110320_001000.nc
${PLOTTER} -var=O2+ -alt=120 3DALM_20110320_001000.nc
# into test directory
cd ../../.. 

# run 59 wavelengt, neuvac
RUN=59_neuvac_nopei
cp -R ../../share/run ./run.${RUN}
cd run.${RUN}
cp ../aether.json.${RUN} ./aether.json
mpirun -np 4 ./aether
cd UA/output
../../../../../srcPython/postAether.py -rm
# This assumes pyitm is installed and the plotter is in the bin directory:
${PLOTTER} -var=Tn -alt=300 3DALG_20110320_001000.nc
${PLOTTER} -var=O2+ -alt=120 3DALM_20110320_001000.nc
# into test directory
cd ../../.. 

# run 59 wavelengths, neuvac
RUN=59_euvac
cp -R ../../share/run ./run.${RUN}
cd run.${RUN}
cp ../aether.json.${RUN} ./aether.json
mpirun -np 4 ./aether
cd UA/output
../../../../../srcPython/postAether.py -rm
# This assumes pyitm is installed and the plotter is in the bin directory:
${PLOTTER} -var=Tn -alt=300 3DALG_20110320_001000.nc
${PLOTTER} -var=O2+ -alt=120 3DALM_20110320_001000.nc
# into test directory
cd ../../.. 

# run 22 wavelengths, euvac
RUN=22_euvac
cp -R ../../share/run ./run.${RUN}
cd run.${RUN}
cp ../aether.json.${RUN} ./aether.json
mpirun -np 4 ./aether
cd UA/output
../../../../../srcPython/postAether.py -rm
# This assumes pyitm is installed and the plotter is in the bin directory:
${PLOTTER} -var=Tn -alt=300 3DALG_20110320_001000.nc
${PLOTTER} -var=O2+ -alt=120 3DALM_20110320_001000.nc
# into test directory
cd ../../.. 

# run 22 wavelengths, HFG
RUN=22_hfg
cp -R ../../share/run ./run.${RUN}
cd run.${RUN}
cp ../aether.json.${RUN} ./aether.json
mpirun -np 4 ./aether
cd UA/output
../../../../../srcPython/postAether.py -rm
# This assumes pyitm is installed and the plotter is in the bin directory:
${PLOTTER} -var=Tn -alt=300 3DALG_20110320_001000.nc
${PLOTTER} -var=O2+ -alt=120 3DALM_20110320_001000.nc
# into test directory
cd ../../.. 

# run 59 wavelengths, euvac
RUN=59_euvac
cp -R ../../share/run ./run.${RUN}
cd run.${RUN}
cp ../aether.json.${RUN} ./aether.json
mpirun -np 4 ./aether
cd UA/output
../../../../../srcPython/postAether.py -rm
# This assumes pyitm is installed and the plotter is in the bin directory:
${PLOTTER} -var=Tn -alt=300 3DALG_20110320_001000.nc
${PLOTTER} -var=O2+ -alt=120 3DALM_20110320_001000.nc
# into test directory
cd ../../.. 

# run 37 wavelengths, euvac
RUN=37_euvac
cp -R ../../share/run ./run.${RUN}
cd run.${RUN}
cp ../aether.json.${RUN} ./aether.json
mpirun -np 4 ./aether
cd UA/output
../../../../../srcPython/postAether.py -rm
# This assumes pyitm is installed and the plotter is in the bin directory:
${PLOTTER} -var=Tn -alt=300 3DALG_20110320_001000.nc
${PLOTTER} -var=O2+ -alt=120 3DALM_20110320_001000.nc
# into test directory
cd ../../.. 


~/Software/PyITM/bin/plot_logfile.py run.*/UA/output/log_geo.txt -vars 17 10

