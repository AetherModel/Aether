#!/bin/sh

PLOTTER=~/Software/PyITM/bin/plot_alt_wpolar.py

# remove old directories:
rm -rf run.*

RUN=cube_cube
PE=6
rm -rf ./run.${RUN}
cp -R ../../share/run ./run.${RUN}
cd run.${RUN}
cp ../aether_${RUN}.json ./aether.json
mpirun -np ${PE} ./aether
../../../srcPython/postAether.py -rm
# This assumes pyitm is installed and the plotter is in the bin directory:
cd UA/output
${PLOTTER} -var=Tn -alt=300 3DALG_20110320_000100.nc
${PLOTTER} -var=O2+ -alt=120 3DALM_20110320_000100.nc
# into test directory
cd ../../..

RUN=sphere_sphere
PE=1
rm -rf ./run.${RUN}
cp -R ../../share/run ./run.${RUN}
cd run.${RUN}
cp ../aether_${RUN}.json ./aether.json
mpirun -np ${PE} ./aether
../../../srcPython/postAether.py -rm
# This assumes pyitm is installed and the plotter is in the bin directory:
cd UA/output
${PLOTTER} -var=Tn -alt=300 3DALG_20110320_000100.nc
${PLOTTER} -var=O2+ -alt=120 3DALM_20110320_000100.nc
# into test directory
cd ../../..

RUN=sphere_sphere
PE=4
rm -rf ./run.${RUN}
cp -R ../../share/run ./run.${RUN}
cd run.${RUN}
cp ../aether_${RUN}.json ./aether.json
mpirun -np ${PE} ./aether
../../../srcPython/postAether.py -rm
# This assumes pyitm is installed and the plotter is in the bin directory:
cd UA/output
${PLOTTER} -var=Tn -alt=300 3DALG_20110320_000100.nc
${PLOTTER} -var=O2+ -alt=120 3DALM_20110320_000100.nc
# into test directory
cd ../../..

RUN=sphere4_sphere4
PE=4
rm -rf ./run.${RUN}
cp -R ../../share/run ./run.${RUN}
cd run.${RUN}
cp ../aether_${RUN}.json ./aether.json
mpirun -np ${PE} ./aether
../../../srcPython/postAether.py -rm
# This assumes pyitm is installed and the plotter is in the bin directory:
cd UA/output
${PLOTTER} -var=Tn -alt=300 3DALG_20110320_000100.nc
${PLOTTER} -var=O2+ -alt=120 3DALM_20110320_000100.nc
# into test directory
cd ../../..

RUN=sphere6_sphere6
PE=6
rm -rf ./run.${RUN}
cp -R ../../share/run ./run.${RUN}
cd run.${RUN}
cp ../aether_${RUN}.json ./aether.json
mpirun -np ${PE} ./aether
../../../srcPython/postAether.py -rm
# This assumes pyitm is installed and the plotter is in the bin directory:
cd UA/output
${PLOTTER} -var=Tn -alt=300 3DALG_20110320_000100.nc
${PLOTTER} -var=O2+ -alt=120 3DALM_20110320_000100.nc
# into test directory
cd ../../..

RUN=cube_sphere6
PE=6
rm -rf ./run.${RUN}
cp -R ../../share/run ./run.${RUN}
cd run.${RUN}
cp ../aether_${RUN}.json ./aether.json
mpirun -np ${PE} ./aether
../../../srcPython/postAether.py -rm
# This assumes pyitm is installed and the plotter is in the bin directory:
cd UA/output
${PLOTTER} -var=Tn -alt=300 3DALG_20110320_000100.nc
${PLOTTER} -var=O2+ -alt=120 3DALM_20110320_000100.nc
# into test directory
cd ../../..


RUN=sphere4_dipole4
PE=4
rm -rf ./run.${RUN}
cp -R ../../share/run ./run.${RUN}
cd run.${RUN}
cp ../aether_${RUN}.json ./aether.json
mpirun -np ${PE} ./aether
../../../srcPython/postAether.py -rm
# This assumes pyitm is installed and the plotter is in the bin directory:
cd UA/output
${PLOTTER} -var=Tn -alt=300 3DALG_20110320_000100.nc
${PLOTTER} -var=O2+ -alt=120 3DALM_20110320_000100.nc
# into test directory
cd ../../..

RUN=cube_dipole6
PE=6
rm -rf ./run.${RUN}
cp -R ../../share/run ./run.${RUN}
cd run.${RUN}
cp ../aether_${RUN}.json ./aether.json
mpirun -np ${PE} ./aether
../../../srcPython/postAether.py -rm
# This assumes pyitm is installed and the plotter is in the bin directory:
cd UA/output
${PLOTTER} -var=Tn -alt=300 3DALG_20110320_000100.nc
${PLOTTER} -var=O2+ -alt=120 3DALM_20110320_000100.nc
# into test directory
cd ../../..

RUN=sphere6_dipole6
PE=6
rm -rf ./run.${RUN}
cp -R ../../share/run ./run.${RUN}
cd run.${RUN}
cp ../aether_${RUN}.json ./aether.json
mpirun -np ${PE} ./aether
../../../srcPython/postAether.py -rm
# This assumes pyitm is installed and the plotter is in the bin directory:
cd UA/output
${PLOTTER} -var=Tn -alt=300 3DALG_20110320_000100.nc
${PLOTTER} -var=O2+ -alt=120 3DALM_20110320_000100.nc
# into test directory
cd ../../..

~/Software/PyITM/bin/plot_logfile.py run.*/UA/output/log_geo.txt -vars 17 10
