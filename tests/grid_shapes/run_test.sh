#!/bin/sh

RUN=sphere_sphere
PE=1
rm -rf ./run.${RUN}
cp -R ../../share/run ./run.${RUN}
cd run.${RUN}
cp ../aether_${RUN}.json ./aether.json
mpirun -np ${PE} ./aether
../../../srcPython/postAether.py -rm
cd ..

RUN=sphere_sphere
PE=4
rm -rf ./run.${RUN}
cp -R ../../share/run ./run.${RUN}
cd run.${RUN}
cp ../aether_${RUN}.json ./aether.json
mpirun -np ${PE} ./aether
../../../srcPython/postAether.py -rm
cd ..

RUN=sphere4_sphere4
PE=4
rm -rf ./run.${RUN}
cp -R ../../share/run ./run.${RUN}
cd run.${RUN}
cp ../aether_${RUN}.json ./aether.json
mpirun -np ${PE} ./aether
../../../srcPython/postAether.py -rm
cd ..

RUN=sphere6_sphere6
PE=6
rm -rf ./run.${RUN}
cp -R ../../share/run ./run.${RUN}
cd run.${RUN}
cp ../aether_${RUN}.json ./aether.json
mpirun -np ${PE} ./aether
../../../srcPython/postAether.py -rm
cd ..

RUN=cube_cube
PE=6
rm -rf ./run.${RUN}
cp -R ../../share/run ./run.${RUN}
cd run.${RUN}
cp ../aether_${RUN}.json ./aether.json
mpirun -np ${PE} ./aether
../../../srcPython/postAether.py -rm
cd ..

RUN=cube_sphere6
PE=6
rm -rf ./run.${RUN}
cp -R ../../share/run ./run.${RUN}
cd run.${RUN}
cp ../aether_${RUN}.json ./aether.json
mpirun -np ${PE} ./aether
../../../srcPython/postAether.py -rm
cd ..

RUN=sphere4_dipole4
PE=4
rm -rf ./run.${RUN}
cp -R ../../share/run ./run.${RUN}
cd run.${RUN}
cp ../aether_${RUN}.json ./aether.json
mpirun -np ${PE} ./aether
../../../srcPython/postAether.py -rm
cd ..

RUN=cube_dipole6
PE=6
rm -rf ./run.${RUN}
cp -R ../../share/run ./run.${RUN}
cd run.${RUN}
cp ../aether_${RUN}.json ./aether.json
mpirun -np ${PE} ./aether
../../../srcPython/postAether.py -rm
cd ..

RUN=sphere6_dipole6
PE=6
rm -rf ./run.${RUN}
cp -R ../../share/run ./run.${RUN}
cd run.${RUN}
cp ../aether_${RUN}.json ./aether.json
mpirun -np ${PE} ./aether
../../../srcPython/postAether.py -rm
cd ..

