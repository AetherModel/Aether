#!/bin/sh

# remove old directories:
rm -rf build run.*

# make code with gradient unit test on:
mkdir build
cd build
cmake -DTEST_GRADIENT=Y ../../..
make -j4
cd ..
cp -R ../../share/run ./run.gradient_test
cd run.gradient_test

cp ../aether.json.gradient_test ./aether.json
mpirun -np 6 ./aether
