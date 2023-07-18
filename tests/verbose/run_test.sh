#!/bin/sh

rm -rf run.test1 run.test2 run.test3

cp -R ../../share/run ./run.test1
cd run.test1
cp ../aether1.json ./aether.json
./aether

cp -R ../../share/run ./run.test2
cd run.test2
cp ../aether2.json ./aether.json
./aether

cp -R ../../share/run ./run.test3
cd run.test3
cp ../aether3.json ./aether.json
./aether