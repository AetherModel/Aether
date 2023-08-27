#!/bin/sh

rm -rf run.test1

cp -R ../../share/run ./run.test1
cd run.test1
cp ../aether1.json ./aether.json
./aether