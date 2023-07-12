#!/bin/sh

rm -rf run.logfile

cp -R ../../share/run ./run.logfile
cd run.logfile
cp ../aether.logfile.json ./aether.json
./aether
../../../srcPython/logfile_read.py
