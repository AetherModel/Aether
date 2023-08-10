# Tests for checking for nan
# !/bin/sh 

rm -rf ./run.test.nan
cp -R ../../share/run ./run.test.nan

cd ./run.test.nan

# create & copy aether files 
cp ../aether.json.check_nan ./aether.json

./aether