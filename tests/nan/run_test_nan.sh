# Tests for checking for nan
# !/bin/sh 

rm -rf ./run.test.insert_nan
cp -R ../../share/run ./run.test.insert_nan

cd ./run.test.insert_nan

# create & copy aether files 
cp ../aether.json.insert_nan ./aether.json

./aether

cd ..

rm -rf ./run.test.check_nan
cp -R ../../share/run ./run.test.check_nan

cd ./run.test.check_nan

# create & copy aether files 
cp ../aether.json.check_nan ./aether.json

./aether
