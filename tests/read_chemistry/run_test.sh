#!/bin/bash

# remove old directory:
rm -rf run.no_errors run.errors

# run without any errors to check it works:
cp -R ../../share/run ./run.no_errors
cd run.no_errors
./aether
cd ../

#now we get into the real testing:
cp -R ../../share/run ./run.errors
cd run.errors

cp -f ../aether.json.errors ./aether.json

for i in 0 1 2 3 4 5; 
do 
    echo
    echo Test $i: checking for successful crash 
    echo ......................
    rm -f ./chemistry_check.csv
    ln -s ../chemistry_check_$i.csv ./chemistry_check.csv

    ./aether;

done;
