#! /usr/bin/bash

# remove old directory:
rm -rf run.no_errors

# run without any errors to check it works:
cp -R ../../share/run ./run.no_errors
cd run.no_errors
./aether

#now we get into the real testing
cd ../
cd run.errors/UA/inputs
for i in 0 1 2 3 4 5; 
do 
    echo
    echo Test $i: checking for successful crash 
    echo ......................
    sed -i "s/chemistry_earth_richards.csv/chemistry_check_$i.csv/g" defaults.json; 
    cd ../../
    ./aether;
    cd UA/inputs;
    sed -i "s/chemistry_check_$i.csv/chemistry_earth_richards.csv/g" defaults.json; 
done;