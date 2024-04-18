#!/bin/sh

# In theory, you should be able to run this script and everything
# should just work.  There are two places where you need to change
# things...

rm -rf Aether Aether.test_from_scratch
git clone https://github.com/AetherModel/Aether
mv Aether Aether.test_from_scratch
cd Aether.test_from_scratch

# 1. Change the line below to set the BRANCH to test:
git checkout develop

# 2. Change the line before to set the OS to test:
mkdir build
cd build
cmake ..

make -j16
cd ..
cp -R share/run ./run.test_from_scratch
cd run.test_from_scratch
./aether

../python/plot_model_results.py -var=24 -alt=110 3DALL_20110320_010000.nc
../python/plot_model_results.py -var=14 -alt=300 3DALL_20110320_010000.nc
