#!/bin/sh

# In theory, you should be able to run this script and everything
# should just work.  There are two places where you need to change
# things...

rm -rf Aether Aether.test_from_scrath
git clone https://github.com/AetherModel/Aether
mv Aether Aether.test_from_scrath
cd Aether.test_from_scrath

# 1. Change the line below to set the BRANCH to test:
git checkout indices

# 2. Change the line before to set the OS to test:
cd src ; rm -f Makefile.OS ; cp Makefile.ubuntu Makefile.OS ; cd ..

make
make rundir ; mv run run.test_from_scratch
cd run.test_from_scratch
./aether.exe

../python/plot_model_results.py -var=24 -alt=110 3DALL_20110320_010000.nc
../python/plot_model_results.py -var=14 -alt=300 3DALL_20110320_010000.nc
