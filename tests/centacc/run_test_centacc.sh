# Tests for centripetal acceleration 
# !/bin/sh 

rm -rf ./run.test.centacc
cp -R ../../share/run ./run.test.centacc

cd ./run.test.centacc

# create & copy aether files 
cp ../aether.json.true ./aether.json

./aether

# post-processing
cd UA/output 
python3 ../../../../../srcPython/postAether.py -rm

# plot latitudinal centripetal acceleration
python3 ~/Aether/aetherpy/bin/run_plot_block_model_results.py 3DMMT_20110320_002000.nc -alt=30 -var=Latitudinal

# plot longitudinal centripetal acceleration 
python3 ~/Aether/aetherpy/bin/run_plot_block_model_results.py 3DMMT_20110320_002000.nc -alt=30 -var=Longitudinal

# plot radial centripetal acceleration 
python3 ~/Aether/aetherpy/bin/run_plot_block_model_results.py 3DMMT_20110320_002000.nc -alt=30 -var=Radial

cd ..
mv output output.centacc
mkdir output

cd ..

cp ../aether.json.false ./aether.json

./aether

# post-processing
cd UA/output
python3 ../../../../../srcPython/postAether.py -rm

# plot latitudinal centripetal acceleration
python3 ~/Aether/aetherpy/bin/run_plot_block_model_results.py 3DMMT_20110320_002000.nc -alt=30 -var=Latitudinal

# plot longitudinal centripetal acceleration 
python3 ~/Aether/aetherpy/bin/run_plot_block_model_results.py 3DMMT_20110320_002000.nc -alt=30 -var=Longitudinal

# plot radial centripetal acceleration 
python3 ~/Aether/aetherpy/bin/run_plot_block_model_results.py 3DMMT_20110320_002000.nc -alt=30 -var=Radial

cd ..
mv output output.no_centacc
