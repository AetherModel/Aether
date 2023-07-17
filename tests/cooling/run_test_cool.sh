# Tests for NO & O cooling. Plots for O_cooling, NO_cooling, both, or none. 
# !/bin/sh 

rm -rf ./run.test.cooling
cp -R ../../share/run ./run.test.cooling

cd ./run.test.cooling

# create & copy aether files 
cp ../aether.json.no_true ./aether.json

./aether

# post-processing
cd UA/output
~/bin/postAether.py -rm

# plot O cooling
~/bin/run_plot_block_model_results.py 3DTHR_20110320_001000.nc -var="O Rad Cooling" -alt=10

# plot NO cooling
~/bin/run_plot_block_model_results.py 3DTHR_20110320_001000.nc -var="NO Rad Cooling" -alt=10


cd ..
mv output output.no.cooling
mkdir output

cd ..

cp ../aether.json.o_true ./aether.json

./aether

# post-processing
cd UA/output
~/bin/postAether.py -rm

# plot O cooling
~/bin/run_plot_block_model_results.py 3DTHR_20110320_001000.nc -var="O Rad Cooling" -alt=10

# plot NO cooling
~/bin/run_plot_block_model_results.py 3DTHR_20110320_001000.nc -var="NO Rad Cooling" -alt=10


cd ..
mv output output.o.cooling
mkdir output

cd ..

cp ../aether.json.all_true ./aether.json

./aether

# post-processing
cd UA/output
~/bin/postAether.py -rm

# plot O cooling
~/bin/run_plot_block_model_results.py 3DTHR_20110320_001000.nc -var="O Rad Cooling" -alt=10

# plot NO cooling
~/bin/run_plot_block_model_results.py 3DTHR_20110320_001000.nc -var="NO Rad Cooling" -alt=10

cd ..
mv output output.all.cooling
mkdir output

cd ..

cp ../aether.json.all_false ./aether.json

./aether

# post-processing
cd UA/output
~/bin/postAether.py -rm

# plot O cooling
~/bin/run_plot_block_model_results.py 3DTHR_20110320_001000.nc -var="O Rad Cooling" -alt=10

# plot NO cooling
~/bin/run_plot_block_model_results.py 3DTHR_20110320_001000.nc -var="NO Rad Cooling" -alt=10

cd ..
mv output output.none.cooling
mkdir output

# To open the png files, go into individual output file (ex aether.json.no_true) 
# open 3DTHR_20110320_001000_ORadCooling.png
# open 3DTHR_20110320_001000_NORadCooling.png
