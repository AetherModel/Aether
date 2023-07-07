#!/bin/sh


rm -rf run.grav
 
cp -r ../../share/run ./run.grav
cd run.grav
cp ../aether_first.json aether.json

./aether

cd UA/output

../../../../../srcPython/postAether.py -rm
pwd

../../../../../srcPython/aetherpy/bin/run_plot_block_model_results.py -var=Galtitude 3DGRA_20110320_000000.nc
../../../../../srcPython/aetherpy/bin/run_plot_block_model_results.py -var=Glatitude 3DGRA_20110320_000000.nc
../../../../../srcPython/aetherpy/bin/run_plot_block_model_results.py -var=Glongitude 3DGRA_20110320_000000.nc
pwd
mv 3DGRA_20110320_000000_Galtitude.png  ../../../../output_pngs
mv 3DGRA_20110320_000000_Glatitude.png ../../../../output_pngs
mv 3DGRA_20110320_000000_Glongitude.png ../../../../output_pngs

cd ../..

cp ../aether_second.json aether.json

./aether

cd UA/output

../../../../../srcPython/postAether.py -rm
pwd

../../../../../srcPython/aetherpy/bin/run_plot_block_model_results.py -var=Galtitude 3DGRA_20110320_000000.nc
../../../../../srcPython/aetherpy/bin/run_plot_block_model_results.py -var=Glatitude 3DGRA_20110320_000000.nc
../../../../../srcPython/aetherpy/bin/run_plot_block_model_results.py -var=Glongitude 3DGRA_20110320_000000.nc

mv 3DGRA_20110320_000000_Galtitude.png ../../../../output_pngs/Oblate_Galtitude.png
mv 3DGRA_20110320_000000_Glatitude.png ../../../../output_pngs/Oblate_Glatitude.png
mv 3DGRA_20110320_000000_Glongitude.png ../../../../output_pngs/Oblate_Glongitude.png

cd ../..

cp ../aether_third.json aether.json

./aether

cd UA/output

../../../../../srcPython/postAether.py -rm

../../../../../srcPython/aetherpy/bin/run_plot_block_model_results.py -var=Galtitude 3DGRA_20110320_000000.nc
../../../../../srcPython/aetherpy/bin/run_plot_block_model_results.py -var=Glatitude 3DGRA_20110320_000000.nc
../../../../../srcPython/aetherpy/bin/run_plot_block_model_results.py -var=Glongitude 3DGRA_20110320_000000.nc

mv 3DGRA_20110320_000000_Galtitude.png ../../../../output_pngs/Oblate_withJ2_Galtitude.png
mv 3DGRA_20110320_000000_Glatitude.png ../../../../output_pngs/Oblate_withJ2_Glatitude.png
mv 3DGRA_20110320_000000_Glongitude.png ../../../../output_pngs/Oblate_withJ2_Glongitude.png
