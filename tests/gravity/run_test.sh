#!/bin/sh


rm -rf run.grav
 
cp -r ../../share/run ./run.grav
cd run.grav
cp ../aether_first.json aether.json

./aether

cd UA/output

../../../../../srcPython/postAether.py -rm
pwd

~/bin/run_plot_block_model_results.py -var=Gvertical 3DGRA_20110320_000000.nc
~/bin/run_plot_block_model_results.py -var=Gnorth 3DGRA_20110320_000000.nc
~/bin/run_plot_block_model_results.py -var=Geast 3DGRA_20110320_000000.nc
pwd
mv 3DGRA_20110320_000000_Gvertical.png  ../../../output_pngs
mv 3DGRA_20110320_000000_Gnorth.png ../../../output_pngs
mv 3DGRA_20110320_000000_Geast.png ../../../output_pngs

rm -rf *png *.bin *.json *.nc

cd ../..

cp ../aether_second.json aether.json

./aether

cd UA/output

../../../../../srcPython/postAether.py -rm
pwd

~/bin/run_plot_block_model_results.py -var=Gvertical 3DGRA_20110320_000000.nc
~/bin/run_plot_block_model_results.py -var=Gnorth 3DGRA_20110320_000000.nc
~/bin/run_plot_block_model_results.py -var=Geast 3DGRA_20110320_000000.nc

mv 3DGRA_20110320_000000_Gvertical.png ../../../output_pngs/Oblate_Gvertical.png
mv 3DGRA_20110320_000000_Gnorth.png ../../../output_pngs/Oblate_Gnorth.png
mv 3DGRA_20110320_000000_Geast.png ../../../output_pngs/Oblate_Geast.png

rm -rf *png *.bin *.json *.nc

cd ../..

cp ../aether_third.json aether.json

./aether

cd UA/output

../../../../../srcPython/postAether.py -rm

~/bin/run_plot_block_model_results.py -var=Gvertical 3DGRA_20110320_000000.nc
~/bin/run_plot_block_model_results.py -var=Gnorth 3DGRA_20110320_000000.nc
~/bin/run_plot_block_model_results.py -var=Geast 3DGRA_20110320_000000.nc

mv 3DGRA_20110320_000000_Gvertical.png ../../../output_pngs/Oblate_withJ2_Gvertical.png
mv 3DGRA_20110320_000000_Gnorth.png ../../../output_pngs/Oblate_withJ2_Gnorth.png
mv 3DGRA_20110320_000000_Geast.png ../../../output_pngs/Oblate_withJ2_Geast.png

rm -rf *png *.bin *.json *.nc
