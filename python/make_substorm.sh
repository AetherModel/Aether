#!/bin/sh

# base version, no substorm
./electrodynamics_write.py -startdate=19980106 -enddate=19980107 -dt=1 -real -outfile=b980106_vNoSSn.bin
./electrodynamics_write.py -startdate=19980106 -enddate=19980107 -dt=1 -real -outfile=b980106_vNoSSs.bin -south
./amie_read_binary.py -start 140 -end 250 -step 5 b980106_vNoSSn.bin

exit

# base version, no motion
./electrodynamics_write.py -startdate=19980106 -enddate=19980107 -dt=1 -real -outfile=b980106_v0n.bin -substorm
./electrodynamics_write.py -startdate=19980106 -enddate=19980107 -dt=1 -real -outfile=b980106_v0s.bin -south -substorm
./amie_read_binary.py -start 140 -end 250 -step 5 b980106_v0n.bin

# move the onset region
./electrodynamics_write.py -startdate=19980106 -enddate=19980107 -dt=1 -real -outfile=b980106_v1n.bin -substorm -move
./electrodynamics_write.py -startdate=19980106 -enddate=19980107 -dt=1 -real -outfile=b980106_v1s.bin -south -substorm -move
./amie_read_binary.py -start 140 -end 250 -step 5 b980106_v1n.bin
