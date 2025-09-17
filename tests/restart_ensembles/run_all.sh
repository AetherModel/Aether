#!/bin/sh

# cubesphere has 6 blocks:
NBLOCKS=6
# run with 3 members:
NMEMBERS=2
# run for a total of 180s
TOTALTIME=180
# this is the mpi command
MPI=/usr/bin/mpirun
# stop this many times
NTIMES=2

# include -dowhole to run the whole simulation as comparison:

../../srcPython/run_restarts.py -totaltime=${TOTALTIME} -mpi=${MPI} -rundir=../../share/run -ensembles=${NMEMBERS} -restarts=${NTIMES} -blocks=${NBLOCKS}

