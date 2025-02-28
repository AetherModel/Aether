#!/bin/sh

NPROC=5
TOTALTIME=180
MPI=/usr/bin/mpirun
NTIMES=2

# include -dowhole to run the whole simulation as comparison:

../../srcPython/run_restarts.py -totaltime=${TOTALTIME} -mpi=${MPI} -rundir=../../share/run -ensembles=${NPROC} -restarts=${NTIMES}

