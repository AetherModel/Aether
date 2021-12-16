// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_PARALLEL_H_
#define INCLUDE_PARALLEL_H_

/// Need MPI (message passing interface) to do parallel stuff: 
#include "mpi.h"

extern int nProcs;
extern int iProc;
extern std::string cProc;

extern int nMembers;
extern int iMember;
extern std::string cMember;

extern int nGrids;
extern int iGrid;
extern std::string cGrid;

extern MPI_Comm aether_communicator;

int init_parallel(Inputs &input, Report &report);

#endif  // INCLUDE_PARALLEL_H_
