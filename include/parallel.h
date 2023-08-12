// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_PARALLEL_H_
#define INCLUDE_PARALLEL_H_

/// Need MPI (message passing interface) to do parallel stuff: 
#include "mpi.h"

/// number of processors in whole simulation
extern int nProcs;
/// processor number in whole simulation
extern int iProc;
/// processor number as a string
extern std::string cProc;

/// number of ensemble members
extern int nMembers;
/// ensemble member number
extern int iMember;
/// ensemble member number as a character
extern std::string cMember;

/// number of grid blocks in each ensemble
extern int nGrids;
/// number of the grid block
extern int iGrid;
/// number of the grid block as a string
extern std::string cGrid;

/// communicator for all of aether
extern MPI_Comm aether_comm;

/**********************************************************************
  \brief initialize mpi and figure out ensembles and grid blocks
**/
bool init_parallel(Quadtree &quadtree);

/**********************************************************************
  \brief Pack variables for message passing
  \param value is variable to pack
  \param packed is the buffer to pack it in
  \param iCounter is the location in the buffer to start
  \param nG is the number of ghost cells
  \param iDir is the direction of the message pass:
     0 - left -> right
     1 - bottom -> top
     2 - right -> left
     3 - top -> bottom
**/

bool pack_border(const arma_cube &value,
		 precision_t *packed,
		 int64_t *iCounter,
		 int64_t nG,
		 int iDir);

/**********************************************************************
  \brief Unpack variable buffer after message pass
  \param value is variable to unpack into
  \param packed is the buffer that contains the values
  \param iCounter is the location in the buffer to start
  \param nG is the number of ghost cells
  \param iDir is the direction of the message pass:
     0 - left -> right
     1 - bottom -> top
     2 - right -> left
     3 - top -> bottom
  \param IsPole indicates whether the message passing happened across the pole,
     in which case the longitudes (Xs) have to be mirrored
**/

bool unpack_border(arma_cube &value,
		   precision_t *packed,
		   int64_t *iCounter,
		   int64_t nG,
		   int iDir,
		   bool DoReverseX,
		   bool DoReverseY,
		   bool XbecomesY);

#endif  // INCLUDE_PARALLEL_H_
