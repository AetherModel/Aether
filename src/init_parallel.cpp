// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "../include/aether.h"

// -----------------------------------------------------------------------------
// Initialize the parallel aspect of Aether:
//   1. Ensembles
//   2. Domain decomposition
// -----------------------------------------------------------------------------

int nProcs;
int iProc;

int nMembers;
int iMember;
int nGrids;
int iGrid;

std::string cProc;
std::string cMember;
std::string cGrid;

int init_parallel(Inputs input,
		  Report &report) {

  int iErr = 0;

  MPI_Init(NULL, NULL);

  // Get the number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

  // Get the rank of the process
  MPI_Comm_rank(MPI_COMM_WORLD, &iProc);  

  // Modify the verbosity of the code by turning of verbose on all
  // processors except specified processor:
  if (iProc != input.get_verbose_proc())
    report.set_verbose(-1);
  
  nMembers = input.get_nMembers();

  // Check to see if we have enough processors to do this stuff:

  int nBlocksLonGeo = input.get_nBlocksLonGeo();
  int nBlocksLatGeo = input.get_nBlocksLatGeo();
  nGrids = nBlocksLonGeo * nBlocksLatGeo;
  int nProcsNeeded = nMembers * nGrids;

  if (nProcsNeeded <= nProcs) {
    // Get Ensemble member number and grid number:
    iMember = iProc / nGrids;
    iGrid = iProc % nGrids;
    if (report.test_verbose(2))
      std::cout << "iProc : " << iProc
		<< "; iMember : " << iMember
		<< "; iGrid : " << iGrid << "\n";
  } else {
    std::cout << "Number of processors needed is not set correctly!\n";
    std::cout << "nProcs needs to be greater that multiplication of :\n";
    std::cout << "   nMembers : " << nMembers << "\n";
    std::cout << "   nBlocksLonGeo : " << nBlocksLonGeo << "\n";
    std::cout << "   nBlocksLatGeo : " << nBlocksLatGeo << "\n";
    std::cout << "   total needed : " << nProcsNeeded << "\n";
    std::cout << "   which is greater than nProcs : " << nProcs << "\n";
    iErr = 1;
  }
  // Create strings to allow for easier output filename creation:
  cProc = "p" + tostr(iProc, 4);
  cMember = "m" + tostr(iMember, 4);
  cGrid = "g" + tostr(iGrid, 4);
  return iErr;  
}
