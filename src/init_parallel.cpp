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

MPI_Comm aether_comm;

bool init_parallel(Inputs &input,
                   Quadtree &quadtree,
                   Report &report) {

  bool DidWork = true;

  MPI_Init(NULL, NULL);

  aether_comm = MPI_COMM_WORLD;

  // Get the number of processes
  MPI_Comm_size(aether_comm, &nProcs);

  // Get the rank of the process
  MPI_Comm_rank(aether_comm, &iProc);

  // Modify the verbosity of the code by turning of verbose on all
  // processors except specified processor:
  if (iProc != input.get_verbose_proc())
    report.set_verbose(-1);

  nMembers = input.get_nMembers();
  nGrids = nProcs / nMembers;

  int64_t nProcsPerNode = nGrids / quadtree.nRootNodes;

  if (report.test_verbose(2))
    std::cout << "Number of PEs per root node available: " << nProcsPerNode << "\n";

  quadtree.max_depth = round(log(nProcsPerNode) / log(4));

  if (report.test_verbose(2))
    std::cout << "Quadtree max depth : " << quadtree.max_depth << "\n";

  // Check to see if we have enough processors to do this stuff:
  int nBlocksLonGeo = pow(2, quadtree.max_depth); // input.get_nBlocksLonGeo();
  int nBlocksLatGeo = pow(2, quadtree.max_depth); // input.get_nBlocksLatGeo();
  nGrids = nBlocksLonGeo * nBlocksLatGeo * quadtree.nRootNodes;
  int nProcsNeeded = nMembers * nGrids;

  if (nProcsNeeded == nProcs) {

    // Get Ensemble member number and grid number:
    iMember = iProc / nGrids;
    iGrid = iProc % nGrids;

    if (report.test_verbose(2))
      std::cout << "iProc : " << iProc
                << "; iMember : " << iMember
                << "; iGrid : " << iGrid << "\n";

    // Create strings to allow for easier output filename creation:
    cProc = "p" + tostr(iProc, 4);
    cMember = "m" + tostr(iMember, 4);
    cGrid = "g" + tostr(iGrid, 4);

    // Need to initialize the random number seeds:
    int seed = input.get_original_seed();

    if (seed == 0) {
      // need to generate a real seed and pass it to all processors:
      if (iProc == 0)
        seed = int(std::chrono::system_clock::now().time_since_epoch().count());

      MPI_Bcast(&seed, 1, MPI_INT, 0, aether_comm);
    }

    // Make each seed unique for the ensemble member:
    seed = seed + iMember;
    input.set_seed(seed);

    if (report.test_verbose(2))
      std::cout << "seed : " << seed << "\n";

    quadtree.build(input, report);

  } else {
    if (iProc == 0) {
      std::cout << "Number of processors needed is not set correctly!\n";
      std::cout << "nProcs needs to be equal to the multiplication of :\n";
      std::cout << "   nMembers : " << nMembers << "\n";
      std::cout << "   nBlocksLonGeo : " << nBlocksLonGeo << "\n";
      std::cout << "   nBlocksLatGeo : " << nBlocksLatGeo << "\n";
      std::cout << "      total needed : " << nProcsNeeded << "\n";
      std::cout << "      which is not equal to nProcs : " << nProcs << "\n";
    }

    DidWork = false;
  }

  return DidWork;
}
