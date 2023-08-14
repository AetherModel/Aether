// Copyright 2023, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md
//
// initial version - A. Ridley - July 28, 2023

#include "aether.h"

// -----------------------------------------------------------------------------
// This is where we will call the different advection schemes
// -----------------------------------------------------------------------------

bool Neutrals::advect_vertical(Grid grid, Times time) {

  bool didWork = true;
  
  std::string function = "Neutrals::advance_vertical";
  static int iFunction = -1;
  report.enter(function, iFunction);

  if (input.get_advection_neutrals_vertical() == "hydro")
    fill_with_hydrostatic(1, grid.get_nZ(), grid);
  else if (input.get_advection_neutrals_vertical() == "rusanov")
    solver_vertical_rusanov(grid, time);
  else {
    std::cout << "Vertical solver not found!\n";
    std::cout << "  ==> Requested : "
	      << input.get_advection_neutrals_vertical()
	      << "\n";
    didWork = false;
  }
  report.exit(function);
  return didWork;
}

