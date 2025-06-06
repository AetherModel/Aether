// Copyright 2023, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md
//
// initial version - A. Ridley - Nov. 1, 2024

#include "aether.h"

// -----------------------------------------------------------------------------
// This is where we will call the different advection schemes
// -----------------------------------------------------------------------------

bool Ions::advect_vertical(Grid grid, Times time) {

  bool didWork = true;

  std::string function = "Ions::advance_vertical";
  static int iFunction = -1;
  report.enter(function, iFunction);

  if (input.get_advection_ions_along() == "rusanov")
    solver_vertical_rusanov(grid, time);
  else {
    std::cout << "Parallel Ion Advection solver not found!\n";
    std::cout << "  ==> Requested : "
              << input.get_advection_ions_along()
              << "\n";
    didWork = false;
  }

  report.exit(function);
  return didWork;
}

