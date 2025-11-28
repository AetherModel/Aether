// Copyright 2023, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md
//
// initial version - A. Ridley - July 28, 2023

#include "aether.h"

// -----------------------------------------------------------------------------
// This is where we will call the different advection schemes
// -----------------------------------------------------------------------------

bool Neutrals::advect_vertical(Grid &grid, Times time) {

  bool didWork = true;

  std::string function = "Neutrals::advance_vertical";
  static int iFunction = -1;
  report.enter(function, iFunction);

  if (grid.get_IsDipole() ||
      input.get_advection_neutrals_vertical() == "hydro")
    fill_with_hydrostatic(0, grid.get_nZ(), grid);

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

bool Neutrals::advect_horizontal(Grid & grid, Times & time) {
  bool didWork = true;

  std::string function = "Neutrals::advance_horizontal";
  static int iFunction = -1;
  report.enter(function, iFunction);

  if (grid.iGridShape_ == iCubesphere_) {
    if (input.get_advection_neutrals_horizontal() == "advect_test")
      solver_horizontal_RK4_advection(grid, time);
    else if (input.get_advection_neutrals_horizontal() == "fv")
      solver_horizontal_RK1_rochi(grid, time);

    else {
      std::cout << "Horizontal solver not found!\n";
      std::cout << "  ==> Requested : "
                << input.get_advection_neutrals_horizontal()
                << "\n";
      didWork = false;
    }
  } else
    advect_sphere(grid, time);


  report.exit(function);
  return didWork;
}

bool Neutrals::advect_horizontal_advection(Grid & grid, Times & time) {
  bool didWork = true;

  std::string function = "Neutrals::advance_horizontal_advection";
  static int iFunction = -1;
  report.enter(function, iFunction);

  //solver_horizontal_rusanov_advection(grid, time);
  solver_horizontal_RK4_advection(grid, time);

  report.exit(function);
  return didWork;
}
