// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>

#include "../include/aether.h"

// -----------------------------------------------------------------------------
// Main file for the Aether model.  This is needed when Aether is not used
// as a library in another code, such as the SWMF.
// -----------------------------------------------------------------------------

int main() {

  int iErr = 0;

  Times time;
  Report report;

  // Define the function and report:
  std::string function = "test_mgrid";
  static int iFunction = -1;
  report.enter(function, iFunction);

  // Create inputs (reading the input file):
  Inputs input(time, report);

  // Initialize the EUV system:
  // Euv euv(input, report);

  // Initialize the planet:
  Planets planet(input, report);

  // Initialize the indices (and read the files):
  Indices indices(input);
  iErr = read_and_store_indices(indices, input, report);

  // Initialize Geographic grid:
  Grid gGrid(input.get_nLonsGeo(),
       input.get_nLatsGeo(),
       input.get_nAltsGeo(), nGeoGhosts);
  gGrid.init_geo_grid(planet, input, report);
  gGrid.fill_grid(planet, report);

  // Initialize Magnetic grid:
  cout << "input.get_nAltsMag() " <<input.get_nAltsMag()<<endl;
  Grid mGrid(input.get_nLonsMag(),
       input.get_nLatsMag(),
       input.get_nAltsMag(),nMagGhosts);
  mGrid.init_mag_grid(planet, input, report);

  report.exit(function);
  report.times();
  return iErr;
}
