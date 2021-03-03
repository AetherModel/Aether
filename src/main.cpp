// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>

#include "aether.h"

int main() {

  int iErr = 0;

  Times time;
  Report report;

  std::string function = "main";
  static int iFunction = -1;
  report.enter(function, iFunction);

  Inputs input(time, report);
  Euv euv(input, report);
  Planets planet(input, report);
  Indices indices(input);
  iErr = read_and_store_indices(indices, input, report);
  
  // Geo grid stuff:
  Grid gGrid(input.get_nLonsGeo(),
	     input.get_nLatsGeo(),
	     input.get_nAltsGeo(), nGeoGhosts);
  gGrid.init_geo_grid(planet, input, report);
  gGrid.fill_grid(planet, report);

  // Magnetic grid stuff:
  Grid mGrid(nMagLonsG, nMagLatsG, nMagAltsG, nMagGhosts);

  Neutrals neutrals(gGrid, input, report);
  Ions ions(gGrid, input, report);
  euv.pair_euv(neutrals, ions, report);

  Chemistry chemistry(neutrals, ions, input, report);

  // This is for the initial output.  If it is not a restart, this will go:
  if (time.check_time_gate(input.get_dt_output(0))) {
    iErr = output(neutrals, ions, gGrid, time, planet, input, report);
  }

  // This is advancing now...

  double dt_couple = 1800.0;

  // The way most codes are set up in the SWMF is that there are two
  // times, an end time which ends the simulation, and an intermediate
  // time, which allows coupling or something to happen.  So, typically
  // the advance functions should only go to this intermediate time,
  // then a loop around that goes to the end time.  Then, the code can
  // be made into a library and run externally.

  while (time.get_current() < time.get_end()) {

    time.increment_intermediate(dt_couple);

    while (time.get_current() < time.get_intermediate())
      iErr = advance(planet,
                     gGrid,
                     time,
                     euv,
                     neutrals,
                     ions,
                     chemistry,
                     indices,
                     input,
                     report);

    // Do some coupling here. But we have no coupling to do. Sad.
  }
  report.exit(function);
  report.times();
  return iErr;
}
