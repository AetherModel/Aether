// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

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
  std::string function = "main";
  static int iFunction = -1;
  report.enter(function, iFunction);

  // Create inputs (reading the input file):
  Inputs input(time, report);

  // Initialize the EUV system:
  Euv euv(input, report);

  // Initialize the planet:
  Planets planet(input, report);

  // Initialize the indices (and read the files):
  Indices indices(input);
  iErr = read_and_store_indices(indices, input, report);

  // Initialize Geographic grid:
  Grid gGrid(input.get_nLonsGeo(),
             input.get_nLatsGeo(),
             input.get_nAltsGeo(),
	     nGeoGhosts);
  gGrid.init_geo_grid(planet, input, report);
  gGrid.fill_grid(planet, report);

  // Initialize Magnetic grid:
  Grid mGrid(nMagLonsG, nMagLatsG, nMagAltsG, nMagGhosts);

  // Initialize Neutrals on geographic grid:
  Neutrals neutrals(gGrid, input, report);

  // Initialize Ions on geographic grid:
  Ions ions(gGrid, input, report);

  // Once EUV, neutrals, and ions have been defined, pair cross sections
  euv.pair_euv(neutrals, ions, report);

  // Initialize Chemical scheme (including reading file):
  Chemistry chemistry(neutrals, ions, input, report);

  // Read in the collision frequencies and other diffusion coefficients:
  read_collision_file(neutrals, ions, input, report);

  // Initialize electrodynamics and check if electrodynamics times
  // works with input time
  Electrodynamics electrodynamics(input, report);
  bool times_are_aligned = electrodynamics.check_times(time.get_current(),
                                                       time.get_end());

  if (!times_are_aligned) {
    iErr = 1;
    std::cout << "Times don't align with electrodynamics file! ";
    std::cout << "Please check this!\n";
    return iErr;
  }

  if (input.get_do_restart()) {
    report.print(1, "Restarting! Reading time file!");
    bool DidWork = time.restart_file(input.get_restartin_dir(), DoRead);
    if (!DidWork)
      std::cout << "Reading Restart for time Failed!!!\n";
  }
  
  // This is for the initial output.  If it is not a restart, this will go:
  if (time.check_time_gate(input.get_dt_output(0)))
    iErr = output(neutrals, ions, gGrid, time, planet, input, report);

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

    // Increment until the intermediate time:
    while (time.get_current() < time.get_intermediate())
      iErr = advance(planet,
                     gGrid,
                     time,
                     euv,
                     neutrals,
                     ions,
                     chemistry,
                     electrodynamics,
                     indices,
                     input,
                     report);

    // Do some coupling here. But we have no coupling to do. Sad.
  }

  report.exit(function);
  report.times();
  return iErr;
}
