// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "../include/aether.h"

// -----------------------------------------------------------------------------
// Main file for the Aether model.  This is needed when Aether is not used
// as a library in another code, such as the SWMF.
// -----------------------------------------------------------------------------

int main() {

  int iErr = 0;
  std::string sError;
  bool didWork = true;

  Times time;

  // Define the function and report:
  std::string function = "main";
  static int iFunction = -1;
  report.enter(function, iFunction);

  try {
    // Create inputs (reading the input file):
    input = Inputs(time);

    if (!input.is_ok())
      throw std::string("input initialization failed!");

    if (input.get_is_student())
      report.print(-1, "Hello " +
                   input.get_student_name() + " - welcome to Aether!");

    Quadtree quadtree;

    if (!quadtree.is_ok())
      throw std::string("quadtree initialization failed!");

    // Initialize MPI and parallel aspects of the code:
    didWork = init_parallel(quadtree);

    if (!didWork)
      throw std::string("init_parallel failed!");

    // Everything should be set for the inputs now, so write a restart file:
    didWork = input.write_restart();

    if (!didWork)
      throw std::string("input.write_restart failed!");

    // Initialize the EUV system:
    Euv euv;

    if (!euv.is_ok())
      throw std::string("EUV initialization failed!");

    // Initialize the planet:
    Planets planet;
    MPI_Barrier(aether_comm);

    if (!planet.is_ok())
      throw std::string("planet initialization failed!");

    // Initialize the indices, read the files, and perturb:
    Indices indices;
    didWork = read_and_store_indices(indices);
    MPI_Barrier(aether_comm);

    if (!didWork)
      throw std::string("read_and_store_indices failed!");

    // Perturb the inputs if user has asked for this
    indices.perturb();
    MPI_Barrier(aether_comm);

    // Initialize Geographic grid:
    Grid gGrid(input.get_nLonsGeo(),
               input.get_nLatsGeo(),
               input.get_nAltsGeo(),
               nGeoGhosts);
    didWork = gGrid.init_geo_grid(quadtree, planet);
    MPI_Barrier(aether_comm);

    if (!didWork)
      throw std::string("init_geo_grid failed!");

    // Find interpolation coefs for the ghostcells if cubesphere grid
    didWork = find_ghostcell_interpolation_coefs(gGrid);

    // Calculate centripetal acceleration, since this is a constant
    // vector on the grid:
    if (input.get_cent_acc())
      gGrid.calc_cent_acc(planet);

    // Initialize Magnetic grid:
    Grid mGrid(nMagLonsG, nMagLatsG, nMagAltsG, nMagGhosts);

    // Initialize Neutrals on geographic grid:
    Neutrals neutrals(gGrid, planet, time, indices);

    // Initialize Ions on geographic grid:
    Ions ions(gGrid, planet);

    // -----------------------------------------------------------------
    // This is a unit test for checking for nans and infinities.
    // Is simply adds nans and infinities in a few places, then
    // checks for them to make sure the checking is working
    // -----------------------------------------------------------------

    if (input.get_nan_test()) {
      neutrals.nan_test(input.get_nan_test_variable());
      ions.nan_test(input.get_nan_test_variable());
    }

    if (input.get_check_for_nans()) {
      didWork = neutrals.check_for_nonfinites();
      didWork = ions.check_for_nonfinites();
    }

    // -----------------------------------------------------------------

    // Once EUV, neutrals, and ions have been defined, pair cross sections
    euv.pair_euv(neutrals, ions);

    // Initialize Chemical scheme (including reading file):
    Chemistry chemistry(neutrals, ions);

    // Read in the collision frequencies and other diffusion coefficients:
    read_collision_file(neutrals, ions);

    // Initialize ion temperatures from neutral temperature
    ions.init_ion_temperature(neutrals, gGrid);

    // Initialize electrodynamics and check if electrodynamics times
    // works with input time
    Electrodynamics electrodynamics(time);

    if (!electrodynamics.is_ok())
      throw std::string("electrodynamics initialization failed!");

    // If the user wants to restart, then get the time of the restart
    if (input.get_do_restart()) {
      report.print(1, "Restarting! Reading time file!");
      didWork = time.restart_file(input.get_restartin_dir(), DoRead);

      if (!didWork)
        throw std::string("Reading Restart for time Failed!!!\n");
    }

    // This is for the initial output.  If it is not a restart, this will go:
    if (time.check_time_gate(input.get_dt_output(0)))
      didWork = output(neutrals, ions, gGrid, time, planet);

    if (!didWork)
      throw std::string("output failed!");

    // This is advancing now... We are not coupling, so set dt_couple to the
    // end of the simulation

    double dt_couple = time.get_end() - time.get_current();

    // The way most codes are set up in the SWMF is that there are two
    // times, an end time which ends the simulation, and an intermediate
    // time, which allows coupling or something to happen.  So, typically
    // the advance functions should only go to this intermediate time,
    // then a loop around that goes to the end time.  Then, the code can
    // be made into a library and run externally.

    Logfile logfile(indices);

    time.set_start_time_loop();

    while (time.get_current() < time.get_end()) {

      time.increment_intermediate(dt_couple);

      // Increment until the intermediate time:
      while (time.get_current() < time.get_intermediate()) {
        didWork = advance(planet,
                          gGrid,
                          time,
                          euv,
                          neutrals,
                          ions,
                          chemistry,
                          electrodynamics,
                          indices,
                          logfile);

        if (!didWork)
          throw std::string("Error in advance!");
      }

      // Should write out some restart files every time we are done with
      // intermediate times.  Just so when we restart, we know that we can
      // couple first thing and everything should be good. (Not sure if
      // restart should be before or after the coupling, but since we are
      // not coupling, it doesn't matter.  Once we do coupling to something,
      // need to figure it out.
      //
      // The odd thing here is that in advance, we most likely JUST
      // wrote out restart files, so we only need to do this if we
      // didn't just do it.  So, check the negative here:
      if (!time.check_time_gate(input.get_dt_write_restarts())) {
        report.print(3, "Writing restart files");

        didWork = neutrals.restart_file(input.get_restartout_dir(), DoWrite);

        if (!didWork)
          throw std::string("Writing Restart for Neutrals Failed!!!\n");

        didWork = ions.restart_file(input.get_restartout_dir(), DoWrite);

        if (!didWork)
          throw std::string("Writing Restart for Ions Failed!!!\n");

        didWork = time.restart_file(input.get_restartout_dir(), DoWrite);

        if (!didWork)
          throw std::string("Writing Restart for time Failed!!!\n");
      }

      // Do some coupling here. But we have no coupling to do. Sad.

    } // End of outer time loop - done with run!

    report.exit(function);
    report.times();

  } catch (std::string error) {
    report.report_errors();

    if (iProc == 0) {
      std::cout << error << "\n";
      std::cout << "---- Must Exit! ----\n";
    }
  }


  // End parallel tasks:
  iErr = MPI_Finalize();

  return iErr;
}
