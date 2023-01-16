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
  bool DidWork = true;

  Times time;
  Report report;

  // Define the function and report:
  std::string function = "main";
  static int iFunction = -1;
  report.enter(function, iFunction);

  try {
  
    // Create inputs (reading the input file):
    Inputs input(time, report);
    if (!input.is_ok())
      throw std::string("input initialization failed!");

    if (input.get_is_student())
      report.print(-1, "Hello " +
		   input.get_student_name() + " - welcome to Aether!");
    
    Quadtree quadtree(input, report);
    if (!quadtree.is_ok())
      throw std::string("quadtree initialization failed!");
    
    // Initialize MPI and parallel aspects of the code:
    DidWork = init_parallel(input, quadtree, report);
    if (!DidWork)
      throw std::string("init_parallel failed!");
  
    // Everything should be set for the inputs now, so write a restart file:
    DidWork = input.write_restart();
    if (!DidWork)
      throw std::string("input.write_restart failed!");
    
    // Initialize the EUV system:
    Euv euv(input, report);
    if (!euv.is_ok())
      throw std::string("EUV initialization failed!");
    
    // Initialize the planet:
    Planets planet(input, report);
    if (!planet.is_ok())
      throw std::string("planet initialization failed!");

    // Initialize the indices, read the files, and perturb:
    Indices indices(input);
    DidWork = read_and_store_indices(indices, input, report);
    if (!DidWork)
      throw std::string("read_and_store_indices failed!");

    // Perturb the inputs if user has asked for this
    indices.perturb(input, report);
    
    // Initialize Geographic grid:
    Grid gGrid(input.get_nLonsGeo(),
	       input.get_nLatsGeo(),
	       input.get_nAltsGeo(),
	       nGeoGhosts);
    DidWork = gGrid.init_geo_grid(quadtree, planet, input, report);
    if (!DidWork)
      throw std::string("init_geo_grid failed!");
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

    // Initialize ion temperatures from neutral temperature
    ions.init_ion_temperature(neutrals, gGrid, report);

    // Initialize electrodynamics and check if electrodynamics times
    // works with input time
    Electrodynamics electrodynamics(time, input, report);
    if (!electrodynamics.is_ok())
      throw std::string("electrodynamics initialization failed!");

    if (input.get_do_restart()) {
      report.print(1, "Restarting! Reading time file!");
      DidWork = time.restart_file(input.get_restartin_dir(), DoRead);
      if (!DidWork)
	throw std::string("Reading Restart for time Failed!!!\n");
    }

    // This is for the initial output.  If it is not a restart, this will go:
    if (time.check_time_gate(input.get_dt_output(0)))
      iErr = output(neutrals, ions, gGrid, time, planet, input, report);

    // This is advancing now... We are not coupling, so set dt_couple to the
    // end of the simulation

    double dt_couple = time.get_end() - time.get_current();

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

	DidWork = neutrals.restart_file(input.get_restartout_dir(), DoWrite);
	if (!DidWork)
	  throw std::string("Writing Restart for Neutrals Failed!!!\n");	

	DidWork = ions.restart_file(input.get_restartout_dir(), DoWrite);
	if (!DidWork)
	  throw std::string("Writing Restart for Ions Failed!!!\n");	

	DidWork = time.restart_file(input.get_restartout_dir(), DoWrite);
	if (!DidWork)
	  throw std::string("Writing Restart for time Failed!!!\n");	
      }

      // Do some coupling here. But we have no coupling to do. Sad.

    } // End of outer time loop - done with run!


    report.exit(function);
    report.times();

  } catch (std::string error) {
    if (iProc == 0) {
      std::cout << error << "\n";
      std::cout << "---- Must Exit! ----\n";
    }
  }

    
  // End parallel tasks:
  iErr = MPI_Finalize();

  return iErr;
}
