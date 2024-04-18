// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "../include/aether.h"
#include <cassert>

// Run with the following settings:
// IsCubeSphere    nMembers      mpi -np _
//     True           1              6
//     False          1              4
//     False          2              8

// -----------------------------------------------------------------------------
// First make copies of the physical cells of gGrid.geoLon_scgc, geoLat_scgc and
// geoAlt.scgc, and set the copies' ghost cells to be -inf
// 
// Then exchange the message
//
// Finally compare the copies' ghost cells and originals' ghost cells to see if
// the messages get exchagned correctly
// -----------------------------------------------------------------------------

int64_t nGCs;
int64_t nLons;
int64_t nLats;
int64_t nAlts;

void wrap_point(precision_t &lon, precision_t &lat) {
    int inverse = 1;
    if (lat < -0.5 * cPI) {
        lat = -cPI - lat;
        inverse = -1;
        lon += cPI;
    } else if (lat > 0.5 * cPI) {
        lat = cPI - lat;
        inverse = -1;
        lon += cPI;
    }
    while (lon < 0) {
        lon += cTWOPI;
    }
    while (lon > cTWOPI) {
        lon -= cTWOPI;
    }
}

// Return true if the inputs are checked, false they are ignored
// assert(false) if the check fails
bool check(precision_t lon, precision_t lat,
           precision_t est_lon, precision_t est_lat) {
    wrap_point(lon, lat);
    // Don't check when lon is at the 360-0 hop
    if (lon < cPI / 1.5 / nLons || lon > cTWOPI - cPI / 1.5 / nLons) {
        // std::cout << (std::abs((est_lon - lon) / lon) < 0.05) << '\n';
        return false;
    }
    assert(std::abs((est_lon - lon) / lon) < 0.05);
    if (lat == 0) {
        assert(std::abs(lat) < 0.001);
    } else {
        assert(std::abs((est_lat - lat) / lat) < 0.05);
    }
    return true;
}

int main() {
  // This test needs precision_t to be double
  assert(sizeof(precision_t) == sizeof(double));

  int iErr = 0;
  std::string sError;
  bool DidWork = true;

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
    DidWork = init_parallel(quadtree);
    if (!DidWork)
      throw std::string("init_parallel failed!");

    // Everything should be set for the inputs now, so write a restart file:
    DidWork = input.write_restart();
    if (!DidWork)
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
    DidWork = read_and_store_indices(indices);
    MPI_Barrier(aether_comm);
    if (!DidWork)
      throw std::string("read_and_store_indices failed!");
    
    // Perturb the inputs if user has asked for this
    indices.perturb();
    MPI_Barrier(aether_comm);
    
    // Initialize Geographic grid:
    Grid gGrid(input.get_nLonsGeo(),
	       input.get_nLatsGeo(),
	       input.get_nAltsGeo(),
	       nGeoGhosts);
    DidWork = gGrid.init_geo_grid(quadtree, planet);
    MPI_Barrier(aether_comm);
    if (!DidWork)
      throw std::string("init_geo_grid failed!");  

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
      DidWork = neutrals.check_for_nonfinites();
      DidWork = ions.check_for_nonfinites();
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
      DidWork = time.restart_file(input.get_restartin_dir(), DoRead);
      if (!DidWork)
	throw std::string("Reading Restart for time Failed!!!\n");
    }

    // This is for the initial output.  If it is not a restart, this will go:
    if (time.check_time_gate(input.get_dt_output(0)))
      iErr = output(neutrals, ions, gGrid, time, planet);

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


    nGCs = gGrid.get_nGCs();
    nLons = gGrid.get_nLons();
    nLats = gGrid.get_nLats();
    nAlts = gGrid.get_nAlts();
    // Initialize arma cubes with -inf
    arma_cube testLon(nLons, nLats, nAlts, fill::value(cNinf));
    arma_cube testLat(nLons, nLats, nAlts, fill::value(cNinf));
    arma_cube testAlt(nLons, nLats, nAlts, fill::value(cNinf));

    // Make copies of physical cells and the innermost ghost cells
    testLon.subcube(nGCs - 1,
                    nGCs - 1,
                    nGCs - 1,
                    nLons - nGCs,
                    nLats - nGCs,
                    nAlts - nGCs)
        = gGrid.geoLon_scgc.subcube(nGCs - 1,
                                    nGCs - 1,
                                    nGCs - 1,
                                    nLons - nGCs,
                                    nLats - nGCs,
                                    nAlts - nGCs);
    testLat.subcube(nGCs - 1,
                    nGCs - 1,
                    nGCs - 1,
                    nLons - nGCs,
                    nLats - nGCs,
                    nAlts - nGCs)
        = gGrid.geoLat_scgc.subcube(nGCs - 1,
                                    nGCs - 1,
                                    nGCs - 1,
                                    nLons - nGCs,
                                    nLats - nGCs,
                                    nAlts - nGCs);
    testAlt.subcube(nGCs - 1,
                    nGCs - 1,
                    nGCs - 1,
                    nLons - nGCs,
                    nLats - nGCs,
                    nAlts - nGCs)
        = gGrid.geoAlt_scgc.subcube(nGCs - 1,
                                    nGCs - 1,
                                    nGCs - 1,
                                    nLons - nGCs,
                                    nLats - nGCs,
                                    nAlts - nGCs);

    // They are all scalars rather than vectors, so do not inverse
    // when crossing the pole
    gGrid.exchange(testLon, false);
    gGrid.exchange(testLat, false);
    gGrid.exchange(testAlt, false);
    
    for (int64_t iAlt = nGCs; iAlt < nGCs + 1; ++iAlt) {
        for (int64_t iLat = 0; iLat < nLats; ++iLat) {
            for (int64_t iLon = 0; iLon < nLons; ++iLon) {
                // Check
                check(gGrid.geoLon_scgc(iLon, iLat, iAlt),
                      gGrid.geoLat_scgc(iLon, iLat, iAlt),
                      testLon(iLon, iLat, iAlt),
                      testLat(iLon, iLat, iAlt));
            }
        }
    }

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
