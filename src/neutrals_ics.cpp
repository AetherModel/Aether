// Copyright 2023, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md
//
// initial version - A. Ridley - May 27, 2023

#include "aether.h"

// -----------------------------------------------------------------------------
//  Set initial conditions for the neutrals.
//    If user wants to restart, then read the restart files.
//    Otherwise, initialize the neutrals.
//
//    Two methods of initialization implemented so far:
//      - Planet: Use fixed density values in the planet.in file and the
//                temperature profile to set the densities and temperature.
//                Densities are filled with hydrostatic solution.
//      - Msis: Use NRL MSIS to set the densities and temperatures.  If the
//              densities are not found, then set to density in planet.in
//              file and fill with hydrostatic.
// -----------------------------------------------------------------------------

bool Neutrals::initial_conditions(Grid &grid,
                                  Times time,
                                  Indices indices) {

  std::string function = "Neutrals::initial_conditions";
  static int iFunction = -1;
  report.enter(function, iFunction);

  // initialize didWork to false, so that we can catch if the initial
  // conditions are not actually set
  bool didWork = false;

  int64_t iLon, iLat, iAlt, iA;
  precision_t alt, r;
  int64_t nAlts = grid.get_nZ(true);
  int64_t nGCs = grid.get_nGCs();
  int64_t nLons = grid.get_nLons();
  int64_t nLats = grid.get_nLats();

  report.print(3, "Creating Neutrals initial_condition");

  // ----------------------------------------------------------
  // Restart file:

  if (input.get_do_restart()) {
    report.print(1, "Restarting! Reading neutral files!");
    didWork = restart_file(input.get_restartin_dir(), grid.get_gridtype(), DoRead);

    if (!didWork)
      report.error("Reading Restart for Neutrals Failed!!!");
  } else {

    json ics = input.get_initial_condition_types();
    std::string icsType = mklower(ics["type"]);

    // ----------------------------------------------------------
    // MSIS:

    if (icsType == "msis") {
      report.print(2, "Using MSIS for Initial Conditions");

      Msis msis;

      if (!msis.is_ok()) {
        didWork = false;
        report.error("MSIS Initial Conditions asked for, ");
        report.error("but MSIS is not compiled/working!");
      } else {
        didWork = true;
        msis.set_time(time);
        precision_t f107 = indices.get_f107(time.get_current());
        precision_t f107a = indices.get_f107a(time.get_current());
        msis.set_f107(f107, f107a);
        msis.set_ap(10.0);
        msis.set_locations(grid.geoLon_scgc,
                           grid.geoLat_scgc,
                           grid.geoAlt_scgc);

        // This is just to check if MSIS is actually working:
        if (msis.is_valid_species("Tn"))
          // if it is, fill will temperature:
          temperature_scgc = msis.get_cube("Tn");
        else
          // if it is not, then fill with a value:
          temperature_scgc.fill(initial_temperatures[0]);

        for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
          if (report.test_verbose(3))
            std::cout << "Setting Species : "
                      << species[iSpecies].cName << "\n";

          if (msis.is_valid_species(species[iSpecies].cName)) {
            if (report.test_verbose(3))
              std::cout << "  Found in MSIS!\n";

            species[iSpecies].density_scgc =
              msis.get_cube(species[iSpecies].cName);
          } else {
            if (report.test_verbose(3))
              std::cout << "  NOT Found in MSIS - setting constant\n";

            species[iSpecies].density_scgc.slice(0).
            fill(species[iSpecies].lower_bc_density);
            fill_with_hydrostatic(iSpecies, 1, nAlts, grid);
          }

        } // for species
      } // msis init worked ok
    } // type = msis

    // ----------------------------------------------------------
    // Planet:

    if (icsType == "planet") {
      report.print(2, "Using planet for Initial Conditions");

      didWork = true;

      // ---------------------------------------------------------------------
      // This section assumes we want a hydrostatic solution given the
      // temperature profile in the planet.in file.
      // ---------------------------------------------------------------------

      // Let's assume that the altitudes are dependent on lat/lon:

      arma_vec alt1d(nAlts);
      arma_vec temp1d(nAlts);

      if (nInitial_temps > 0) {
        for (iLon = 0; iLon < nLons; iLon++) {
          for (iLat = 0; iLat < nLats; iLat++) {
            alt1d = grid.geoAlt_scgc.tube(iLon, iLat);

            for (iAlt = 0; iAlt < nAlts; iAlt++) {
              alt = alt1d(iAlt);

              // Find temperatures:
              if (alt <= initial_altitudes[0])
                temp1d[iAlt] = initial_temperatures[0];

              else {
                if (alt >= initial_altitudes[nInitial_temps - 1])
                  temp1d[iAlt] = initial_temperatures[nInitial_temps - 1];

                else {
                  // Linear interpolation!
                  iA = 0;

                  while (alt > initial_altitudes[iA])
                    iA++;

                  iA--;
                  // alt will be between iA and iA+1
                  r = (alt - initial_altitudes[iA]) /
                      (initial_altitudes[iA + 1] - initial_altitudes[iA]);
                  temp1d[iAlt] =
                    (1.0 - r) * initial_temperatures[iA] +
                    (r) * initial_temperatures[iA + 1];
                }
              }
            }

            temperature_scgc.tube(iLon, iLat) = temp1d;
          }
        }
      } else
        temperature_scgc.fill(200.0);

      // Make the initial condition in the lower ghost cells to be consistent
      // with the actual lower BC:
      for (iLon = 0; iLon < nLons; iLon ++) {
        for (iLat = 0; iLat < nLats; iLat++) {
          for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {

            species[iSpecies].density_scgc.subcube(
              iLon, iLat, 0, iLon, iLat, grid.first_lower_gc(iLon, iLat) + 1).fill(
                species[iSpecies].lower_bc_density);

          }
        }
      }

      report.print(2, "Calculating scale height");
      calc_scale_height(grid);
      report.print(2, "setting lower BCs");
      set_lower_bcs(grid, time, indices);
      report.print(2, "Filling with hydrostatic");

      for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++)
        fill_with_hydrostatic(iSpecies, nGCs - 1, nAlts, grid);
    } // type = planet
  }

  /*
  This section is for putting an initial blob into the simulation
  to test the advection solver.
  precision_t lon_0 = 0.0;
  precision_t lat_0 = 0.0;
  precision_t r_0 = 150.0 * 1000.0 * 10.0;

  for (iAlt = 0; iAlt < nAlts; iAlt++) {

    for (int64_t iLat = 0; iLat < nLats; iLat++) {
      for (int64_t iLon = 0; iLon < nLons; iLon++) {
        precision_t curr_lat = grid.geoLat_scgc(iLon, iLat, iAlt);
        precision_t curr_lon = grid.geoLon_scgc(iLon, iLat, iAlt);
        precision_t R = grid.radius_scgc(iLon, iLat, iAlt);

        // Calculate great circle distance
        precision_t dlon_2 = (curr_lon - lon_0) / 2.0;
        precision_t dlat_2 = (curr_lat - lat_0) / 2.0;

        precision_t r_d = 2.0 * R * asin(sqrt(sin(dlat_2) * sin(dlat_2) + sin(
                                                dlon_2) * sin(dlon_2) * cos(curr_lat) * cos(lat_0)));

        if (r_d < r_0) {
          for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
            species[iSpecies].density_scgc(iLon, iLat,
                                           iAlt) = species[iSpecies].density_scgc(iLon, iLat, iAlt) * 10.;
            std::cout << "increasing density!\n";
          }
        }
      }
    }
  }
  */


  // ensure that the densities are all within bounds:
  clamp_density();

  if (!didWork)
    report.error("Issue with initial conditions!");

  report.exit(function);
  return didWork;
}

bool Neutrals::cosine_bell_ic(Grid grid,
                              Times time,
                              Indices indices,
                              Planets planet) {
  std::string function = "Neutrals::cosine_bell_ic";
  static int iFunction = -1;
  report.enter(function, iFunction);

  // Planet.get_radius() takes in latitude
  // but at current stage is unimplemented
  // Anyway, we use equator radius as assumption for CubeSphere
  // CubeSphere must be a perfect sphere!!
  precision_t planet_R = planet.get_radius(0);

  // radius of planet + altitude
  // just pick alt at (0,0) loction
  arma_vec R_Alts = grid.geoAlt_scgc.tube(0, 0) + planet_R;


  /** Get a bunch of constants for setting up the ic **/
  precision_t R = R_Alts(2); // select R in the middle
  //7.37128e+06 meters Earth radius + 1000km height
  //std::cout << R << std::endl;

  // Determine flow direction
  // 0 - Equatorial
  // pi/2 - Meridional
  // pi/4 - NE direction
  precision_t alpha_0 = 0.;

  // Scaling factor for physical velocity
  // 12 day period in miliseconds
  precision_t u_0 = cTWOPI * R / (12.*24.*60.*60.);

  // Radius of the cosine bell
  precision_t r_0 = R / 3.;

  // Center of the cosine bell
  precision_t lon_0 = 7.*cPI / 4. + cPI / 8;
  precision_t lat_0 = 0.;

  // Maximum height for the cosine bell
  precision_t h_0 = 1000.;

  // Some grid dimensions and coordinates
  int64_t nLats = grid.get_nLats();
  int64_t nLons = grid.get_nLons();
  int64_t nAlts = grid.get_nAlts();
  arma_mat lat_grid = grid.geoLat_scgc.slice(2);
  arma_mat lon_grid = grid.geoLon_scgc.slice(2);

  // Calculate for physical velocity for every altitude
  // First we prepare velocities for one slice
  arma_mat slice_u = velocity_vcgc[0].slice(2);
  arma_mat slice_v = velocity_vcgc[1].slice(2);

  // Fill velocities in one slice
  for (int64_t iLat = 0; iLat < nLats; iLat++) {
    for (int64_t iLon = 0; iLon < nLons; iLon++) {
      precision_t curr_lat = lat_grid(iLon, iLat);
      precision_t curr_lon = lon_grid(iLon, iLat);
      slice_u(iLon, iLat) = u_0 * (cos(alpha_0) * cos(curr_lat) + sin(alpha_0) * cos(
                                     curr_lon) * sin(curr_lat));
      slice_v(iLon, iLat) = -u_0 * sin(alpha_0) * sin(curr_lon);
    }
  }

  // Update this slice of velocity to all slices (for completeness)
  for (int64_t iAlt = 0; iAlt < nAlts; iAlt++) {
    velocity_vcgc[0].slice(iAlt) = slice_u;
    velocity_vcgc[1].slice(iAlt) = slice_v;
  }

  // Calculate the cosine bell or rho_scgc
  // First, again take a slice
  arma_mat slice_rho = rho_scgc.slice(2);

  // Fill rho in one slice
  for (int64_t iLat = 0; iLat < nLats; iLat++) {
    for (int64_t iLon = 0; iLon < nLons; iLon++) {
      precision_t curr_lat = lat_grid(iLon, iLat);
      precision_t curr_lon = lon_grid(iLon, iLat);

      // Calculate great circle distance
      precision_t dlon_2 = (curr_lon - lon_0) / 2.0;
      precision_t dlat_2 = (curr_lat - lat_0) / 2.0;

      precision_t r_d = 2.0 * R * asin(sqrt(sin(dlat_2) * sin(dlat_2) + sin(
                                              dlon_2) * sin(dlon_2) * cos(curr_lat) * cos(lat_0)));

      if (r_d < r_0)
        slice_rho(iLon, iLat) = (h_0 / 2) * (1 + cos(cPI * r_d / r_0));

      else
        slice_rho(iLon, iLat) = 0.;
    }
  }

  // Update this slice of rho to all slices (for completeness)
  for (int64_t iAlt = 0; iAlt < nAlts; iAlt++) {
    rho_scgc.slice(iAlt) = slice_rho;

    // Do zero concentration conversion
    for (int64_t iSpec = 0; iSpec < nSpecies; iSpec++)
      species[iSpec].density_scgc.slice(iAlt) = slice_rho / species[iSpec].mass;
  }

  // Add some velocity pertubation
  //std::cout <<  velocity_vcgc[0].slice(2) << std::endl;
  //std::cout << rho_scgc.slice(2) << std::endl;

  return 1;
}

bool Neutrals::blob_ic(Grid grid,
                       Times time,
                       Indices indices,
                       Planets planet) {
  std::string function = "Neutrals::blob_ic";
  static int iFunction = -1;
  report.enter(function, iFunction);

  // Planet.get_radius() takes in latitude
  // but at current stage is unimplemented
  // Anyway, we use equator radius as assumption for CubeSphere
  // CubeSphere must be a perfect sphere!!
  precision_t planet_R = planet.get_radius(0);

  // radius of planet + altitude
  // just pick alt at (0,0) loction
  arma_vec R_Alts = grid.geoAlt_scgc.tube(0, 0) + planet_R;


  /** Get a bunch of constants for setting up the ic **/
  precision_t R = R_Alts(2); // select R in the middle
  //7.37128e+06 meters Earth radius + 1000km height
  //std::cout << R << std::endl;

  // Radius of the blob
  // Hardcoded
  precision_t r_0 = 111321 * 10;

  // Center of the blob
  precision_t lon_0 = 7.*cPI / 4. - cPI / 8.;
  precision_t lat_0 = 0.;

  // Some grid dimensions and coordinates
  int64_t nLats = grid.get_nLats();
  int64_t nLons = grid.get_nLons();
  int64_t nAlts = grid.get_nAlts();
  arma_mat lat_grid = grid.geoLat_scgc.slice(2);
  arma_mat lon_grid = grid.geoLon_scgc.slice(2);

  // Calculate for physical velocity for every altitude
  // First we prepare velocities for one slice
  arma_mat slice_u = velocity_vcgc[0].slice(2);
  arma_mat slice_v = velocity_vcgc[1].slice(2);

  // Fill velocities in one slice
  for (int64_t iLat = 0; iLat < nLats; iLat++) {
    for (int64_t iLon = 0; iLon < nLons; iLon++) {
      slice_u(iLon, iLat) = 0.;
      slice_v(iLon, iLat) = 0.;
    }
  }

  // Update this slice of velocity to all slices (for completeness)
  for (int64_t iAlt = 0; iAlt < nAlts; iAlt++) {
    velocity_vcgc[0].slice(iAlt) = slice_u;
    velocity_vcgc[1].slice(iAlt) = slice_v;
  }

  // Calculate the cosine bell or rho_scgc
  // First, again take a slice
  arma_mat slice_rho = rho_scgc.slice(2);

  // Fill rho in one slice
  for (int64_t iLat = 0; iLat < nLats; iLat++) {
    for (int64_t iLon = 0; iLon < nLons; iLon++) {
      precision_t curr_lat = lat_grid(iLon, iLat);
      precision_t curr_lon = lon_grid(iLon, iLat);

      // Calculate great circle distance
      precision_t dlon_2 = (curr_lon - lon_0) / 2.0;
      precision_t dlat_2 = (curr_lat - lat_0) / 2.0;

      precision_t r_d = 2.0 * R * asin(sqrt(sin(dlat_2) * sin(dlat_2) + sin(
                                              dlon_2) * sin(dlon_2) * cos(curr_lat) * cos(lat_0)));

      if (r_d < r_0)
        slice_rho(iLon, iLat) = 5e-12;

      else
        slice_rho(iLon, iLat) = 1e-12;
    }
  }

  // Update this slice of rho to all slices (for completeness)
  for (int64_t iAlt = 0; iAlt < nAlts; iAlt++) {
    rho_scgc.slice(iAlt) = slice_rho;

    // Do zero concentration conversion
    for (int64_t iSpec = 0; iSpec < nSpecies; iSpec++)
      species[iSpec].density_scgc.slice(iAlt) = slice_rho / species[iSpec].mass;
  }

  // Temperature setup
  for (int64_t iAlt = 0; iAlt < nAlts; iAlt++)
    temperature_scgc.slice(iAlt) = 600.*arma_mat(nLons, nLats, fill::ones);

  // Add some velocity pertubation
  //std::cout <<  velocity_vcgc[0].slice(2) << std::endl;
  //std::cout << rho_scgc.slice(2) << std::endl;

  return 1;
}


