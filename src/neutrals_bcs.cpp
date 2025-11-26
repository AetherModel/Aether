// Copyright 2023, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md
//
// initial version - A. Ridley - May 27, 2023

#include "aether.h"

// -----------------------------------------------------------------------------
//  Set initial conditions for the neutrals.
//    Two methods implemented so far:
//      - Planet: Use fixed density values in the planet.in file and the
//                temperature profile to set the densities and temperature.
//                Densities are filled with hydrostatic solution.
//      - Msis: Use NRL MSIS to set the densities and temperatures.  If the
//              densities are not found, then set to density in planet.in
//              file and fill with hydrostatic.
// -----------------------------------------------------------------------------

//----------------------------------------------------------------------
// set_bcs - This is for setting the vertical BCs
//----------------------------------------------------------------------

bool Neutrals::set_bcs(Grid &grid,
                       Times time,
                       Indices indices) {

  std::string function = "Neutrals::set_bcs";
  static int iFunction = -1;
  report.enter(function, iFunction);

  bool didWork = true;

  if (grid.get_nAlts(false) > 1) {
    didWork = set_lower_bcs(grid, time, indices);

    if (didWork)
      didWork = set_upper_bcs(grid);

    if (didWork)
      calc_mass_density();
  }

  if (!didWork)
    report.error("issue with BCs!");

  report.exit(function);
  return didWork;
}

//----------------------------------------------------------------------
// set lower boundary conditions for the neutrals
//----------------------------------------------------------------------

bool Neutrals::set_upper_bcs(Grid &grid) {

  std::string function = "Neutrals::set_upper_bcs";
  static int iFunction = -1;
  report.enter(function, iFunction);

  bool didWork = true;

  int64_t nAlts = grid.get_nZ();
  int64_t nX = grid.get_nX(), iX;
  int64_t nY = grid.get_nY(), iY;
  int64_t nGCs = grid.get_nGCs();
  int64_t iAlt;
  arma_mat h;

  for (iAlt = nAlts - nGCs; iAlt < nAlts; iAlt++) {

    // Bulk Quantities:
    temperature_scgc.slice(iAlt) = temperature_scgc.slice(iAlt - 1);
    velocity_vcgc[0].slice(iAlt) = velocity_vcgc[0].slice(iAlt - 1);
    velocity_vcgc[1].slice(iAlt) = velocity_vcgc[1].slice(iAlt - 1);

    // For each species:
    for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      // Horizontal velocities - zero gradient:
      species[iSpecies].velocity_vcgc[0].slice(iAlt) =
        species[iSpecies].velocity_vcgc[0].slice(iAlt - 1);
      species[iSpecies].velocity_vcgc[1].slice(iAlt) =
        species[iSpecies].velocity_vcgc[1].slice(iAlt - 1);

      // Allow upflow, but not downflow:
      for (iX = nGCs; iX < nX - nGCs; iX++)
        for (iY = nGCs; iY < nY - nGCs; iY++)
          if (species[iSpecies].velocity_vcgc[2](iX, iY, iAlt - 1) > 0)
            species[iSpecies].velocity_vcgc[2](iX, iY, iAlt) =
              species[iSpecies].velocity_vcgc[2](iX, iY, iAlt - 1);
          else
            species[iSpecies].velocity_vcgc[2](iX, iY, iAlt) = 0.0;

      h = species[iSpecies].scale_height_scgc.slice(iAlt);
      species[iSpecies].density_scgc.slice(iAlt) =
        species[iSpecies].density_scgc.slice(iAlt - 1) %
        exp(-grid.dk_edge_m.slice(iAlt) / h);
    }
  }

  report.exit(function);
  return didWork;
}

//----------------------------------------------------------------------
// set lower boundary conditions for the neutrals
//----------------------------------------------------------------------

bool Neutrals::set_lower_bcs(Grid &grid,
                             Times time,
                             Indices indices) {

  std::string function = "Neutrals::set_lower_bcs";
  static int iFunction = -1;
  report.enter(function, iFunction);

  bool didWork = false;

  json bcs = input.get_boundary_condition_types();
  int64_t nGCs = grid.get_nGCs();
  int64_t iSpecies, iAlt, iDir;
  int64_t nLats = grid.get_nLats();
  int64_t nLons = grid.get_nLons();

  std::string bcsType = mklower(bcs["type"]);

  //-----------------------------------------------
  // MSIS BCs - only works if FORTRAN is enabled!
  //-----------------------------------------------

  // ALB changes to lower BCs only really work now for dipole grid. Don't use msis
  // if we are handed a dipole grid.

  if (bcsType == "msis" && !grid.IsDipole) {

    report.print(2, "Using MSIS for Boundary Conditions");

    Msis msis;

    if (!msis.is_ok()) {
      didWork = false;
      report.error("MSIS initialization not ok");

      if (report.test_verbose(0)) {
        std::cout << "MSIS Boundary Conditions asked for, ";
        std::cout << "but MSIS is not compiled! Yikes!\n";
      }
    } else
      didWork = true;

    msis.set_time(time);
    precision_t f107 = indices.get_f107(time.get_current());
    precision_t f107a = indices.get_f107a(time.get_current());
    msis.set_f107(f107, f107a);
    msis.set_ap(10.0);
    msis.set_locations(grid.geoLon_scgc.slice(0),
                       grid.geoLat_scgc.slice(0),
                       grid.geoAlt_scgc.slice(0));

    // This is just to check if MSIS is actually working:
    if (msis.is_valid_species("Tn"))
      // if it is, fill will temperature:
      temperature_scgc.slice(0) = msis.get_mat("Tn");
    else
      // if it is not, then fill with a value:
      temperature_scgc.slice(0).fill(initial_temperatures[0]);

    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      if (report.test_verbose(3))
        std::cout << "Setting Species : " << species[iSpecies].cName << "\n";

      if (msis.is_valid_species(species[iSpecies].cName)) {
        if (report.test_verbose(3))
          std::cout << "  Found in MSIS!\n";

        species[iSpecies].density_scgc.slice(0) =
          msis.get_mat(species[iSpecies].cName);
      } else {
        if (report.test_verbose(3))
          std::cout << "  NOT Found in MSIS - setting constant\n";

        species[iSpecies].density_scgc.slice(0).
        fill(species[iSpecies].lower_bc_density);
      }

    }

  } // type == Msis

  precision_t sh_ave;

  //-----------------------------------------------
  // Fill the lower+ ghost cells
  //-----------------------------------------------
  // - Planet BCs are in here too, can be refactored out
  // - Dipole grid must use planet BCs, for now.
  // - This kind-of assumes nGCs=2, so may need to be updated.
  // - If the first_lower_gc is at iAlt = 1, this may cause issues.
  // - The equator-most (j-hat) grid cell will be entirely below min_alt!
  for (int iLon = 0; iLon < nLons; iLon++) {
    for (int iLat = 0; iLat < nLats; iLat++) {

      // k-index of 1st lower ghost cell is not constant on the dipole grid.
      // On the latlon grid with nGCS=2, this will be 1
      iAlt = grid.first_lower_gc(iLon, iLat);
      temperature_scgc(iLon, iLat, iAlt) = initial_temperatures[0];
      // Set all lower ghost cells to bottom temperature:
      temperature_scgc.subcube(iLon, iLat, 0, iLon, iLat, iAlt - 1).fill(
        temperature_scgc(iLon, iLat, iAlt));

      precision_t t = temperature_scgc(iLon, iLat, 0);
      precision_t g = abs(grid.gravity_vcgc[2](iLon, iLat, iAlt));

      precision_t alt1 = grid.geoAlt_scgc(iLon, iLat, iAlt);
      precision_t alt0 = grid.altitude_lower_bc;
      precision_t dz = alt1 - alt0;

      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {

        precision_t m = mean_major_mass_scgc(iLon, iLat, iAlt);

        //if (m == 0)
        m = species[iSpecies].mass;

        precision_t h = cKB * t / (m * g);
        precision_t factor = exp(-dz / h);

        //-----------------------------------------------
        // Planet BCs - set to fixed constant values.
        //-----------------------------------------------
        if (bcsType == "planet" || grid.IsDipole) {

          // Fill all lower ghost cells density with lower boundary condition:
          species[iSpecies].density_scgc.subcube(iLon, iLat, 0,
                                                 iLon, iLat, iAlt - 1).fill(
                                                   factor *
                                                   species[iSpecies].lower_bc_density);
        }  // planet bc type

        // 1st ghost cell density is filled with a hydrostatic solution.
        sh_ave = (species[iSpecies].scale_height_scgc(iLon, iLat, iAlt)
                  + species[iSpecies].scale_height_scgc(iLon, iLat, iAlt + 1)) / 2;

        species[iSpecies].density_scgc(iLon, iLat, iAlt) =
          temperature_scgc(iLon, iLat, iAlt + 1)
          / temperature_scgc(iLon, iLat, iAlt)
          * species[iSpecies].density_scgc(iLon, iLat, iAlt + 1)
          * exp(-grid.dr_edge(iLon, iLat, iAlt) / sh_ave);

        // Vertical velocities: (In GITM this projected down with mesh coeffs)
        // Take lowest physical cell's vertical velocity and project it down nGCs cells.
        // All "GCs" lower than that have 0 vertical velocity since they're nonphysical.
        species[iSpecies].velocity_vcgc[2].subcube(
          iLon, iLat, iAlt - 1, size(1, 1, nGCs)).fill(
            species[iSpecies].velocity_vcgc[2](iLon, iLat, iAlt + 1));
        //project_onesided_alt_3rd(species[iSpecies].velocity_vcgc[2], grid, iAlt);

        if (iAlt > nGCs - 1) { // Fill all lower GCs w/ zero vertical velocity
          species[iSpecies].velocity_vcgc[2].subcube(iLon, iLat, 0,
                                                     size(1, 1, iAlt - 1)).zeros();

        }
      }
    }
  }

  // This section is for the dipole grid.  If the field-lines are
  // closed, then we will treat the N/S ghostcells as LOWER boundaries.
  // If thr grid is in the south, then treat the north bounday as the
  // lower boundary.  If the grid is in the north, treat the south boundary
  // as the lower boundary.
  // Because we are expecting to be chemically dominant, the lower BCs don't
  // matter as much for the ions.  We really just want to fill them with some
  // reasonable values.

  int64_t nX = grid.get_nX();
  int64_t nY = grid.get_nY();
  int64_t iX, iY, iYs, iYe, iFirst;

  if (grid.setNorthAsDown) {
    // First physical cell:
    iFirst = nY - nGCs - 1;
    iYs = nY - nGCs;
    iYe = nY;
  }

  if (grid.setSouthAsDown) {
    // First physical cell:
    iFirst = nGCs;
    iYs = 0;
    iYe = nGCs;
  }

  if (grid.setNorthAsDown || grid.setSouthAsDown) {

    for (iX = 0; iX < nX; iX++) {
      for (int iY = iYs; iY < iYe; iY++) {
        // Bulk Quantities:
        temperature_scgc.tube(iX, iY) = temperature_scgc.tube(iX, iFirst);

        // For each species:
        for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
          // Assume each species falls off a bit.
          // this BC shouldn't matter, since we are not going to do
          // horizontal advection on neutrals:
          species[iSpecies].density_scgc.tube(iX, iY) =
            0.95 * species[iSpecies].density_scgc.tube(iX, iFirst);
        }

        for (iAlt = 0; iAlt <= grid.first_lower_gc(iX, iY); iAlt++) {
          temperature_scgc(iX, iY, iAlt) =
            temperature_scgc(iX, iFirst, grid.first_lower_gc(iX, iFirst));
        }

      }
    }
  }

  didWork = true;

  calc_bulk_velocity();

  if (!didWork) {
    report.error("issue with lower BCs!");
    report.error("maybe check boundaryconditions type : " + bcsType);
  }

  report.exit(function);
  return didWork;
}

//----------------------------------------------------------------------
// set_horizontal_bcs
//   iDir tells which direction to set:
//      iDir = 0 -> +x
//      iDir = 1 -> +y
//      iDir = 2 -> -x
//      iDir = 3 -> -y
//----------------------------------------------------------------------

bool Neutrals::set_horizontal_bcs(int64_t iDir, Grid &grid) {

  std::string function = "Neutrals::set_horizontal_bcs";
  static int iFunction = -1;
  report.enter(function, iFunction);

  bool didWork = true;

  int64_t nX = grid.get_nX(), iX;
  int64_t nY = grid.get_nY(), iY;
  int64_t nAlts = grid.get_nAlts(true), iAlt;
  int64_t nGCs = grid.get_nGCs();
  int64_t iV;

  // iDir = 0 is right BC:
  if (iDir == 0) {
    for (iX = nX - nGCs; iX < nX; iX++) {
      for (iY = 0; iY < nY; iY++) {
        // Constant Gradient for Temperature:
        temperature_scgc.tube(iX, iY) =
          2 * temperature_scgc.tube(iX - 1, iY) -
          temperature_scgc.tube(iX - 2, iY);

        // Constant Value for Velocity:
        for (iV = 0; iV < 3; iV++)
          velocity_vcgc[iV].tube(iX, iY) = velocity_vcgc[iV].tube(iX - 1, iY);

        // Constant Gradient for densities:
        for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++)
          species[iSpecies].density_scgc.tube(iX, iY) =
            2 * species[iSpecies].density_scgc.tube(iX - 1, iY) -
            species[iSpecies].density_scgc.tube(iX - 2, iY);
      }
    }
  }

  // iDir = 2 is left BC:
  if (iDir == 2) {
    for (iX = nGCs - 1; iX >= 0; iX--) {
      for (iY = 0; iY < nY; iY++) {
        // Constant Gradient for Temperature:
        temperature_scgc.tube(iX, iY) =
          2 * temperature_scgc.tube(iX + 1, iY) -
          temperature_scgc.tube(iX + 2, iY);

        // Constant Value for Velocity:
        for (iV = 0; iV < 3; iV++)
          velocity_vcgc[iV].tube(iX, iY) = velocity_vcgc[iV].tube(iX + 1, iY);

        // Constant Gradient for densities:
        for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++)
          species[iSpecies].density_scgc.tube(iX, iY) =
            2 * species[iSpecies].density_scgc.tube(iX + 1, iY) -
            species[iSpecies].density_scgc.tube(iX + 2, iY);
      }
    }
  }

  // iDir = 1 is upper BC:
  if (iDir == 1) {
    for (iX = 0; iX < nX; iX++) {
      for (iY = nX - nGCs; iY < nY; iY++) {
        // Constant Gradient for Temperature:
        temperature_scgc.tube(iX, iY) =
          2 * temperature_scgc.tube(iX, iY - 1) -
          temperature_scgc.tube(iX, iY - 2);

        // Constant Value for Velocity:
        for (iV = 0; iV < 3; iV++)
          velocity_vcgc[iV].tube(iX, iY) = velocity_vcgc[iV].tube(iX, iY - 1);

        // Constant Gradient for densities:
        for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++)
          species[iSpecies].density_scgc.tube(iX, iY) =
            2 * species[iSpecies].density_scgc.tube(iX, iY - 1) -
            species[iSpecies].density_scgc.tube(iX, iY - 2);
      }
    }
  }

  // iDir = 3 is lower BC:
  if (iDir == 3) {
    for (iX = 0; iX < nX; iX++) {
      for (iY = nGCs - 1; iY >= 0; iY--) {
        // Constant Gradient for Temperature:
        temperature_scgc.tube(iX, iY) =
          2 * temperature_scgc.tube(iX, iY + 1) -
          temperature_scgc.tube(iX, iY + 2);

        // Constant Value for Velocity:
        for (iV = 0; iV < 3; iV++)
          velocity_vcgc[iV].tube(iX, iY) = velocity_vcgc[iV].tube(iX, iY + 1);

        // Constant Gradient for densities:
        for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++)
          species[iSpecies].density_scgc.tube(iX, iY) =
            2 * species[iSpecies].density_scgc.tube(iX, iY + 1) -
            species[iSpecies].density_scgc.tube(iX, iY + 2);
      }
    }
  }

  report.exit(function);
  return didWork;
}
