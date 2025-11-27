// Copyright 2023, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md
//
// initial version - A. Ridley - May 27, 2023

#include "aether.h"

// -----------------------------------------------------------------------------
// Set boundary conditions for the ions.
// The dipolar grid is fundamentally different than the sphere/cubesphere grids.
// We need to treat them differently.
// -----------------------------------------------------------------------------

//----------------------------------------------------------------------
// set_bcs - This is for setting the vertical BCs
//----------------------------------------------------------------------

bool Ions::set_bcs(Grid &grid,
                   Times time,
                   Indices indices) {

  std::string function = "Ions::set_bcs";
  static int iFunction = -1;
  report.enter(function, iFunction);

  bool didWork = true;

  if (grid.get_nZ(false) > 1) {
    didWork = set_lower_bcs(grid, time, indices);

    if (didWork)
      didWork = set_upper_bcs(grid);

    if (didWork)
      fill_electrons();
  }

  if (!didWork)
    report.error("issue with ion BCs!");

  report.exit(function);
  return didWork;
}

//----------------------------------------------------------------------
// set upper boundary conditions for the ions
//----------------------------------------------------------------------

bool Ions::set_upper_bcs(Grid &grid) {

  std::string function = "Ions::set_upper_bcs";
  static int iFunction = -1;
  report.enter(function, iFunction);

  bool didWork = true;

  int64_t nAlts = grid.get_nZ();
  int64_t nX = grid.get_nX(), iX;
  int64_t nY = grid.get_nY(), iY;
  int64_t nGCs = grid.get_nGCs();
  int64_t iAlt;
  arma_mat h;
  arma_mat aveT;

  // If we are on the dipole grid and our field-lines are closed, then
  // we don't want to set upper boundary conditions, since we will
  // be message passing them.

  if (!grid.IsClosed) {

    for (iAlt = nAlts - nGCs; iAlt < nAlts; iAlt++) {
      // Bulk Quantities:
      // Constant gradient (ignoring grid spacing...)
      temperature_scgc.slice(iAlt) =
        2 * temperature_scgc.slice(iAlt - 1) - temperature_scgc.slice(iAlt - 2);

      // For each species:
      for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
        // Constant gradient (ignoring grid spacing...)
        species[iSpecies].temperature_scgc.slice(iAlt) =
          2 * species[iSpecies].temperature_scgc.slice(iAlt - 1) -
          species[iSpecies].temperature_scgc.slice(iAlt - 2);

        aveT = (species[iSpecies].temperature_scgc.slice(iAlt) +
                electron_temperature_scgc.slice(iAlt));
        // Calculate scale height for the species:
        h = cKB / species[iSpecies].mass *
            species[iSpecies].temperature_scgc.slice(iAlt) /
            abs(grid.gravity_vcgc[2].slice(iAlt));
        // Assume each species falls of with (modified) hydrostatic:
        species[iSpecies].density_scgc.slice(iAlt) =
          species[iSpecies].temperature_scgc.slice(iAlt) /
          species[iSpecies].temperature_scgc.slice(iAlt - 1) %
          species[iSpecies].density_scgc.slice(iAlt - 1) %
          exp(-grid.dk_edge_m.slice(iAlt) / h);
        species[iSpecies].velocity_vcgc[2].slice(iAlt).zeros();
      }
    }

  }

  report.exit(function);
  return didWork;
}

//----------------------------------------------------------------------
// set lower boundary conditions for the ions
//----------------------------------------------------------------------

bool Ions::set_lower_bcs(Grid &grid, Times time, Indices indices) {

  std::string function = "Ions::set_lower_bcs";
  static int iFunction = -1;
  report.enter(function, iFunction);

  bool didWork = true;

  int64_t nAlts = grid.get_nZ();
  int64_t nX = grid.get_nX(), iX;
  int64_t nY = grid.get_nY(), iY, iYs, iYe;
  int64_t nGCs = grid.get_nGCs();
  int64_t iAlt, iFirst;
  arma_mat h;
  arma_mat aveT;

  // This is true for all grids:
  for (iX = 0; iX < nX; iX++) {
    for (int iY = 0; iY < nY; iY++) {
      iFirst = grid.first_lower_gc(iX, iY);

      for (iAlt = iFirst; iAlt >= 0; iAlt--) {
        // Bulk Quantities:
        temperature_scgc.slice(iAlt) = temperature_scgc.slice(iFirst + 1);

        // For each species:
        for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
          // assign all species temperatures the bulk temperature:
          species[iSpecies].temperature_scgc.slice(iAlt) =
            temperature_scgc.slice(iAlt);
          // Assume each species falls off a bit.
          // this BC shouldn't matter, since the bottom of the code
          // should be in chemical equalibrium:
          species[iSpecies].density_scgc.slice(iAlt) =
            0.95 * species[iSpecies].density_scgc.slice(iFirst + 1);
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

  if (grid.setNorthAsDown) {
    // First physical cell:
    iFirst = nY - nGCs - 2;
    iYs = nY - nGCs - 1;
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
          // assign all species temperatures the bulk temperature:
          species[iSpecies].temperature_scgc.tube(iX, iY) =
            temperature_scgc.tube(iX, iFirst);
          // Assume each species falls off a bit.
          // this BC shouldn't matter, since the bottom of the code
          // should be in chemical equalibrium:
          species[iSpecies].density_scgc.tube(iX, iY) =
            0.95 * species[iSpecies].density_scgc.tube(iX, iFirst);
        }

        for (iAlt = 0; iAlt <= grid.first_lower_gc(iX, iY); iAlt++) {
          //std::cout << "ion bcs, dipole, setnorth iAlt : "
          //          << iAlt << " "
          //          << grid.first_lower_gc(iX, iFirst) << " "
          //          << temperature_scgc(iX, iFirst, grid.first_lower_gc(iX, iFirst) + 1) << " "
          //          << temperature_scgc(iX, iFirst, iAlt) << " "
          //          << grid.geoAlt_scgc(iX, iFirst, grid.first_lower_gc(iX,
          //                              iFirst) + 1) / 1000.0 << "\n";
          temperature_scgc(iX, iY, iAlt) =
            temperature_scgc(iX, iFirst, grid.first_lower_gc(iX, iFirst) + 1);
        }

      }
    }
  }

  report.exit(function);
  return didWork;
}

