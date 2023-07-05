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

bool Neutrals::set_bcs(Grid grid,
                       Times time,
                       Indices indices,
                       Inputs input,
                       Report &report) {

  std::string function = "Neutrals::set_bcs";
  static int iFunction = -1;
  report.enter(function, iFunction);

  bool didWork;

  didWork = set_lower_bcs(grid, time, indices, input, report);
  didWork = set_upper_bcs(grid, input, report);
  calc_mass_density(report);
  
  report.exit(function);
  return didWork;
}

//----------------------------------------------------------------------
// set lower boundary conditions for the neutrals
//----------------------------------------------------------------------

bool Neutrals::set_upper_bcs(Grid grid,
                             Inputs input,
                             Report &report) {

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
  
  for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    for (iAlt = nAlts - nGCs; iAlt < nAlts; iAlt++) {

      // Allow upflow, but not downflow:
      for (iX = nGCs; iX < nX - nGCs; iX++)
	for (iY = nGCs; iY < nY - nGCs; iY++)
	  if (species[iSpecies].velocity_vcgc[2](iX, iY, iAlt-1) > 0)
	    species[iSpecies].velocity_vcgc[2](iX, iY, iAlt) =
	      species[iSpecies].velocity_vcgc[2](iX, iY, iAlt-1);
	  else
	    species[iSpecies].velocity_vcgc[2](iX, iY, iAlt) = 0.0;
      
      h = cKB * temperature_scgc.slice(iAlt) / (species[iSpecies].mass * grid.gravity_scgc.slice(iAlt));
      species[iSpecies].density_scgc.slice(iAlt) = 
	species[iSpecies].density_scgc.slice(iAlt-1) %
	exp(-grid.dalt_lower_scgc.slice(iAlt) / h);
    }
  }

  temperature_scgc.slice(nAlts - 2) = temperature_scgc.slice(nAlts - 3);
  temperature_scgc.slice(nAlts - 1) = temperature_scgc.slice(nAlts - 2);
  
  report.exit(function);
  return didWork;
}

//----------------------------------------------------------------------
// set lower boundary conditions for the neutrals
//----------------------------------------------------------------------

bool Neutrals::set_lower_bcs(Grid grid,
                             Times time,
                             Indices indices,
                             Inputs input,
                             Report &report) {

  std::string function = "Neutrals::set_lower_bcs";
  static int iFunction = -1;
  report.enter(function, iFunction);

  bool didWork = true;

  json bcs = input.get_boundary_condition_types();

  //-----------------------------------------------
  // MSIS BCs - only works if FORTRAN is enabled!
  //-----------------------------------------------

  if (bcs["type"] == "Msis") {

    report.print(2, "Using MSIS for Boundary Conditions");

    Msis msis(input);

    if (!msis.is_ok()) {
      didWork = false;

      if (report.test_verbose(0)) {
        std::cout << "MSIS Boundary Conditions asked for, ";
        std::cout << "but MSIS is not compiled! Yikes!\n";
      }
    }

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

    for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
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

  //-----------------------------------------------
  // Planet BCs - set to fixed constant values.
  //-----------------------------------------------

  if (bcs["type"] == "Planet") {

    report.print(2, "setting lower bcs to planet");
    arma_mat h;
    
    // Set the lower boundary condition:
    for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      species[iSpecies].density_scgc.slice(0).
	fill(species[iSpecies].lower_bc_density);
      // fill the second grid cell with a hydrostatic solution:
      h = cKB * temperature_scgc.slice(0) / (species[iSpecies].mass * grid.gravity_scgc.slice(0));
      species[iSpecies].density_scgc.slice(1) = 
	species[iSpecies].density_scgc.slice(0) %
	exp(-grid.dalt_lower_scgc.slice(1) / h);
      species[iSpecies].velocity_vcgc[2].slice(0).zeros();
      species[iSpecies].velocity_vcgc[2].slice(1).zeros();
    }
    temperature_scgc.slice(0).fill(initial_temperatures[0]);
    temperature_scgc.slice(1).fill(initial_temperatures[0]);
  } // type == Planet

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

bool Neutrals::set_horizontal_bcs(int64_t iDir, Grid grid, Report &report) {

  std::string function = "Neutrals::set_horizontal_bcs";
  static int iFunction = -1;
  report.enter(function, iFunction);

  bool didWork = true;

  int64_t nX = grid.get_nX(), iX;
  int64_t nY = grid.get_nY(), iY;
  int64_t nAlts = grid.get_nAlts(), iAlt;
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

  // iDir = 2 is left BC:
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
