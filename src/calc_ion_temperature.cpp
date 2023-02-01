// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <math.h>
#include <iostream>
#include <fstream>

#include "aether.h"


// --------------------------------------------------------------------------
// Initialize the ion temperature - set equal to the neutral temperature
// --------------------------------------------------------------------------

void Ions::init_ion_temperature(Neutrals neutrals, Grid grid, Report &report) {

  int64_t iIon;

  for (iIon = 0; iIon < nIons; iIon++)
    species[iIon].temperature_scgc = neutrals.temperature_scgc;

  temperature_scgc = neutrals.temperature_scgc;

  return;
}


// --------------------------------------------------------------------------
// Calculate the ion temperature
// --------------------------------------------------------------------------

void Ions::calc_ion_temperature(Neutrals neutrals, Grid grid,
                                Times time, Inputs input, Report &report) {

  std::string function = "Ions::calc_ion_temperature";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t iIon, iLon, iLat, nSpecs;
  int64_t nLons = grid.get_nLons();
  int64_t nLats = grid.get_nLats();
  int64_t nAlts = grid.get_nAlts();

  arma_vec temp1d(nAlts);
  arma_vec lambda1d(nAlts);
  arma_vec front1d(nAlts);
  arma_vec dalt1d(nAlts);
  arma_vec conduction1d(nAlts);

  arma_cube tempT(nLons, nLats, nAlts);
  arma_cube tempD(nLons, nLats, nAlts);

  // Get the time step size
  precision_t dt = time.get_dt();

  // Loop over all species or assume only bulk calculation
  if (input.get_do_calc_bulk_ion_temp())
    // First ion species only, currently is O+
    nSpecs = 1;
  else
    nSpecs = nIons;

  if (report.test_verbose(4)) {
    std::cout << "Bulk ion temp flag: " << input.get_do_calc_bulk_ion_temp()
              << " so 'number of ions' is " << nSpecs << "\n";
  }

  // Loop over all species or assume only bulk calculation
  for (iIon = 0; iIon < nSpecs; iIon++) {

    for (iLon = 0; iLon < nLons; iLon++) {
      for (iLat = 0; iLat < nLats; iLat++) {

        // ---------------------------------------------------------------------
        // Calculate heat flux (conduction) in 1D; loop over all lat,lon
        // ---------------------------------------------------------------------
        temp1d   = species[iIon].temperature_scgc.tube(iLon, iLat);
        lambda1d = 25.0 * pow(cKB, 2) * pow(temp1d, 2.5) / species[iIon].mass
                   / species[iIon].nu_ion_ion[iIon] / 8.0;
        front1d  = 2.0 / species[iIon].density_scgc.tube(iLon, iLat)
                   / cKB / 3.0;
        dalt1d   = grid.dalt_lower_scgc.tube(iLon, iLat);

        conduction1d.zeros();    // reset temp variable to zero

        conduction1d = solver_conduction(temp1d, lambda1d, front1d, dt, dalt1d);

        // The conduction solver gives Tnew-Told, so divide by dt
        conduction_scgc.tube(iLon, iLat) = conduction1d / dt;

      } // Lats
    } // Lons

    // -------------------------------------------------------------------------
    // Add temperature terms together to advance ion temperature
    // As more temperature terms get coded, they are added to the parenthesis
    // for inclusion in the advancement of the ion temperature
    // -------------------------------------------------------------------------
    if (!input.get_do_calc_bulk_ion_temp()) {
      species[iIon].temperature_scgc = species[iIon].temperature_scgc +
                                       dt * (conduction_scgc);
    }
  } // Ions

  if (!input.get_do_calc_bulk_ion_temp()) {
    // Use the density averaged temperature to fill the bulk temperature
    tempT.zeros();
    tempD.zeros();

    for (iIon = 0; iIon < nIons; iIon++) {
      tempT = tempT +
              species[iIon].temperature_scgc % species[iIon].density_scgc;
      tempD = tempD +
              species[iIon].density_scgc;
    }

    temperature_scgc = tempT / tempD;
  }

  if (input.get_do_calc_bulk_ion_temp()) {
    // Add temperature terms together to advance bulk ion temperature
    temperature_scgc = temperature_scgc + dt * (conduction_scgc);

    // Use the bulk ion temperature to fill all ion specie temperatures
    for (iIon = 0; iIon < nIons; iIon++)
      species[iIon].temperature_scgc = temperature_scgc;
  }

  report.exit(function);
  return;
}
