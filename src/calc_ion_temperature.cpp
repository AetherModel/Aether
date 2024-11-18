// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <math.h>
#include <iostream>
#include <fstream>

#include "aether.h"

// --------------------------------------------------------------------------
// Initialize the ion temperature - set equal to the neutral temperature
// --------------------------------------------------------------------------

void Ions::init_ion_temperature(Neutrals neutrals, Grid grid) {

  int64_t iIon;

  for (iIon = 0; iIon < nSpecies; iIon++) {
    species[iIon].temperature_scgc = neutrals.temperature_scgc;

    // This is the first place where we have the ions, neutrals, and grid:
    species[iIon].nu_ion_neutral_vcgc = make_cube_vector(grid.get_nLons(),
                                                         grid.get_nLats(),
                                                         grid.get_nAlts(),
                                                         neutrals.nSpecies);
  }

  temperature_scgc = neutrals.temperature_scgc;

  return;
}

// --------------------------------------------------------------------------
// Calculate the ion temperature
// --------------------------------------------------------------------------

void Ions::calc_ion_temperature(Neutrals neutrals, Grid grid,
                                Times time) {

  std::string function = "Ions::calc_ion_temperature";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t iIon, iLon, iLat, nSpecs, jIon;
  int64_t nLons = grid.get_nLons();
  int64_t nLats = grid.get_nLats();
  int64_t nAlts = grid.get_nAlts();
  int64_t nGCs = grid.get_nGCs();
  precision_t Mi, Mj;

  arma_vec temp1d(nAlts);
  arma_vec lambda1d(nAlts);
  arma_vec front1d(nAlts);
  arma_vec dalt1d(nAlts);
  arma_vec conduction1d(nAlts);
  arma_vec sources1d(nAlts);
  arma_vec ratios(nAlts);
  arma_vec density_ratio(nAlts);

  arma_cube tempT(nLons, nLats, nAlts);
  arma_cube tempD(nLons, nLats, nAlts);

  // Get the time step size
  precision_t dt = time.get_dt();

  //temperature_scgc = 200.0 + sqrt(grid.geoAlt_scgc / 1000.0 - 90.0) * 100.0;

  //for (iIon = 0; iIon < nSpecies; iIon++)
  //  species[iIon].temperature_scgc = temperature_scgc;

  //species[iIon].temperature_scgc = neutrals.temperature_scgc;

  //report.exit(function);
  //return;

  // Loop over all species or assume only bulk calculation
//  if (input.get_do_calc_bulk_ion_temp())
//    // First ion species only, currently is O+
//    nSpecs = 1;
//  else
    nSpecs = nSpecies;

  if (report.test_verbose(4)) {
    std::cout << "Bulk ion temp flag: " << input.get_do_calc_bulk_ion_temp() ?
              "true" : "false";
    std::cout << " so 'number of ions' is " << nSpecs << "\n";
  }

  // Loop over all species or assume only bulk calculation
  for (iIon = 0; iIon < nSpecs; iIon++) {
    //std::cout << "iIon : " << iIon << "\n";
    Mi = species[iIon].mass / cAMU;
    for (iLon = nGCs; iLon < nLons - nGCs; iLon++) {
      for (iLat = nGCs; iLat < nLats - nGCs; iLat++) {
        //std::cout << "iLon, iLat : " << iLon << " " << iLat << "\n";
        // ---------------------------------------------------------------------
        // Calculate heat flux (conduction) in 1D; loop over all lat,lon
        // ---------------------------------------------------------------------
        ratios.zeros();
        for (jIon = 0; jIon < nSpecies; jIon++) {
          if (jIon != iIon) {
            Mj = species[jIon].mass / cAMU;
            density_ratio = species[jIon].density_scgc.tube(iLon, iLat) / 
              species[iIon].density_scgc.tube(iLon, iLat);
            density_ratio.clamp(0.001, 1000.0);
            ratios = ratios + density_ratio * 
              (species[jIon].charge * species[jIon].charge / 
               species[iIon].charge / species[iIon].charge) *
              sqrt(Mj / (Mi + Mj)) *
              (3 * Mi * Mi + 1.6 * Mi * Mj + 1.3 * Mj * Mj)/
              ((Mi + Mj) * (Mi + Mj));
          }
        }
        temp1d = species[iIon].temperature_scgc.tube(iLon, iLat);
        temp1d(0) = neutrals.temperature_scgc(iLon, iLat, 0);
        temp1d(1) = neutrals.temperature_scgc(iLon, iLat, 1);
        lambda1d = 3.1e6 / sqrt(Mi) / pow(species[iIon].charge, 4) *
           pow(temp1d, 2.5) % (1 + 1.75 * ratios) * cE;
        //lambda1d = 25.0 * cKB * pow(temp1d, 2.5) * (cKB / species[iIon].mass)
        //           / species[iIon].nu_ion_ion[iIon] / 8.0;
        front1d  = 3.0 / 2.0 * cKB * species[iIon].density_scgc.tube(iLon, iLat);
        dalt1d   = grid.dalt_lower_scgc.tube(iLon, iLat);
        sources1d = (heating_neutral_friction_scgc.tube(iLon, iLat) +
                    heating_neutral_heat_transfer_scgc.tube(iLon, iLat));
        sources1d = sources1d / front1d;

        //std::cout << "iIon lab : " << iIon << "\n" << lambda1d << "\n source:\n " << sources1d << "\n";

        conduction1d.zeros();    // reset temp variable to zero

        conduction1d = solver_conduction(temp1d, 
                                         lambda1d, 
                                         front1d, 
                                         sources1d, 
                                         dalt1d,
                                         dt/100., 
                                         nGCs, 
                                         false);

        // The conduction solver gives Tnew-Told, so divide by dt
        conduction1d.clamp(200, 5000);
        species[iIon].temperature_scgc.tube(iLon, iLat) = conduction1d;

        //std::cout << "temp : " << conduction1d << "\n";

      } // Lats
    } // Lons
  } // Ions

  //if (!input.get_do_calc_bulk_ion_temp()) {
    // Use the density averaged temperature to fill the bulk temperature
    tempT.zeros();
    tempD.zeros();

    for (iIon = 0; iIon < nSpecies; iIon++) {
      tempT = tempT +
              species[iIon].temperature_scgc % species[iIon].density_scgc;
      tempD = tempD +
              species[iIon].density_scgc;
    }

    temperature_scgc = tempT / tempD;
  //}
/*
  if (input.get_do_calc_bulk_ion_temp()) {
    // Add temperature terms together to advance bulk ion temperature
    temperature_scgc = temperature_scgc + dt * (conduction_scgc);

    // Use the bulk ion temperature to fill all ion specie temperatures
    for (iIon = 0; iIon < nSpecies; iIon++)
      species[iIon].temperature_scgc = temperature_scgc;
  }
*/

  std::cout << "ion temp : " << temperature_scgc(2,2,20) << " " << neutrals.temperature_scgc(2,2,20) << "\n";
  report.exit(function);
  return;
}
