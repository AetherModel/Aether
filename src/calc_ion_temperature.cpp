// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <math.h>
#include <iostream>
#include <fstream>

#include "aether.h"

// --------------------------------------------------------------------------
// Initialize the ion temperature - set equal to the neutral temperature
// --------------------------------------------------------------------------

void Ions::init_ion_temperature(Neutrals neutrals, Grid &grid) {

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

  // For electron temperature, we need to check if some species are present or not.
  // Do this check now & warn if needed:
  if ((neutrals.get_species_id("O") == -1)
      || (neutrals.get_species_id("O2") == -1)
      || (neutrals.get_species_id("N2") == -1)) {
    if (input.get_do_photoelectron_heating()
        || input.get_do_ionization_heating()
        || input.get_do_electron_neutral_elastic_collisional_heating())
      report.error("Your electron temperature sources require neutral O, O2, and N2 to be present.");
  }

  return;
}

// --------------------------------------------------------------------------
// Calculate the ion temperature
// --------------------------------------------------------------------------

void Ions::calc_ion_temperature(const Neutrals &neutrals, Grid &grid,
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

  nSpecs = nSpecies;

  if (report.test_verbose(4)) {
    std::cout << "Bulk ion temp flag: " << input.get_do_calc_bulk_ion_temp() ?
              "true" : "false";
    std::cout << " so 'number of ions' is " << nSpecs << "\n";
  }

  calc_lambda();

  // -------------------------------------------------
  // This is for calculating the bulk temperature:
  // -------------------------------------------------
  if (input.get_do_calc_bulk_ion_temp()) {
    for (iLon = nGCs; iLon < nLons - nGCs; iLon++) {
      for (iLat = nGCs; iLat < nLats - nGCs; iLat++) {
        temp1d = temperature_scgc.tube(iLon, iLat);
        lambda1d = lambda.tube(iLon, iLat);
        lambda1d(1) = lambda1d(2);
        lambda1d(0) = lambda1d(2);
        front1d  = 3.0 / 2.0 * cKB * density_scgc.tube(iLon, iLat);
        dalt1d   = grid.dk_edge_m.tube(iLon, iLat);
        sources1d = (heating_neutral_friction_scgc.tube(iLon, iLat) +
                     heating_neutral_heat_transfer_scgc.tube(iLon, iLat));
        sources1d = sources1d / front1d;
        conduction1d.zeros();    // reset temp variable to zero
        conduction1d = solver_conduction(temp1d,
                                         lambda1d,
                                         front1d,
                                         sources1d,
                                         dalt1d,
                                         dt / 10.,
                                         nGCs,
                                         false);
        // The conduction solver gives Tnew-Told, so divide by dt
        conduction1d.clamp(200, 5000);
        temperature_scgc.tube(iLon, iLat) = conduction1d;
      }
    }

    for (iIon = 0; iIon < nSpecies; iIon++)
      species[iIon].temperature_scgc = temperature_scgc;

  } else {

    // -------------------------------------------------
    // This is for calculating the individual temperature:
    // -------------------------------------------------

    for (iIon = 0; iIon < nSpecies; iIon++) {
      for (iLon = nGCs; iLon < nLons - nGCs; iLon++) {
        for (iLat = nGCs; iLat < nLats - nGCs; iLat++) {
          temp1d = species[iIon].temperature_scgc.tube(iLon, iLat);
          temp1d(0) = neutrals.temperature_scgc(iLon, iLat, 0);
          temp1d(1) = neutrals.temperature_scgc(iLon, iLat, 1);
          lambda1d = species[iIon].lambda.tube(iLon, iLat);
          lambda1d(1) = lambda1d(2);
          lambda1d(0) = lambda1d(2);
          front1d  = 3.0 / 2.0 * cKB * species[iIon].density_scgc.tube(iLon, iLat);
          dalt1d   = grid.dk_edge_m.tube(iLon, iLat);
          sources1d = (species[iIon].heating_neutral_friction_scgc.tube(iLon, iLat) +
                       species[iIon].heating_neutral_heat_transfer_scgc.tube(iLon, iLat));
          sources1d = sources1d / front1d;

          conduction1d.zeros();    // reset temp variable to zero
          conduction1d = solver_conduction(temp1d,
                                           lambda1d,
                                           front1d,
                                           sources1d,
                                           dalt1d,
                                           dt / 10.,
                                           nGCs,
                                           false);

          conduction1d.clamp(200, 5000);
          species[iIon].temperature_scgc.tube(iLon, iLat) = conduction1d;
        } // Lats
      } // Lons
    } // Ions

    tempT.zeros();
    tempD.zeros();

    for (iIon = 0; iIon < nSpecies; iIon++) {
      tempT = tempT +
              species[iIon].temperature_scgc % species[iIon].density_scgc;
      tempD = tempD +
              species[iIon].density_scgc;
    }

    temperature_scgc = tempT / tempD;
  }

  report.exit(function);
  return;
}
