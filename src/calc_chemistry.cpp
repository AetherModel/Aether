// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>

#include "aether.h"

// -----------------------------------------------------------------------------
// calculates the chemical reactions in the model by taking into account
// EUV ionization and chemistry as this time
// -----------------------------------------------------------------------------

void Chemistry::calc_chemistry(Neutrals &neutrals,
                               Ions &ions,
                               Times time,
                               Grid grid,
                               Report &report) {

  int iSpecies;

  std::string function = "Chemistry::calc_chemistry";
  static int iFunction = -1;
  report.enter(function, iFunction);

  precision_t dt = time.get_dt();

  // ------------------------------------
  // Calculate electron densities
  // ------------------------------------

  ions.fill_electrons(report);

  // ----------------------------------------------------------
  // Initialize the sources and losses with EUV stuff:
  // ----------------------------------------------------------

  // Neutrals have losses due to ionization
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    neutrals.species[iSpecies].losses_scgc =
      neutrals.species[iSpecies].ionization_scgc;
    neutrals.species[iSpecies].sources_scgc.zeros();
  }

  // Ions have sources due to ionization
  for (iSpecies = 0; iSpecies < nIons; iSpecies++) {
    ions.species[iSpecies].losses_scgc.zeros();
    ions.species[iSpecies].sources_scgc =
      ions.species[iSpecies].ionization_scgc;
  }

  // ----------------------------------------------------
  // Calculate the chemical sources and losses
  // ----------------------------------------------------

  calc_chemical_sources(neutrals, ions, report);

  // ---------------------------------------------------------
  // Once sources and losses are done, solve for new densities
  // ---------------------------------------------------------

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    neutrals.species[iSpecies].density_scgc =
      solver_chemistry(neutrals.species[iSpecies].density_scgc,
                       neutrals.species[iSpecies].sources_scgc,
                       neutrals.species[iSpecies].losses_scgc,
                       dt);
  }

  for (iSpecies = 0; iSpecies < nIons; iSpecies++) {
    ions.species[iSpecies].density_scgc =
      solver_chemistry(ions.species[iSpecies].density_scgc,
                       ions.species[iSpecies].sources_scgc,
                       ions.species[iSpecies].losses_scgc,
                       dt);
  }

  // ---------------------------------------------------------
  // Recalculate electrons
  // ---------------------------------------------------------

  ions.fill_electrons(report);

  report.exit(function);
  return;
}
