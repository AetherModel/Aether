// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>

#include "aether.h"

// -----------------------------------------------------------------------------
// loops through all of the chemical reactions doing 3 things:
//   1. Determine change (per unit time) of particles of loss
//   2. Add this to sources (keeping track of ions v. neutrals)
//   3. Add this to losses (keeping track of ions v. neutrals)
// -----------------------------------------------------------------------------

void Chemistry::calc_chemical_sources(Neutrals &neutrals,
                                      Ions &ions,
                                      Report &report) {

  std::string function = "Chemistry::calc_chemical_sources";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t iReaction, iLoss, iSource;
  precision_t rate;
  bool IsNeutral;
  int id_;

  // This is make change the same size as the grid:
  arma_cube change3d = neutrals.temperature_scgc;
  // Then zero is out:
  change3d.zeros();

  for (iReaction = 0; iReaction < nReactions; iReaction++) {

    if (report.test_verbose(3)) {
      std::cout << "===> Reaction : " << iReaction
                << " of " << nReactions << "\n";
      display_reaction(reactions[iReaction]);
    }

    // Zero calculate reaction rate:

    rate = reactions[iReaction].rate;

    // First calculate the amount of change:
    //    Change is calculated as
    //    reaction rate * loss den 1 * loss den 2 (* loss den 3 if needed)

    change3d.fill(rate);

    for (iLoss = 0; iLoss < reactions[iReaction].nLosses; iLoss++) {
      IsNeutral = reactions[iReaction].losses_IsNeutral[iLoss];
      id_ = reactions[iReaction].losses_ids[iLoss];

      if (IsNeutral)
        change3d = change3d % neutrals.species[id_].density_scgc;
      else
        change3d = change3d % ions.species[id_].density_scgc;
    }

    // Second add change to the different consituents:
    for (iLoss = 0; iLoss < reactions[iReaction].nLosses; iLoss++) {
      IsNeutral = reactions[iReaction].losses_IsNeutral[iLoss];
      id_ = reactions[iReaction].losses_ids[iLoss];

      if (IsNeutral)
        neutrals.species[id_].losses_scgc =
          neutrals.species[id_].losses_scgc + change3d;
      else
        ions.species[id_].losses_scgc =
          ions.species[id_].losses_scgc + change3d;
    }

    // Third add change to the difference constituents:
    for (iSource = 0; iSource < reactions[iReaction].nSources; iSource++) {
      IsNeutral = reactions[iReaction].sources_IsNeutral[iSource];
      id_ = reactions[iReaction].sources_ids[iSource];

      if (IsNeutral)
        neutrals.species[id_].sources_scgc =
          neutrals.species[id_].sources_scgc + change3d;
      else
        ions.species[id_].sources_scgc =
          ions.species[id_].sources_scgc + change3d;
    }  // for iSource
  }  // for iReaction

  report.exit(function);
    //cout << "Found it: calc_chemical_sources in calc_chemical_sources.cpp\n";
}
