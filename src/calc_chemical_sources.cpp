// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>

#include "../include/report.h"
#include "../include/chemistry.h"

void Chemistry::calc_chemical_sources(Neutrals &neutrals,
                                      Ions &ions,
                                      Report &report) {

  std::string function = "Chemistry::calc_chemical_sources";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t iReaction, iLoss, iSource;
  float change, rate;
  int IsNeutral, id_;

  // This is make change the same size as the grid:
  fcube change3d = neutrals.temperature_scgc;
  // Then zero is out:
  change3d.zeros();

  for (iReaction = 0; iReaction < nReactions; iReaction++) {

    if (report.test_verbose(8)) display_reaction(reactions[iReaction]);
    
    // First calculate reaction rate:

    rate = reactions[iReaction].rate;

    // Second calculate the amount of change:

    change3d.fill(rate);
    for (iLoss = 0; iLoss < reactions[iReaction].nLosses; iLoss++) {
      IsNeutral = reactions[iReaction].losses_IsNeutral[iLoss];
      id_ = reactions[iReaction].losses_ids[iLoss];
      if (IsNeutral)
        change3d = change3d % neutrals.neutrals[id_].density_scgc;
      else
        change3d = change3d % ions.species[id_].density_scgc;
    }

    // Third add change to the different consituents:

    //  (a) First to losses:
    for (iLoss = 0; iLoss < reactions[iReaction].nLosses; iLoss++) {
      IsNeutral = reactions[iReaction].losses_IsNeutral[iLoss];
      id_ = reactions[iReaction].losses_ids[iLoss];
      if (IsNeutral)
        neutrals.neutrals[id_].losses_scgc =
          neutrals.neutrals[id_].losses_scgc + change3d;
      else
        ions.species[id_].losses_scgc =
          ions.species[id_].losses_scgc + change3d;
    }

    //  (b) Second to sources:
    for (iSource = 0; iSource < reactions[iReaction].nSources; iSource++) {
      IsNeutral = reactions[iReaction].sources_IsNeutral[iSource];
      id_ = reactions[iReaction].sources_ids[iSource];
      if (IsNeutral)
        neutrals.neutrals[id_].sources_scgc =
          neutrals.neutrals[id_].sources_scgc + change3d;
      else
        ions.species[id_].sources_scgc =
          ions.species[id_].sources_scgc + change3d;
    }  // for iSource
  }  // for iReaction
  report.exit(function);
}
