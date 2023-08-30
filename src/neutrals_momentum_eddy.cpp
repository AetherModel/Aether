// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "aether.h"

void Neutrals::vertical_momentum_eddy(Grid &gGrid) {

  std::string function = "Neutrals::vertical_momentum_eddy";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t iSpecies, iSpecies_;

  // Initialize all of the accelerations to zero:
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    species[iSpecies].acc_eddy.zeros();
  
  if (input.get_use_eddy_momentum()) {

    arma_cube log_cons;
    arma_cube grad_cons;

    // Only consider the advected species:
    for (iSpecies = 0; iSpecies < nSpeciesAdvect; iSpecies++) {
      iSpecies_ = species_to_advect[iSpecies];
      // Take the natural log of each concentration cube:
      log_cons = log(species[iSpecies_].concentration_scgc);

      // calculate gradient:
      grad_cons = calc_gradient_alt(log_cons, gGrid);
      species[iSpecies_].acc_eddy = - grad_cons % kappa_eddy_scgc;
    }
  }

  report.exit(function);
  return;
} // neutral_friction_momentum_eddy
