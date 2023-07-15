// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "aether.h"

void Neutrals::vertical_momentum_eddy(Grid &gGrid, Inputs inputs, Report &report) {

  std::string function = "Neutrals::vertical_momentum_eddy";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t iSpecies;

  if (inputs.get_use_eddy_momentum()) {

    arma_cube log_cons;
    arma_cube grad_cons;

    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      // Take the natural log of each concentration cube:
      log_cons = log(species[iSpecies].concentration_scgc);

      // calculate gradient:
      grad_cons = calc_gradient_alt(log_cons, gGrid);
      species[iSpecies].acc_eddy = - grad_cons % kappa_eddy_scgc;
    }
  } else {
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      species[iSpecies].acc_eddy.zeros();
  }

  report.exit(function);
  return;
} // neutral_friction_momentum_eddy
