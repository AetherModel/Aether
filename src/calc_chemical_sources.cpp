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
  // Then zero it out:
  change3d.zeros();

  // likewise for chemical heating:
  arma_cube chemical_heating = neutrals.temperature_scgc;
  chemical_heating.zeros();

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

    // check for type of temperature dependence and adjust
    if (reactions[iReaction].type > 0) {
      // use Ti by default
      arma_cube temp = ions.temperature_scgc;
      std::string denom = reactions[iReaction].denominator;

      if (denom == "Te")
        temp = ions.electron_temperature_scgc;

      else if (denom == "Tn")
        temp = neutrals.temperature_scgc;

      if (reactions[iReaction].numerator &&
          reactions[iReaction].type == 1) {
        change3d =
          change3d %
          pow(reactions[iReaction].numerator / temp,
              reactions[iReaction].exponent);
      } else if (reactions[iReaction].numerator &&
                 reactions[iReaction].type == 2) {
        change3d =
          change3d %
          temp %
          exp(reactions[iReaction].numerator / temp);
      } else if (reactions[iReaction].numerator &&
                 reactions[iReaction].type == 3) {

        temp = temp + 0.33 *
               pow(ions.efield_vcgc[0], 2) %
               pow(ions.efield_vcgc[1], 2) %
               pow(ions.efield_vcgc[2], 2); //.33 * E'^2

        precision_t coeff_a, coeff_b, coeff_c;

        if (denom == "12.9a") {
          coeff_a = 1.533  * 0.000000000001;   //10^-12
          coeff_b = -5.92  * 0.0000000000001;  //10^-13
          coeff_c = 8.60   * 0.00000000000001; //10^-14
        } else if (denom == "12.9b") {
          coeff_a = 2.73   * 0.000000000001;   //10^-12
          coeff_b = -1.155 * 0.000000000001;   //10^-12
          coeff_c = 1.483  * 0.0000000000001;  //10^-13
        }

        change3d.fill(coeff_a);
        change3d += coeff_b * temp / 300;
        change3d += coeff_c * pow(temp / 300, 2);

      }
    }

    // if temperature dependence is piecewise, only operate on cells within range
    if (reactions[iReaction].min || reactions[iReaction].max) {
      arma_cube temp = ions.temperature_scgc;
      std::string piecewiseTemp = reactions[iReaction].piecewiseVar;

      if (piecewiseTemp == "Te")
        temp = ions.electron_temperature_scgc;

      else if (piecewiseTemp == "Tn")
        temp = neutrals.temperature_scgc;

      change3d = change3d % (change3d > reactions[iReaction].min);

      if (reactions[iReaction].max > 0)
        change3d = change3d % (change3d <= reactions[iReaction].max);
    }

    for (iLoss = 0; iLoss < reactions[iReaction].nLosses; iLoss++) {
      IsNeutral = reactions[iReaction].losses_IsNeutral[iLoss];
      id_ = reactions[iReaction].losses_ids[iLoss];

      if (IsNeutral)
        change3d = change3d % neutrals.species[id_].density_scgc;
      else
        change3d = change3d % ions.species[id_].density_scgc;
    }

    // calculate heat change
    chemical_heating += change3d * reactions[iReaction].energy;

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

  // convert chemical heating to temperature change in K / s
  neutrals.heating_chemical_scgc =
    chemical_heating * cE /
    neutrals.Cv_scgc / neutrals.rho_scgc;

  report.exit(function);
}
