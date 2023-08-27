// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>

#include "aether.h"

// -----------------------------------------------------------------------------
// loops through all of the chemical reactions doing 4 things:
//   1. Determine change (per unit time) of particles of loss
//   2. Add this to sources (keeping track of ions v. neutrals)
//   3. Add this to losses (keeping track of ions v. neutrals)
//   4. Figures out the chemical heating through exothermic reactions
//
// 2022 - A. Ridley
// 2023/03 - M. Rinaldi
// -----------------------------------------------------------------------------

void Chemistry::calc_chemical_sources(Neutrals &neutrals,
                                      Ions &ions) {

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

    // Grab reaction rate. For temperature dependent rates, this is
    // the multiplicative factor in front of the equation:
    rate = reactions[iReaction].rate;

    // First calculate the amount of change:
    //    Change is calculated as
    //    reaction rate * loss den 1 * loss den 2 (* loss den 3 if needed)

    change3d.fill(rate);

    // Check for type of temperature dependence and calculate
    if (reactions[iReaction].type > 0) {
      // Determined which temperature to use in equation:
      // use Tn by default
      arma_cube temp = neutrals.temperature_scgc;
      std::string denom = reactions[iReaction].denominator;

      if (denom == "Te")
        temp = ions.electron_temperature_scgc;
      else if (denom == "Ti")
        temp = ions.temperature_scgc;

      // Calculate reaction rate:
      if (reactions[iReaction].numerator &&
          reactions[iReaction].type == 1) {
        // Form is RR = R * (num / Temp) ^ exp
        change3d =
          change3d %
          pow(reactions[iReaction].numerator / temp,
              reactions[iReaction].exponent);
      } else if (reactions[iReaction].numerator &&
                 reactions[iReaction].type == 2) {
        // Form is RR = R * exp(num / Temp)

        change3d =
          change3d %
          temp %
          exp(reactions[iReaction].numerator / temp);
      } else if (reactions[iReaction].numerator &&
                 reactions[iReaction].type == 3) {
        // This is a placeholder for more complicated reaction rates,
        // such as the charge exchange for O+ + N2 at Earth. Specifically,
        // this is what is outlined in Schunk and Nagy:
        temp = temp + 0.33 * (
                 pow(ions.efield_vcgc[0], 2) +
                 pow(ions.efield_vcgc[1], 2) +
                 pow(ions.efield_vcgc[2], 2)); //.33 * E'^2

        precision_t coeff_a, coeff_b, coeff_c;

        if (denom == "12.9a") {
          coeff_a = 1.533e-12;   //10^-12
          coeff_b = -5.92e-13 ;  //10^-13
          coeff_c = 8.60e-14; //10^-14
        } else if (denom == "12.9b") {
          coeff_a = 2.73e-12;   //10^-12
          coeff_b = -1.155e-12;   //10^-12
          coeff_c = 1.483e-13;  //10^-13
        }

        change3d.fill(coeff_a);
        change3d += coeff_b * temp / 300;
        change3d += coeff_c * pow(temp / 300, 2);

      }
    }

    // If temperature dependence is piecewise, only operate on cells
    // within temperature range:
    if (reactions[iReaction].min || reactions[iReaction].max) {
      // Figure out which temperature is the limiter.  Default to ions:
      arma_cube temp = ions.temperature_scgc;
      std::string piecewiseTemp = reactions[iReaction].piecewiseVar;

      if (piecewiseTemp == "Te")
        temp = ions.electron_temperature_scgc;
      else if (piecewiseTemp == "Tn")
        temp = neutrals.temperature_scgc;

      // Limit the reagion to where the temperautre is in the range:
      change3d = change3d % (change3d > reactions[iReaction].min);

      if (reactions[iReaction].max > 0)
        change3d = change3d % (change3d <= reactions[iReaction].max);
    }

    // Now that the reaction rate is calculated, multiply by the
    // densities on the left side of the equation (loss terms):
    for (iLoss = 0; iLoss < reactions[iReaction].nLosses; iLoss++) {
      // Determine if constituent is a neutral, and grab it's id:
      IsNeutral = reactions[iReaction].losses_IsNeutral[iLoss];
      id_ = reactions[iReaction].losses_ids[iLoss];

      if (IsNeutral)
        change3d = change3d % neutrals.species[id_].density_scgc;
      else
        change3d = change3d % ions.species[id_].density_scgc;
    }

    // calculate heat change
    chemical_heating += change3d * reactions[iReaction].energy;

    // Now that full loss term is calculated, we can then add this
    // value to the losses:
    for (iLoss = 0; iLoss < reactions[iReaction].nLosses; iLoss++) {
      // Once again, figure out if it is a neutral, determine id, and add:
      IsNeutral = reactions[iReaction].losses_IsNeutral[iLoss];
      id_ = reactions[iReaction].losses_ids[iLoss];

      if (IsNeutral)
        neutrals.species[id_].losses_scgc =
          neutrals.species[id_].losses_scgc + change3d;
      else
        ions.species[id_].losses_scgc =
          ions.species[id_].losses_scgc + change3d;
    }

    // Then add this to the sources:
    for (iSource = 0; iSource < reactions[iReaction].nSources; iSource++) {
      // Once again, figure out if it is a neutral, determine id, and add:
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
