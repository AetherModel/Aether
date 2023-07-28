// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "aether.h"

// ---------------------------------------------------------------------
// This calculates the acceleration due to neutral-neutral
// friction between species.
// dt was removed, since we can assume that it is 1.0 and then
// multiply by dt later.
// ---------------------------------------------------------------------

arma_vec neutral_friction_one_cell(int64_t iLon, int64_t iLat, int64_t iAlt,
                                   arma_vec &vels,
                                   Neutrals &neutrals) {

  std::string function = "neutral_friction_one_cell";
  static int iFunction = -1;
  report.enter(function, iFunction);
  precision_t ktom, temp_dij;
  int64_t iSpecies, jSpecies;

  static arma_mat matrix(neutrals.nSpecies, neutrals.nSpecies, fill::zeros);
  static arma_mat coefmatrix(neutrals.nSpecies, neutrals.nSpecies, fill::zeros);

  for (iSpecies = 0; iSpecies < neutrals.nSpecies; iSpecies++) {

    // ktom = boltzmann's constant * temperature / mass
    ktom =
      cKB *
      neutrals.temperature_scgc(iLon, iLat, iAlt) /
      neutrals.species[iSpecies].mass;

    for (jSpecies = 0; jSpecies < neutrals.nSpecies; jSpecies++) {
      if (iSpecies == jSpecies)
        continue;

      // temp_dij holds the dij binary coefficients based upon the formulation by Banks and Kokarts.
      // These coefficients demand that:
      // (1) density be in cm^-3 (hence the 1.0e-06) factor below
      // (2) Additionally, the Dij's are in cm^2/s, thus the 1.0e-04 factor
      temp_dij =
        1e-4 * neutrals.species[iSpecies].diff0[jSpecies] *
        pow(neutrals.temperature_scgc(iLon, iLat, iAlt),
            neutrals.species[iSpecies].diff_exp[jSpecies]) /
        (neutrals.density_scgc(iLon, iLat, iAlt) * 1e-6);

      coefmatrix(iSpecies, jSpecies) =
        ktom * neutrals.species[jSpecies].density_scgc(iLon, iLat, iAlt) /
        (temp_dij * neutrals.density_scgc(iLon, iLat, iAlt));
    } // jSpec loop
  } // iSpec loop

  matrix = -1 * coefmatrix;

  // Fill in diagonal of matrix:
  for (iSpecies = 0; iSpecies < neutrals.nSpecies; iSpecies++)
    matrix(iSpecies, iSpecies) = 1 + sum(coefmatrix.row(iSpecies));

  // initialize array of each neutral species' accelerations at (iLon, iLat, iAlt):
  arma_vec accs(neutrals.nSpecies, fill::zeros);

  // Solve system of equations:
  arma_vec new_vels = arma::solve(matrix, vels, solve_opts::fast);

  // put the new values into the velocity cubes:
  for (iSpecies = 0; iSpecies < neutrals.nSpecies; iSpecies++)
    accs(iSpecies) = new_vels(iSpecies) - vels(iSpecies);

  report.exit(function);
  return accs;
}

// ---------------------------------------------------------------------
//
// ---------------------------------------------------------------------

void calc_neutral_friction(Neutrals &neutrals) {

  std::string function = "calc_neutral_friction";
  static int iFunction = -1;
  report.enter(function, iFunction);

  arma_vec vels(neutrals.nSpecies, fill::zeros);
  arma_vec acc(neutrals.nSpecies, fill::zeros);
  int64_t iAlt, iLat, iLon, iDir, iSpecies;
  int64_t nXs = neutrals.temperature_scgc.n_rows;
  int64_t nYs = neutrals.temperature_scgc.n_cols;
  int64_t nZs = neutrals.temperature_scgc.n_slices;

  for (iAlt = 0; iAlt < nZs; iAlt++) {
    for (iLat = 0; iLat < nYs; iLat++) {
      for (iLon = 0; iLon < nXs; iLon++) {
        for (iDir = 0; iDir < 3; iDir++) {
          vels.zeros();

          //Put the old velocities into vels:
          for (iSpecies = 0; iSpecies < neutrals.nSpecies; iSpecies++)
            vels(iSpecies) = neutrals.species[iSpecies].velocity_vcgc[iDir](iLon, iLat,
                             iAlt);

          acc = neutral_friction_one_cell(iLon, iLat, iAlt, vels, neutrals);

          for (iSpecies = 0; iSpecies < neutrals.nSpecies; iSpecies++)
            neutrals.species[iSpecies].acc_neutral_friction[iDir](iLon, iLat,
                                                                  iAlt) = acc(iSpecies);
        } // for direction
      } // for long
    } // for lat
  } // for alt

  report.exit(function);
  return;
} //calc_neutral_friction
