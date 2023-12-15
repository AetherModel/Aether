// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "aether.h"

// ---------------------------------------------------------------------
// This calculates the acceleration due to neutral-neutral
// friction between species.
// dt was removed, since we can assume that it is 1.0 and then
// multiply by dt later.
// ---------------------------------------------------------------------

arma_vec Neutrals::calc_friction_one_cell(int64_t iLon, int64_t iLat, int64_t iAlt,
                                   arma_vec &vels) {
  std::string function = "neutral_friction_one_cell";
  static int iFunction = -1;
  report.enter(function, iFunction);
  precision_t ktom, temp_dij;
  int64_t iSpecies, jSpecies, iSpecies_, jSpecies_;

  static arma_mat matrix(nSpeciesAdvect, nSpeciesAdvect, fill::zeros);
  static arma_mat coefmatrix(nSpeciesAdvect, nSpeciesAdvect, fill::zeros);

  for (iSpecies = 0; iSpecies < nSpeciesAdvect; iSpecies++) {
    iSpecies_ = species_to_advect[iSpecies];

    Neutrals::species_chars & advected_neutral = species[iSpecies_];

    // ktom = boltzmann's constant * temperature / mass
    ktom =
      cKB *
      temperature_scgc(iLon, iLat, iAlt) /
      advected_neutral.mass;

    for (jSpecies = 0; jSpecies < nSpeciesAdvect; jSpecies++) {
      jSpecies_ = species_to_advect[jSpecies];

      if (iSpecies_ == jSpecies_)
        continue;

      // temp_dij holds the dij binary coefficients based upon the formulation by Banks and Kokarts.
      // These coefficients demand that:
      // (1) density be in cm^-3 (hence the 1.0e-06) factor below
      // (2) Additionally, the Dij's are in cm^2/s, thus the 1.0e-04 factor
      temp_dij =
        1e-4 * advected_neutral.diff0[jSpecies_] *
        pow(temperature_scgc(iLon, iLat, iAlt),
            advected_neutral.diff_exp[jSpecies_]) /
        (density_scgc(iLon, iLat, iAlt) * 1e-6);

      coefmatrix(iSpecies, jSpecies) =
        ktom * advected_neutral.density_scgc(iLon, iLat, iAlt) /
        (temp_dij * density_scgc(iLon, iLat, iAlt));
    } // jSpec loop
  } // iSpec loop

  matrix = -1 * coefmatrix;

  // Fill in diagonal of matrix:
  for (iSpecies = 0; iSpecies < nSpeciesAdvect; iSpecies++)
    matrix(iSpecies, iSpecies) = 1 + sum(coefmatrix.row(iSpecies));

  // initialize array of each neutral species' accelerations at (iLon, iLat, iAlt):
  arma_vec accs(nSpeciesAdvect, fill::zeros);

  // Solve system of equations:
  arma_vec new_vels = arma::solve(matrix, vels, solve_opts::fast);

  /*
  // put the new values into the velocity cubes:
  for (iSpecies = 0; iSpecies < nSpeciesAdvect; iSpecies++)
    accs(iSpecies) = new_vels(iSpecies) - vels(iSpecies);
  */
  report.exit(function);
  return new_vels;
}

// ---------------------------------------------------------------------
//
// ---------------------------------------------------------------------

void Neutrals::calc_neutral_friction() {

  std::string function = "calc_neutral_friction";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t iAlt, iLat, iLon, iDir, iSpecies, iSpecies_;

  // Initialize all of the accelerations to zero:
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    for (iDir = 0; iDir < 3; iDir++)
      species[iSpecies].acc_neutral_friction[iDir].zeros();

  if (input.get_advection_neutrals_vertical() != "hydro") {

    arma_vec vels(nSpeciesAdvect, fill::zeros);
    arma_vec acc(nSpeciesAdvect, fill::zeros);
    arma_vec new_vels;
    int64_t nXs = temperature_scgc.n_rows;
    int64_t nYs = temperature_scgc.n_cols;
    int64_t nZs = temperature_scgc.n_slices;

    // Calculate friction terms for only species that advect.
    //   - If only 1 species is advected, then it will have no friction
    if (nSpeciesAdvect > 1) {
      for (iAlt = 0; iAlt < nZs; iAlt++) {
        for (iLat = 0; iLat < nYs; iLat++) {
          for (iLon = 0; iLon < nXs; iLon++) {
            iDir = 2;
            //for (iDir = 0; iDir < 3; iDir++) {
              vels.zeros();

              //Put the old velocities into vels:
              for (iSpecies = 0; iSpecies < nSpeciesAdvect; iSpecies++) {
                iSpecies_ = species_to_advect[iSpecies];
                vels(iSpecies) =
                  species[iSpecies_].velocity_vcgc[iDir](iLon, iLat, iAlt);
              }

              //acc = neutral_friction_one_cell(iLon, iLat, iAlt, vels);
              new_vels = calc_friction_one_cell(iLon, iLat, iAlt, vels);

              for (iSpecies = 0; iSpecies < nSpeciesAdvect; iSpecies++) {
                iSpecies_ = species_to_advect[iSpecies];
                species[iSpecies_].velocity_vcgc[iDir](iLon, iLat, iAlt) = new_vels(iSpecies);
                //species[iSpecies_].acc_neutral_friction[iDir](iLon, iLat, iAlt) =
                //  acc(iSpecies);
              } // iSpeciesAdvect
            //} // for direction
          } // for long
        } // for lat
      } // for alt
    } // if nSpecies > 1
  } // if !hydro

  report.exit(function);
  return;
} //calc_neutral_friction
