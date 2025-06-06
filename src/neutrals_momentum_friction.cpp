// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "aether.h"

// ---------------------------------------------------------------------
// This calculates the velocities in one cell. The inputs are the
// updated velocities given all of the acceleration terms. This is the right
// side of the equation (source terms * dt + Vold).  The velocities
// of all species are solved for together, since they are interdependent.
// This is done with a matrix inversion, and is the cell-by-cell
// implicit solver.  There is a semi-implicit solver below.
// ---------------------------------------------------------------------

arma_vec Neutrals::calc_friction_one_cell(int64_t iLon, int64_t iLat,
                                          int64_t iAlt,
                                          precision_t dt, arma_vec &vels) {
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
      // It is not clear where the 1000 comes from. Need to figure this out.
      temp_dij =
        advected_neutral.diff0[jSpecies_] *
        pow(temperature_scgc(iLon, iLat, iAlt),
            advected_neutral.diff_exp[jSpecies_]) * 1000.0;
      coefmatrix(iSpecies, jSpecies) =
        ktom * species[jSpecies_].density_scgc(iLon, iLat, iAlt) / temp_dij;
    } // jSpec loop
  } // iSpec loop

  matrix = -1 * coefmatrix * dt;

  // Fill in diagonal of matrix:
  for (iSpecies = 0; iSpecies < nSpeciesAdvect; iSpecies++)
    matrix(iSpecies, iSpecies) = 1 - sum(coefmatrix.row(iSpecies));

  // initialize array of each neutral species' accelerations at (iLon, iLat, iAlt):
  arma_vec accs(nSpeciesAdvect, fill::zeros);
  // Solve system of equations:
  arma_vec new_vels = arma::solve(matrix, vels, solve_opts::fast);

  report.exit(function);
  return new_vels;
}

// ---------------------------------------------------------------------
// This solves for the velocities using a fully implicit solver, where
// all of the velocities are solved for at the same time. This is slower
// since there needs to be a full matrix solve in each cell, but it ties
// the velocities together.
// ---------------------------------------------------------------------

void Neutrals::calc_neutral_friction_implicit(precision_t dt) {

  std::string function = "calc_neutral_friction";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t iAlt, iLat, iLon, iDir, iSpecies, iSpecies_;
  int64_t jSpecies, jSpecies_;

  arma_cube ktom;
  arma_cube coef;
  arma_cube tpower;

  // Initialize all of the accelerations to zero:
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    for (iDir = 0; iDir < 3; iDir++)
      species[iSpecies].acc_neutral_friction[iDir].zeros();
  }

  if (input.get_advection_neutrals_vertical() != "hydro") {

    arma_vec vels(nSpeciesAdvect, fill::zeros);
    arma_vec acc(nSpeciesAdvect, fill::zeros);
    arma_vec new_vels(nSpeciesAdvect, fill::zeros);
    int64_t nXs = temperature_scgc.n_rows;
    int64_t nYs = temperature_scgc.n_cols;
    int64_t nZs = temperature_scgc.n_slices;
    int64_t nGCs = 2;

    // Calculate friction terms for only species that advect.
    //   - If only 1 species is advected, then it will have no friction
    if (nSpeciesAdvect > 1) {
      for (iAlt = nGCs; iAlt < nZs - nGCs; iAlt++) {
        for (iLat = nGCs; iLat < nYs - nGCs; iLat++) {
          for (iLon = nGCs; iLon < nXs - nGCs; iLon++) {
            // Only worry about the vertical direction for now:
            iDir = 2;
            //for (iDir = 0; iDir < 3; iDir++) {
            vels.zeros();

            // The velocities are just after the vertical solver, so the velocities are
            // the source terms for the friction solver.
            for (iSpecies = 0; iSpecies < nSpeciesAdvect; iSpecies++) {
              iSpecies_ = species_to_advect[iSpecies];
              vels(iSpecies) =
                species[iSpecies_].newVelocity_vcgc[iDir](iLon, iLat, iAlt);
            }

            // = neutral_friction_one_cell(iLon, iLat, iAlt, vels);
            new_vels = calc_friction_one_cell(iLon, iLat, iAlt, dt, vels);

            for (iSpecies = 0; iSpecies < nSpeciesAdvect; iSpecies++) {
              iSpecies_ = species_to_advect[iSpecies];
              species[iSpecies_].newVelocity_vcgc[iDir](iLon, iLat,
                                                        iAlt) = new_vels(iSpecies);
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


// ---------------------------------------------------------------------
// If we are solving this system semi-implicitly, we don't need to do a
// matrix solve in each cell, so we can calculate all of the variables
// before doing the time-stepping in the vertical solver.
// semi-implicit means using the current velocity on both the left and
// right sides of the equation, but using the old velocities for the
// other species.  This means that the velocities are not consistent
// in any time-step unless they are close to steady-state, but they will
// move together.
// ---------------------------------------------------------------------

void Neutrals::calc_neutral_friction_coefs() {

  std::string function = "calc_neutral_friction_coefs";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t iAlt, iLat, iLon, iDir, iSpecies, iSpecies_;
  int64_t jSpecies, jSpecies_;

  arma_cube ktom;
  arma_cube coef;
  arma_cube tpower;

  // Initialize all of the accelerations to zero:
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    species[iSpecies_].neutral_friction_coef.zeros();

    for (iDir = 0; iDir < 3; iDir++)
      species[iSpecies].acc_neutral_friction[iDir].zeros();
  }

  if (input.get_advection_neutrals_vertical() != "hydro") {

    for (iSpecies = 0; iSpecies < nSpeciesAdvect; iSpecies++) {
      iSpecies_ = species_to_advect[iSpecies];
      // ktom = boltzmann's constant * temperature / mass
      ktom =
        cKB *
        temperature_scgc /
        species[iSpecies_].mass;

      for (jSpecies = 0; jSpecies < nSpeciesAdvect; jSpecies++) {
        jSpecies_ = species_to_advect[jSpecies];

        if (iSpecies_ == jSpecies_)
          continue;

        tpower = pow(temperature_scgc, species[iSpecies_].diff_exp[jSpecies_]);
        // NEED TO REMOVE /100!!!
        coef = ktom % species[jSpecies_].density_scgc /
               (species[iSpecies_].diff0[jSpecies_] * tpower) / 100.0;
        species[iSpecies_].neutral_friction_coef =
          species[iSpecies_].neutral_friction_coef + coef;
        iDir = 2;
        species[iSpecies_].acc_neutral_friction[iDir] =
          species[iSpecies_].acc_neutral_friction[iDir] +
          coef % species[iSpecies_].velocity_vcgc[iDir];
      }
    }
  } // if !hydro

  report.exit(function);
  return;
} //calc_neutral_friction
