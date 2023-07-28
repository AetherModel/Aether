// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "../include/aether.h"

void calc_ion_collisions(Neutrals &neutrals,
                         Ions &ions) {

  std::string function = "calc_ion_collisions";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t nX = neutrals.density_scgc.n_rows;
  int64_t nY = neutrals.density_scgc.n_cols;
  int64_t nZ = neutrals.density_scgc.n_slices;
  int64_t nSpecies = neutrals.nSpecies, iSpecies;
  int64_t iDir, iIon, iNeutral;

  arma_cube rho_n(nX, nY, nZ);
  arma_cube rho_i(nX, nY, nZ);
  arma_cube rho_sum(nX, nY, nZ);

  // Calculate acceleration due to ion drag. Based on Formula 4.124b
  // in Ionospheres text.
  for (iNeutral = 0; iNeutral < nSpecies; iNeutral++) {
    rho_n =
      neutrals.species[iNeutral].mass *
      neutrals.species[iNeutral].density_scgc;

    for (iDir = 0; iDir < 3; iDir++) {
      rho_sum.zeros();

      for (iIon = 0; iIon < ions.nSpecies; iIon++) {
        rho_i = ions.species[iIon].mass * ions.species[iIon].density_scgc;
        rho_sum = rho_sum +
                  rho_i % ions.species[iIon].nu_ion_neutral_vcgc[iNeutral] %
                  (ions.species[iIon].par_velocity_vcgc[iDir] +
                   ions.species[iIon].perp_velocity_vcgc[iDir] -
                   neutrals.species[iNeutral].velocity_vcgc[iDir]);
      } // for each ion

      neutrals.species[iNeutral].acc_ion_drag[iDir] = rho_sum / rho_n;
    } // for each direction
  } // for each neutral

  report.exit(function);
  return;
} // calc_ion_collisions
