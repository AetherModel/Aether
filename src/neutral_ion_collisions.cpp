// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "../include/aether.h"

void calc_ion_collisions(Neutrals &neutrals,
			 Ions &ions,
			 Report &report) {
    
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
    
  // Calculate acceleration due to ion drag. Based on Formula 4.124b in Ionospheres text.
  for (iNeutral = 0; iNeutral < neutrals.nSpeciesAdvect; iNeutral++) {
    Neutrals::species_chars & advected_neutral = neutrals.species[neutrals.species_to_advect[iNeutral]];
    rho_n = advected_neutral.mass * advected_neutral.density_scgc;
    for (iDir = 0; iDir < 3; iDir++) {
      rho_sum.zeros();
      for (iIon = 0; iIon < ions.nSpeciesAdvect; iIon++) {
          Ions::species_chars & advected_ion = ions.species[ions.species_to_advect[iIon]];
	rho_i = advected_ion.mass * advected_ion.density_scgc;
	rho_sum = rho_sum +
	  rho_i % advected_ion.nu_ion_neutral_vcgc[iNeutral] %
	  (advected_ion.par_velocity_vcgc[iDir] +
	   advected_ion.perp_velocity_vcgc[iDir] -
	   advected_neutral.velocity_vcgc[iDir]);
      } // for each ion
      advected_neutral.acc_ion_drag[iDir] = rho_sum / rho_n;
    } // for each direction
  } // for each neutral
    
  report.exit(function);
    
} // calc_ion_collisions
