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
  int64_t iDir, iIon, iIon_, iNeutral, iNeutral_;

  arma_cube rho_n(nX, nY, nZ);
  arma_cube rho_i(nX, nY, nZ);
  arma_cube rho_sum(nX, nY, nZ);
  
  //  energy is the total energy transfered from ions to neutrals
  arma_cube energy(nX, nY, nZ);
  // beta is the sum of the collision frequencies * mass density of ions
  arma_cube beta(nX, nY, nZ);
  // velocity difference between ions and neutrals
  arma_cube vDiff(nX, nY, nZ);
  // momentum so we can divide by the mass density later
  std::vector<arma_cube> momentum;
  momentum =  make_cube_vector(nX, nY, nZ, 3);

  beta.zeros();

  // If we are using the bulk (horizontal, primarily) neutral winds
  // then approximate some of the collisional quantities

  neutrals.heating_ion_friction_scgc.zeros();
  neutrals.heating_ion_heat_transfer_scgc.zeros();

  if (input.get_advection_neutrals_bulkwinds()) {
    for (iIon = 0; iIon < ions.nSpeciesAdvect; iIon++) {
      iIon_ = ions.species_to_advect[iIon];
      Ions::species_chars & advected_ion = ions.species[iIon_];
      rho_i = advected_ion.mass * advected_ion.density_scgc;
      for (iNeutral = 0; iNeutral < neutrals.nSpeciesAdvect; iNeutral++) {
        iNeutral_ = neutrals.species_to_advect[iNeutral];
        beta = beta + rho_i % advected_ion.nu_ion_neutral_vcgc[iNeutral_];
      }
    }

    // Now use the bulk quantities for the collisions
    // (beta is included in the last step)
    // heat transfer between ions and neutrals:
    neutrals.heating_ion_heat_transfer_scgc = 3 * cKB / ions.mean_major_mass_scgc %
      (ions.temperature_scgc - neutrals.temperature_scgc);
    for (iDir = 0; iDir < 3; iDir++) {
      // need the velocity difference for momentum and energy eqns:
      vDiff = (ions.velocity_vcgc[iDir] - neutrals.velocity_vcgc[iDir]);
      // ion - neutral drag (acceleration):
      neutrals.acc_ion_collisions[iDir] = 
        beta / neutrals.rho_scgc % vDiff;
      // Frictional heating between ions and neutrals:
      neutrals.heating_ion_friction_scgc = neutrals.heating_ion_friction_scgc + vDiff % vDiff;
    }
    // multiply by collision frequencies and convert
    // energy change to temperature change:
    neutrals.heating_ion_friction_scgc = 
      beta % neutrals.heating_ion_friction_scgc / (2 * neutrals.rho_scgc % neutrals.Cv_scgc);
    neutrals.heating_ion_heat_transfer_scgc = 
      beta % neutrals.heating_ion_friction_scgc / (2 * neutrals.rho_scgc % neutrals.Cv_scgc);
  } else {
    energy.zeros();

    // Calculate acceleration due to ion drag. Based on Formula 4.124b in Ionospheres text.
    for (iNeutral = 0; iNeutral < neutrals.nSpeciesAdvect; iNeutral++) {
      Neutrals::species_chars & advected_neutral =
        neutrals.species[neutrals.species_to_advect[iNeutral]];
      rho_n = advected_neutral.mass * advected_neutral.density_scgc;

      for (iDir = 0; iDir < 3; iDir++)
        momentum[iDir].zeros();

      for (iIon = 0; iIon < ions.nSpeciesAdvect; iIon++) {
        Ions::species_chars & advected_ion = ions.species[ions.species_to_advect[iIon]];
        rho_i = advected_ion.mass * advected_ion.density_scgc;
        beta = rho_i % advected_ion.nu_ion_neutral_vcgc[iNeutral];
        precision_t one_over_masses = 1.0 / (advected_ion.mass + advected_neutral.mass);

        // B = rho_i * Nu_in
        // Acc (for each species) = sum_over_ions(B * (Vi - Vn)) / rho
        //     Momentum = sum(B * (Vi - Vn))
        // Energy = sum_neutrals(sum__ions(B/(Mi + Mn) * (Ti - Tn) + Mi * (Vi-Vn)^2))

        neutrals.heating_ion_heat_transfer_scgc = 
          neutrals.heating_ion_heat_transfer_scgc + 
          3 * cKB * one_over_masses * 
          (ions.temperature_scgc - neutrals.temperature_scgc);

        for (iDir = 0; iDir < 3; iDir++) {
          vDiff = (advected_ion.par_velocity_vcgc[iDir] +
                   advected_ion.perp_velocity_vcgc[iDir] -
                   advected_neutral.velocity_vcgc[iDir]);
          neutrals.heating_ion_friction_scgc = 
            neutrals.heating_ion_friction_scgc + 
            (advected_ion.mass * one_over_masses) * vDiff % vDiff;
          momentum[iDir] = momentum[iDir] + beta % vDiff;

        } // for each ion
        // 
        neutrals.heating_ion_friction_scgc = 
          neutrals.heating_ion_friction_scgc % beta;
        neutrals.heating_ion_heat_transfer_scgc = 
          neutrals.heating_ion_heat_transfer_scgc % beta;
      } // for each ion
      // Divide by the mass density to get the acceleration
      for (iDir = 0; iDir < 3; iDir++)
        advected_neutral.acc_ion_drag[iDir] = momentum[iDir]/rho_n;
    } // for each neutral
    // Convert from energy into K/s:
    neutrals.heating_ion_friction_scgc = 
        neutrals.heating_ion_friction_scgc / (neutrals.rho_scgc % neutrals.Cv_scgc);
    neutrals.heating_ion_heat_transfer_scgc = 
        neutrals.heating_ion_heat_transfer_scgc / (neutrals.rho_scgc % neutrals.Cv_scgc);    
  } // bulk neutral winds

  report.exit(function);
  return;
} // calc_ion_collisions
