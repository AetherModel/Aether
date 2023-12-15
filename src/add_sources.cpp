// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>
#include <string>

#include "aether.h"

// -----------------------------------------------------------------------------
// Adds all of the sources to the states. Needs time to get dt.
// -----------------------------------------------------------------------------

void Neutrals::add_sources(Times time, Planets planet, Grid grid) {

  std::string function = "add_sources";
  static int iFunction = -1;
  report.enter(function, iFunction);

  precision_t dt = time.get_dt();

  heating_sources_total = heating_euv_scgc
                        + heating_chemical_scgc
                        + heating_ion_collisions_scgc
                        - O_cool_scgc
                        - NO_cool_scgc;
  
  // Solve the laplace equations using the source terms,
  // updating the neutral temperature:
  update_temperature(grid, time);

  std::vector<arma_cube> acc_coriolis;
  acc_coriolis = make_cube_vector(grid.get_nX(), grid.get_nY(), grid.get_nZ(), 3);

  int64_t iDir, iSpec, iSpecies;
  double tSim = time.get_simulation_time();
  precision_t ramp = tSim / 3600.0;
  if (ramp > 1.0) ramp = 1.0;
  
  // Vertical winds use species winds:
  for (iSpec = 0; iSpec < nSpeciesAdvect; iSpec++) {
    // Pick out the advected neutral species:
    species_chars & advected_neutral = species[species_to_advect[iSpec]];
    // Calculate Coriolis:
    if (input.get_use_coriolis()) 
      acc_coriolis = coriolis(advected_neutral.velocity_vcgc, 
			      planet.get_omega(), 
			      grid.geoLat_scgc);

    iDir = 2;
    // update velocities based on acceleration:
    // reduce neutral friction until solver is added
    advected_neutral.velocity_vcgc[iDir] =
      advected_neutral.velocity_vcgc[iDir] +
      dt * (ramp * grid.cent_acc_vcgc[iDir] + 
	    ramp * acc_coriolis[iDir] + 
	    advected_neutral.acc_neutral_friction[iDir] / 4.0 +
	    advected_neutral.acc_ion_drag[iDir] +
	    advected_neutral.acc_eddy);
  }

  calc_mass_density();
  // Calculate bulk vertical winds:
  velocity_vcgc[2].zeros();
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    if (species[iSpecies].DoAdvect) {
      velocity_vcgc[2] = velocity_vcgc[2] + 
        species[iSpecies].mass * species[iSpecies].density_scgc % 
        species[iSpecies].velocity_vcgc[2] / rho_scgc;
    }
  
  
  // Horizontal winds use bulk winds:
  if (input.get_use_coriolis()) 
    acc_coriolis = coriolis(velocity_vcgc, planet.get_omega(), grid.geoLat_scgc);
  // Add Velocity sources to bulk winds:
  for (iDir = 0; iDir < 2; iDir++) {
    velocity_vcgc[iDir] =
      velocity_vcgc[iDir] + dt * (
				  ramp * grid.cent_acc_vcgc[iDir] +
				  ramp * acc_coriolis[iDir] + 
				  acc_ion_collisions[iDir]);
    acc_sources_total[iDir].zeros();
  }
  // Apply Viscosity:
  update_horizontal_velocity(grid, time);

  // Assign bulk horizontal velocity to all species:
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    for (iDir = 0; iDir < 2; iDir++)
      species[iSpecies].velocity_vcgc[iDir] = velocity_vcgc[iDir];

  /*  
  // If we only consider the bulk winds in the horizontal direction:
  if (input.get_advection_neutrals_bulkwinds()) {
    // Calculate Coriolis:
    if (input.get_use_coriolis()) 
      acc_coriolis = coriolis(velocity_vcgc, planet.get_omega(), grid.geoLat_scgc);
    // Add Velocity sources to bulk winds:
    for (int iDir = 0; iDir < 3; iDir++) {
      velocity_vcgc[iDir] = velocity_vcgc[iDir] + dt * (
        grid.cent_acc_vcgc[iDir] +
        acc_coriolis[iDir] + 
        acc_ion_collisions[iDir]);
      acc_sources_total[iDir].zeros();
    }
    // Apply Viscosity:
    update_horizontal_velocity(grid, time);
  } else {
    for (int64_t iSpec = 0; iSpec < nSpeciesAdvect; iSpec++) {
      // Pick out the advected neutral species:
      species_chars & advected_neutral = species[species_to_advect[iSpec]];
      // Calculate Coriolis:
      if (input.get_use_coriolis()) 
        acc_coriolis = coriolis(advected_neutral.velocity_vcgc, 
                                planet.get_omega(), 
                                grid.geoLat_scgc);

      for (int iDir = 0; iDir < 2; iDir++) {
        // update velocities based on acceleration:
        // reduce neutral friction until solver is added
        advected_neutral.velocity_vcgc[iDir] =
          advected_neutral.velocity_vcgc[iDir] +
          dt * (grid.cent_acc_vcgc[iDir] + 
                acc_coriolis[iDir] + 
                advected_neutral.acc_neutral_friction[iDir] / 4.0 +
                advected_neutral.acc_ion_drag[iDir]);
        // eddy acceleration is only in the vertical direction:
        if (iDir == 2)
          advected_neutral.velocity_vcgc[iDir] =
            advected_neutral.velocity_vcgc[iDir] +
            dt * advected_neutral.acc_eddy;
      }
    }
    calc_bulk_velocity();
  }
  assign_bulk_velocity();
  */
  report.exit(function);
  return;
}
