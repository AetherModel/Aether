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

  // If we only consider the bulk winds in the horizontal direction:
  if (input.get_advection_neutrals_bulkwinds()) {
    // Calculate Coriolis:
    acc_coriolis = coriolis(velocity_vcgc, planet.get_omega(), grid.geoLat_scgc);
    // Add Velocity sources to bulk winds:
    for (int iDir = 0; iDir < 3; iDir++) {
      velocity_vcgc[iDir] = velocity_vcgc[iDir] + dt * (
        grid.cent_acc_vcgc[iDir] +
        acc_coriolis[iDir] + 
        acc_ion_collisions[iDir]);
    }
  } else {
    for (int64_t iSpec = 0; iSpec < nSpeciesAdvect; iSpec++) {
      // Pick out the advected neutral species:
      species_chars & advected_neutral = species[species_to_advect[iSpec]];
      // Calculate Coriolis:
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

  report.exit(function);
  return;
}
