// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>
#include <string>

#include "aether.h"

// -----------------------------------------------------------------------------
// Adds all of the sources to the states. Needs time to get dt.
// -----------------------------------------------------------------------------

void Neutrals::add_sources(Times time) {

  std::string function = "add_sources";
  static int iFunction = -1;
  report.enter(function, iFunction);

  precision_t dt = time.get_dt();

  temperature_scgc = temperature_scgc +
                     dt * (heating_euv_scgc
                           + heating_chemical_scgc
                           + conduction_scgc
                           - O_cool_scgc
                           - NO_cool_scgc);

  for (int64_t iSpec = 0; iSpec < nSpecies; iSpec++) {
    for (int iDir = 0; iDir < 3; iDir++) {
      // update velocities based on acceleration:
      // reduce neutral friction until solver is added
      species[iSpec].velocity_vcgc[iDir] =
        species[iSpec].velocity_vcgc[iDir] +
        dt * (species[iSpec].acc_neutral_friction[iDir] / 4.0 +
              species[iSpec].acc_ion_drag[iDir]);

      // eddy acceleration is only in the vertical direction:
      if (iDir == 2)
        species[iSpec].velocity_vcgc[iDir] =
          species[iSpec].velocity_vcgc[iDir] +
          dt * species[iSpec].acc_eddy;
    }
  }

  calc_bulk_velocity();

  report.exit(function);
  return;
}
